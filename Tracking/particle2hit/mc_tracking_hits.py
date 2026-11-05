"""Standalone MCParticle -> simulated tracking-hit tables (NumPy/Awkward/Uproot/pandas).

read_mc_tracking_hits() returns a dict with particles, hits and surfaces DataFrames
and metadata. For large inputs, use iter_mc_tracking_hits() and reduce each block.
Particles with no hit are retained. No assumption of one particle/event is made.
Only direct hit.particle references are followed; secondary ancestry is not added.

Join particles/hits on KEY. A detector-layer key is (hit_collection, system,
layer), NOT layer alone. It identifies a detector layer, not an ACTS surface;
decoded module/sensor fields and the full CellID retain finer segmentation.
If a readout has no system/layer field, its table value is NA; use its other
decoded fields (e.g. cell_sector/cell_module) instead of inventing a layer ID.

CellID layouts are read from PODIO frame metadata. If absent, pass encodings
from the matching geometry: {readout_name: DD4hep_bitfield_description}.
load_readout_encodings(xml_files) reads these descriptions without DD4hep.
No project modules, PODIO, ROOT or geometry libraries are imported here.
"""
from pathlib import Path
import xml.etree.ElementTree as ET

import awkward as ak
import numpy as np
import pandas as pd
import uproot

KEY = ['file_id', 'event_entry', 'mc_index']
SURFACE_KEY = ['hit_collection', 'system', 'layer']


def load_readout_encodings(xml_files):
    """Read <readout><id> layouts from the actual XMLs used by the simulation.

    Pass the included detector XMLs, not all alternative geometry variants.
    This reads descriptions only; it does not expand XML includes or expressions.
    """
    layouts = {}
    for filename in xml_files:
        for readout in ET.parse(filename).getroot().iter('readout'):
            description = readout.findtext('id')
            if description is None:
                continue
            name, description = readout.attrib['name'], description.strip()
            if name in layouts and layouts[name] != description:
                raise ValueError(f'Conflicting readout layouts for {name}')
            layouts[name] = description
    return layouts


def _cell_fields(description):
    """Parse DD4hep field offsets/signed widths, checking layout consistency."""
    fields, offset, used = {}, 0, 0
    for token in description.split(','):
        parts = token.strip().split(':')
        if len(parts) not in (2, 3):
            raise ValueError(f'Invalid CellID field: {token}')
        name, width = parts[0], int(parts[-1])
        if len(parts) == 3:
            offset = int(parts[1])
        size = abs(width)
        if not name or not 0 < size <= 64 or offset < 0 or offset+size > 64:
            raise ValueError(f'Invalid CellID field: {token}')
        mask = ((1 << size)-1) << offset
        if name in fields or used & mask:
            raise ValueError(f'Duplicate/overlapping CellID field: {token}')
        fields[name] = (offset, width)
        used |= mask
        offset += size
    return fields


def decode_cell_ids(cells, description):
    """Decode DD4hep name:width / name:offset:width fields, including signed ones."""
    cells = np.asarray(cells)
    if cells.dtype.kind not in 'iu':
        raise TypeError('CellIDs must be exact integers')
    cells = cells.astype(np.uint64, copy=False)
    result = {}
    for name, (offset, width) in _cell_fields(description).items():
        size = abs(width)
        values = (cells >> np.uint64(offset)) & np.uint64((1 << size)-1)
        if width < 0:
            sign = np.uint64(1 << (size-1))
            values = ((values ^ sign)-sign).view(np.int64)
        # Nullable integers keep exact values when different readouts have
        # different fields and their tables are concatenated (no float casts).
        result[name] = pd.array(values, dtype='Int64' if width < 0 else 'UInt64')
    return result


def _field(record, name):
    matches = [key for key in record.fields if key == name or key.endswith('.'+name)]
    if len(matches) != 1:
        raise ValueError('Cannot resolve metadata field: '+name)
    return ak.to_list(record[matches[0]])


def _frame_encodings(root):
    layouts = {}
    for name in ('metadata', 'meta', 'runs'):
        if name not in root or 'GPStringKeys' not in root[name]:
            continue
        keys = ak.to_list(root[name]['GPStringKeys'].array(entry_stop=1))
        values = ak.to_list(root[name]['GPStringValues'].array(entry_stop=1))
        if not keys:
            continue
        for key, value in zip(keys[0], values[0]):
            if key.endswith('__CellIDEncoding'):
                description = value if isinstance(value, str) else value[0]
                collection = key[:-len('__CellIDEncoding')]
                if collection in layouts and layouts[collection] != description:
                    raise ValueError('Conflicting stored encoding: '+collection)
                layouts[collection] = description
    return layouts


def _column_name(field):
    names = {'cellID': 'cell_id', 'eDep': 'energy_deposit_GeV',
             'pathLength': 'path_length_mm', 'time': 'time_ns',
             'PDG': 'pdg', 'generatorStatus': 'generator_status',
             'simulatorStatus': 'simulator_status', 'mass': 'mass_GeV', 'charge': 'charge_e'}
    for axis in 'xyz':
        names['position.'+axis] = axis+'_mm'
        names['momentum.'+axis] = 'p'+axis+'_GeV'
        names['vertex.'+axis] = 'vertex_'+axis+'_mm'
        names['endpoint.'+axis] = 'endpoint_'+axis+'_mm'
        names['momentumAtEndpoint.'+axis] = 'endpoint_p'+axis+'_GeV'
    return names.get(field, field.replace('.', '_'))


def _surface_table(hits):
    """Catalogue observed detector layers; unhit geometry is not inferred."""
    if hits.empty:
        return pd.DataFrame(columns=SURFACE_KEY+['n_hits', 'n_particles'])
    grouped = hits.groupby(SURFACE_KEY, observed=True, dropna=False)
    counts = grouped.size().rename('n_hits')
    particles = hits.drop_duplicates(SURFACE_KEY+KEY).groupby(SURFACE_KEY, observed=True, dropna=False).size()
    return pd.concat([counts, particles.rename('n_particles')], axis=1).reset_index()


def iter_mc_tracking_hits(filename, *, collections=None, primary_only=True, pdgs=None,
                          eta_range=None, p_range=None, encodings=None, stop=None, step_size=5000):
    """Yield dicts of linked DataFrames, reading each selected branch once/block.

    MC selections use generatorStatus==1 when primary_only=True, optional exact
    signed PDG codes, and half-open [low,high) eta / |p| [GeV] ranges. All hits
    directly linked to those particles are returned; there are no hit/fit cuts.
    Every stored scalar SimTrackerHit field is retained, plus decoded CellID
    fields. MC parent/daughter relation arrays are not materialized.
    """
    if step_size < 1 or (stop is not None and stop < 1):
        raise ValueError('step_size and stop must be positive')
    for bounds in (eta_range, p_range):
        if bounds is not None and (len(bounds) != 2 or np.any(np.isnan(bounds)) or bounds[0] >= bounds[1]):
            raise ValueError('Selection ranges must be increasing [low,high) pairs')
    pdgs = None if pdgs is None else tuple(int(value) for value in np.atleast_1d(pdgs))
    source = Path(filename).resolve()
    with uproot.open(source, array_cache=None) as root:
        tree = root['events']
        info = root['podio_metadata']['events___CollectionTypeInfo'].array()[0]
        catalogue = {name: {'collection_id': int(cid), 'type': kind, 'subset': bool(subset)}
                     for name, cid, kind, subset in zip(_field(info, 'name'), _field(info, 'collectionID'),
                                                       _field(info, 'dataType'), _field(info, 'isSubset'))}
        available = [name for name, row in catalogue.items()
                     if row['type'] == 'edm4hep::SimTrackerHitCollection' and not row['subset']]
        selected = available if collections is None else ([collections] if isinstance(collections, str) else list(collections))
        if not selected or len(set(selected)) != len(selected) or any(name not in available for name in selected):
            raise ValueError(f'Choose distinct owning SimTrackerHit collections from {available}')
        layouts = _frame_encodings(root)
        for name, description in (encodings or {}).items():
            if name in layouts and _cell_fields(layouts[name]) != _cell_fields(description):
                raise ValueError('Provided/stored encoding mismatch: '+name)
            layouts[name] = description
        for name in selected:
            if name not in layouts:
                raise ValueError(f'No CellID encoding for {name}; pass encodings from the matching geometry')
            _cell_fields(layouts[name])
        file_id = str(root.file.uuid)
        mc_id = catalogue['MCParticles']['collection_id']
        branch_maps = {name: {key.removeprefix(name+'.'): tree[name][key] for key in tree[name].keys()}
                       for name in ['MCParticles', *selected, *['_'+name+'_particle' for name in selected]]}
        stop_entry = tree.num_entries if stop is None else min(stop, tree.num_entries)
        metadata = {'source': str(source), 'file_id': file_id, 'source_bytes': source.stat().st_size,
                    'source_mtime_ns': source.stat().st_mtime_ns, 'entries': tree.num_entries,
                    'collections': {name: catalogue[name] for name in ['MCParticles', *selected]},
                    'encodings': {name: layouts[name] for name in selected},
                    'selection': {'primary_only': primary_only, 'pdgs': None if pdgs is None else list(pdgs),
                                  'eta_range': eta_range, 'p_range_GeV': p_range},
                    'units': 'momenta/energy GeV; positions/path length mm; times ns; angles rad',
                    'particle_key': KEY, 'surface_key': SURFACE_KEY,
                    'hit_fields': {name: {_column_name(field): field for field in branch_maps[name]} for name in selected}}

        def read(name, field):
            return branch_maps[name][field].array(entry_start=start, entry_stop=end)

        for start in range(0, stop_entry, step_size):
            end = min(start+step_size, stop_entry)
            pdg = read('MCParticles', 'PDG')
            sizes = ak.to_numpy(ak.num(pdg))
            offsets = np.r_[0, np.cumsum(sizes)]
            mc_index = ak.to_numpy(ak.flatten(ak.local_index(pdg)))
            mc_events = np.repeat(np.arange(end-start), sizes)
            mc = {_column_name(field): ak.to_numpy(ak.flatten(pdg if field == 'PDG' else read('MCParticles', field)))
                  for field in branch_maps['MCParticles'] if not field.endswith(('_begin', '_end'))}
            p = np.sqrt(sum(mc['p'+axis+'_GeV'].astype(float)**2 for axis in 'xyz'))
            with np.errstate(divide='ignore', invalid='ignore'):
                eta = np.arcsinh(mc['pz_GeV']/np.hypot(mc['px_GeV'], mc['py_GeV']))
            keep_mc = (mc['generator_status'] == 1) if primary_only else np.ones(len(p), bool)
            if pdgs is not None:
                keep_mc &= np.isin(mc['pdg'], list(pdgs))
            for values, bounds in ((eta, eta_range), (p, p_range)):
                if bounds is not None:
                    keep_mc &= (values >= bounds[0]) & (values < bounds[1])
            particles = pd.DataFrame({name: values[keep_mc] for name, values in mc.items()})
            particles['file_id'] = pd.Categorical([file_id]*len(particles), categories=[file_id])
            particles['event_entry'] = (mc_events+start)[keep_mc]
            particles['mc_index'] = mc_index[keep_mc]
            particles['mc_collection_id'] = np.full(len(particles), mc_id, dtype=np.uint32)
            particles['p_GeV'], particles['eta'] = p[keep_mc], eta[keep_mc]
            particles['phi_rad'] = np.arctan2(mc['py_GeV'], mc['px_GeV'])[keep_mc]
            counts = np.zeros(len(p), dtype=np.int64)
            frames = []
            for name in selected:
                cells = read(name, 'cellID')
                hit_sizes = ak.to_numpy(ak.num(cells))
                events = np.repeat(np.arange(end-start), hit_sizes)
                relation = '_'+name+'_particle'
                target_index, target_id = read(relation, 'index'), read(relation, 'collectionID')
                if not (np.array_equal(ak.num(target_index), hit_sizes) and np.array_equal(ak.num(target_id), hit_sizes)):
                    raise ValueError('Hit/reference counts differ: '+name)
                indices, ids = ak.to_numpy(ak.flatten(target_index)), ak.to_numpy(ak.flatten(target_id))
                if np.any(ids != mc_id) or np.any((indices < 0) | (indices >= sizes[events])):
                    raise ValueError('Invalid hit -> MCParticle ObjectID: '+name)
                target = offsets[events]+indices
                keep = keep_mc[target]
                counts += np.bincount(target[keep], minlength=len(p))
                data = {}
                for field in branch_maps[name]:
                    values = cells if field == 'cellID' else read(name, field)
                    if not np.array_equal(ak.num(values), hit_sizes):
                        raise ValueError('Non-scalar hit attribute or inconsistent counts: '+name+'.'+field)
                    data[_column_name(field)] = ak.to_numpy(ak.flatten(values))[keep]
                hits = pd.DataFrame(data)
                hits['cell_id'] = hits.cell_id.astype(np.uint64)
                hits['file_id'] = pd.Categorical([file_id]*len(hits), categories=[file_id])
                hits['event_entry'] = events[keep]+start
                hits['mc_index'] = indices[keep]
                hits['mc_collection_id'] = ids[keep].astype(np.uint32)
                hits['hit_collection'] = pd.Categorical([name]*len(hits), categories=selected)
                hits['hit_collection_id'] = np.full(len(hits), catalogue[name]['collection_id'], dtype=np.uint32)
                hits['hit_index'] = ak.to_numpy(ak.flatten(ak.local_index(cells)))[keep]
                for field, values in decode_cell_ids(hits.cell_id.to_numpy(), layouts[name]).items():
                    hits[field if field in ('system', 'layer') else 'cell_'+field] = values
                for field in ('system', 'layer'):
                    if field not in hits:
                        hits[field] = pd.array([pd.NA]*len(hits), dtype='UInt64')
                frames.append(hits)
            particles['n_hits'] = counts[keep_mc]
            hits = pd.concat(frames, ignore_index=True)
            yield {'particles': particles, 'hits': hits, 'surfaces': _surface_table(hits),
                   'metadata': dict(metadata, entry_range=[start, end])}


def read_mc_tracking_hits(filename, **options):
    """Return tables for one file/prefix. Use the iterator for bounded-memory batches."""
    blocks = list(iter_mc_tracking_hits(filename, **options))
    if not blocks:
        raise ValueError('Input has no event entries')
    particles = pd.concat([block['particles'] for block in blocks], ignore_index=True)
    hits = pd.concat([block['hits'] for block in blocks], ignore_index=True)
    metadata = dict(blocks[0]['metadata'], entry_range=[0, blocks[-1]['metadata']['entry_range'][1]])
    return {'particles': particles, 'hits': hits, 'surfaces': _surface_table(hits), 'metadata': metadata}
