'''
    Podio-specific helpers for epic tracking analysis.
    Keep this module optional and import only in environments with podio.
    Shujie Li, Aug 2025
'''

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import fnmatch

from epic_analysis_base import (
    TRACK_HIT_COUNT_MIN_MIN,
    TRACK_HIT_COUNT_MIN,
    TRACK_MOM_MIN,
    TRACK_HIT_FRACTION_MIN,
    TRACK_HIT_COUNT_GHOST_MAX,
    VERTEX_CUT_R_MAX,
    VERTEX_CUT_Z_MAX,
    status_to_source,
)

def read_podio(fname, s3_dir="",tname="events"):
    """Read ROOT file with podio. Does NOT work for the metadata tree """
    # from podio import root_io
    import sys
    if 'podio' not in sys.modules:
        print("Loading podio ROOT IO reader(this will take ~2 minutes)...")

    from podio.root_io import Reader  # More specific
    
    server = 'root://dtn-eic.jlab.org//volatile/eic/'
    if len(s3_dir) > 1:
        fname = server + s3_dir + fname
    reader = Reader(fname)
    tree = reader.get("events")
    print(f"read_podio: read {fname}:{tname}. {len(tree)} events in total")
    return tree  

def show_getter_podio(collection):
    """Display basic info about collection or single object"""
    type_check = check_type(collection)
    
    if type_check.is_iterable:
        print(f"Number of objects: {len(collection)}")
        sample_obj = collection[0] if len(collection) > 0 else None
    else:
        print("Single object")
        sample_obj = collection
    
    if sample_obj:
        sample_type = type(sample_obj).__name__
        print(f"Object type: {sample_type}")
        getters = [m for m in dir(sample_obj) if m.startswith('get') and not m.startswith('_')]
        print(f"Available getter methods ({len(getters)}):")
        for getter in getters:
            print(f"  {getter}")

        # Pretty print version
        scalars = {m: getattr(sample_obj, m)() for m in dir(sample_obj) 
                if m.startswith('get') and not m.startswith('_') and callable(getattr(sample_obj, m)) 
                and isinstance(getattr(sample_obj, m)(), (int, float, bool, str))}
        for k, v in scalars.items(): print(f"{k}: {v}")

def show_collections_podio(event, pattern=None):
    """Show collections matching pattern (case-insensitive), or all if no pattern"""
    all_collections = event.getAvailableCollections()
    if pattern:
        collections = [name for name in all_collections 
                      if fnmatch.fnmatch(name.lower(), pattern.lower())]
    else:
        collections = all_collections
    # print(f"Collections ({len(collections)}):")
    # for name in collections:
    #     print(f"  {name}") 
    return collections

def get_collection_member_podio(podio_collection, member_name):
    """
    Access podio object members in uproot/awkward style
    
    Args:
        podio_collection: The podio branch (collection or single object)
        member_name: Member name (like "nHoles", "chi2")
    """
    type_check = check_type(podio_collection)
    if type_check.is_empty:
        raise ValueError("Collection is empty")
    
    # Helper function to process a single object
    def process_single_object(obj):
        # Common getter patterns in podio
        possible_getters = [
            f"get{member_name}",
            f"get{member_name.capitalize()}",
            f"getN{member_name.capitalize()}",  # for count-like members
            member_name
        ]
        
        for getter_name in possible_getters:
            if hasattr(obj, getter_name):
                attr = getattr(obj, getter_name)
                if callable(attr):
                    return attr()
                else:
                    return attr
        
        # If not found, show available members
        all_methods = [method for method in dir(obj) if not method.startswith('_')]
        getters = [method for method in all_methods if method.startswith('get')]
        other_attrs = [attr for attr in all_methods if not attr.startswith('get') and not callable(getattr(obj, attr, None))]
        
        error_msg = f"Cannot find member '{member_name}' in {type(obj).__name__}\n"
        error_msg += f"Available getter methods:\n"
        for getter in getters:
            error_msg += f"  {getter}\n"
        if other_attrs:
            error_msg += f"Available attributes:\n"
            for attr in other_attrs:
                error_msg += f"  {attr}\n"
        
        raise AttributeError(error_msg)
    
    if type_check.is_iterable:
        # Process collection - return list of results
        results = []
        for obj in podio_collection:
            results.append(process_single_object(obj))
        return results
    else: 
        # Process single object - return single result
        return process_single_object(podio_collection)

class PodioCollectionWrapper:
    def __init__(self, podio_collection):
        self.collection = podio_collection
        self.type_check = check_type(podio_collection)
    
    def __getitem__(self, member_name):
        return get_collection_member_podio(self.collection, member_name)
    
    def __len__(self):
        return len(self.collection) if self.type_check.is_iterable else 1
    
    def __iter__(self):
        return iter(self.collection) if self.type_check.is_iterable else iter([self.collection])

def check_type(obj):
    """Check object type with comprehensive categorization"""    
    def get_value_category(value):
        # Check for numbers first
        if isinstance(value, (int, float, complex, bool)):
            return 'number'
        # Check for numpy arrays
        elif isinstance(value, np.ndarray):
            if value.ndim == 1:
                return 'array'
            else:
                return 'nested_array'  # 2D, 3D arrays etc.
        # Check for Python lists/tuples
        elif isinstance(value, (list, tuple)):
            if len(value) > 0:
                # Check if it's a list of arrays (nested structure)
                first_item = value[0]
                if isinstance(first_item, (np.ndarray, list, tuple)):
                    return 'nested_array'
                else:
                    return 'list'
            else:
                return 'list'
        # Check for podio RelationRange or other iterable collections
        elif 'RelationRange' in str(type(value)) or (hasattr(value, '__len__') and hasattr(value, '__iter__') and not isinstance(value, str)):
            return 'range'
        # Everything else is an object
        else:
            return 'object'

    # Create a simple container
    class CategoryValue:
        def __init__(self, val):
            self.value = val
            cat = get_value_category(val)
            self.category = cat
            self.size        = getattr(val, '__len__', lambda: 1)()
            self.is_empty    = getattr(val, '__len__', lambda: 1)() == 0

            self.is_range  = (cat == "range")
            self.is_number = (cat == "number") 
            self.is_object = (cat == "object")
            self.is_array  = (cat == "array")
            self.is_nested_array = (cat == "nested_array")
            self.is_list   = (cat == "list")
            # Convenience groupings
            self.is_iterable = cat in ["range", "array", "nested_array", "list"]
            self.is_numpy    = cat in ["array", "nested_array"]
            self.is_simple   = cat in ["number", "list", "array"]
            
    return CategoryValue(obj)


def build_rawhit_lookup(obj_list):
    """Build fast lookup using cellID + time + charge as unique key"""
    return {(obj.getCellID(), obj.getTimeStamp(), obj.getCharge()): i 
            for i, obj in enumerate(obj_list)}

def get_rawhit_index(lookup,obj):
    key = (obj.getCellID(), obj.getTimeStamp(), obj.getCharge())
    return lookup.get(key, -1)

def is_valid_podio_object(obj):
    """Return False for a null or unresolved PODIO relation."""
    if obj is None:
        return False
    try:
        object_id = obj.id()
    except (AttributeError, ReferenceError, RuntimeError):
        return False
    return object_id.index >= 0 and object_id.collectionID != 0xFFFFFFFF

def build_obj_lookup(obj_list):
    """Build lookup that handles multiple objects with same collectionID+index"""
    lookup = {}
    for i, obj in enumerate(obj_list):
        if not is_valid_podio_object(obj):
            continue
        key = (obj.id().collectionID, obj.id().index)
        if key not in lookup:
            lookup[key] = []
        lookup[key].append(i)
    return lookup

def get_obj_indices(lookup, obj):
    """Return list of all indices matching the object's key"""
    if not is_valid_podio_object(obj):
        return []
    key = (obj.id().collectionID, obj.id().index)
    return lookup.get(key, [])

def get_traj_hits(event,bname="CentralCKFTrajectories",kcombine=0):
    ## ALL raw (one) --> sim hit (many, with weight) associations for a given event. Noise hits won't have association.
    ## ----------------
    # traj-based info
    ## ----------------
    ## keep track of the rawhit(simhit) index in the association, which includes all central tracker hits in one collection.
    asso=event.get("CentralTrackingRawHitAssociations") 
    asso_raw=PodioCollectionWrapper(asso)["RawHit"]
    lookup = build_obj_lookup(asso_raw)
    br    = event.get(bname) ## for one event, can have multiple subentries
    vname = "Measurements_deprecated"
    ltraj=[]
    lhit=[]
    lweight=[]
    # lpos =[]
    lpart=[]
    lrec_hit_col=[]
    lrec_hit_id=[]
    lmeasurement_col=[]
    lmeasurement_id=[]
    lhit_in_measurement=[]
    def add_invalid_hit(traj_id, measurement_col, measurement_id, hit_number,
                        rec_col=-1, rec_id=-1, association_id=-1):
        """Keep an invalid or noise hit in the trajectory accounting."""
        ltraj.append(traj_id)
        lhit.append(association_id)
        lpart.append(-1)
        lrec_hit_col.append(rec_col)
        lrec_hit_id.append(rec_id)
        lmeasurement_col.append(measurement_col)
        lmeasurement_id.append(measurement_id)
        lhit_in_measurement.append(hit_number)

    ## for each traj--> each measurement --> rec hit --> raw hit --> match raw with sim by association index --> particle
    for ii,traj in enumerate(br):
        measurements  = PodioCollectionWrapper(traj)[vname]
        for measurement in measurements:
            if not is_valid_podio_object(measurement):
                add_invalid_hit(ii, -1, -1, -1)
                continue

            measurement_col = measurement.id().collectionID
            measurement_id  = measurement.id().index
            hits = measurement.getHits()
            if len(hits) == 0:
                add_invalid_hit(ii, measurement_col, measurement_id, -1)
                continue

            # A Measurement2D can contain more than one constituent hit. In
            # particular, TOF clusters commonly contain two shared rec hits.
            for hit_number, hit in enumerate(hits):
                if not is_valid_podio_object(hit):
                    add_invalid_hit(
                        ii, measurement_col, measurement_id, hit_number
                    )
                    continue

                rec_col = hit.id().collectionID
                rec_id  = hit.id().index
                raw = hit.getRawHit() # rec2raw is one-to-one
                if not is_valid_podio_object(raw):
                    add_invalid_hit(
                        ii, measurement_col, measurement_id, hit_number,
                        rec_col, rec_id
                    )
                    continue

                # One raw hit can occur more than once when it is associated
                # with several simulated hits.
                indx = get_obj_indices(lookup, raw)
                if len(indx) == 0:
                    # A valid rec hit without a truth association is noise.
                    add_invalid_hit(
                        ii, measurement_col, measurement_id, hit_number,
                        rec_col, rec_id
                    )
                    continue

                for ind in indx:
                    sim = asso[ind].getSimHit()
                    if not is_valid_podio_object(sim):
                        add_invalid_hit(
                            ii, measurement_col, measurement_id, hit_number,
                            rec_col, rec_id, ind
                        )
                        continue

                    part = sim.getParticle()
                    if not is_valid_podio_object(part):
                        add_invalid_hit(
                            ii, measurement_col, measurement_id, hit_number,
                            rec_col, rec_id, ind
                        )
                        continue

                    ltraj.append(ii)
                    lhit.append(ind)
                    lpart.append(part.id().index)
                    lrec_hit_col.append(rec_col)
                    lrec_hit_id.append(rec_id)
                    lmeasurement_col.append(measurement_col)
                    lmeasurement_id.append(measurement_id)
                    lhit_in_measurement.append(hit_number)
    traj_hits=pd.DataFrame({
        "traj_id": ltraj,
        "part_id": lpart,
        "asso_hit": lhit,
        "measurement_col": lmeasurement_col,
        "measurement_id": lmeasurement_id,
        "hit_in_measurement": lhit_in_measurement,
        "rec_hit_col": lrec_hit_col,
        "rec_hit_id": lrec_hit_id,
    })
    #"position":lpos, "hit_weight":lweight})

    ## FIXME: check for overlapped tracks (for now ambiguity solver config won't allow sharing hits)
    # reoccur_hit=traj_hits.groupby('asso_hit').filter(lambda group: len(group) > 1)
    # for row in reoccur_hit.itertuples():
    #     print(f'WARNING: duplihits detected:', row)

    ## if one rec hit is associated to multiple sim hits, but all sim hits go to the same particle, then only keep one entry. 
    ## FIXME: not sure if this will work if we allow overlap traj hit. 
    traj_hits = traj_hits.drop_duplicates(
        subset=[
            "part_id", "traj_id", "measurement_col", "measurement_id",
            "hit_in_measurement", "rec_hit_col", "rec_hit_id",
        ],
        keep="first",
    )
    traj_hits['weight'] = (traj_hits.groupby(['traj_id', 'part_id']).transform('size') / 
                        traj_hits.groupby('traj_id').transform('size'))

    if kcombine:
        traj_hits = traj_hits.groupby(['traj_id', 'part_id'], as_index=False).agg({
            'asso_hit': list,
            'weight': 'first'  # Keep the weight value (should be same for same traj+particle)
        })

    return traj_hits

def get_part_hits(event, traj_hits, ksignal=0):
    ## ----------------
    # particle-based info
    ## ----------------
    ltraj_id=[]
    lpart_id=[]
    lsimhit=[]
    lgenID=[]
    lraw_hit_col=[]
    lraw_hit_id=[]
    # lpart=[]
    # lposition=[]

    ## for fast lookup (instead of traj_hits.asso_hit==ii). Assume no shared hits
    # Invalid or noise hits use asso_hit=-1 and intentionally have no entry
    # in the truth-association lookup.
    traj_base = traj_hits[traj_hits["asso_hit"] >= 0].set_index("asso_hit")
    if not traj_base.index.is_unique:
        raise ValueError("get_part_hits: duplicate asso_hit in traj_hits")
    traj_map = traj_base['traj_id']
    ## for each simhit (that is converted to rec hit-->measurement candidate), find related particle and traj
    asso=event.get("CentralTrackingRawHitAssociations")
    for ii,association in enumerate(asso):
        sim  = association.getSimHit()
        if not is_valid_podio_object(sim):
            continue
        part = sim.getParticle()
        if not is_valid_podio_object(part):
            continue
        raw  = association.getRawHit()
        if not is_valid_podio_object(raw):
            continue
        # cond_vertex   =  (np.sqrt(part.getVertex().x**2+part.getVertex().y**2)<1 )&(abs(part.getVertex().z)<100)
        # cond = cond_vertex
        status = part.getGeneratorStatus()
        cond = (status in (1, 2)) if ksignal == 1 else True
        if cond:
            lsimhit.append(ii) ## as before, use unique index (and unique sim hit) from the association. 
            # lpart.append(part)
            lpart_id.append(part.id().index)
            lgenID.append(status)
            lraw_hit_col.append(raw.id().collectionID)
            lraw_hit_id.append(raw.id().index)
            # lposition.append(sim.getPosition())

            ## find which trajectory used this hit
            ltraj_id.append(int(traj_map.get(ii, -1)))
    part_hits = pd.DataFrame({"part_id":lpart_id, "part_status":lgenID, "asso_hit":lsimhit, "traj_id":ltraj_id, 
                            "raw_hit_col": lraw_hit_col, "raw_hit_id": lraw_hit_id})#, "position":lposition, "particle": lpart})
    part_hits = part_hits.drop_duplicates(subset=["part_id", "traj_id", "raw_hit_col", "raw_hit_id"],  keep="first")

    # df = primary_hits[primary_hits.groupby('particle')['particle'].transform('count') >= 3]
    return part_hits

def get_traj_purity(traj_hits):
    '''
    traj_hits is the output from get_traj_hits()
    returns trajectory and source, max fraction=purity
    '''
    grouped = traj_hits.groupby(['traj_id'])
    # Analyze each group
    result = grouped['part_id'].agg([
        ('total_count', 'count'),
        ('unique_source', lambda x: x.nunique()),
        ('most_common_source', lambda x: x.value_counts().idxmax()),
        ('max_count', lambda x: x.value_counts().max())
    ])

    # Calculate derived columns
    result['max_fraction'] = result['max_count'] / result['total_count']
    # result['all_same'] = result['unique_source'] == 1
    return result.copy().reset_index()


def plot_part_traj_flow(df, params=None, mcpart=None):
    """Create alluvial-style diagram showing particle-trajectory flows"""
    # dict_part = dict(zip(df["part_id"], df["particle"]))

    # Separate used and unused hits
    df_used = df[df['traj_id'] != -1]  # Only hits used in trajectories
    df_all = df  # All hits including unused
    
    # Calculate statistics
    particle_totals = df_all.groupby('part_id')['asso_hit'].apply(len)    # Total hits per particle
    particle_used   = df_used.groupby('part_id')['asso_hit'].apply(len)   # Used hits per particle
    traj_totals     = df_used.groupby('traj_id')['asso_hit'].apply(len)   # Total hits per trajectory
    
    # Only include particles that have some used hits (have flows)
    particles_with_flows = df_used['part_id'].unique()
    trajectories = sorted(df_used['traj_id'].unique())  # Only trajectories with used hits
    len1 = len(particles_with_flows)
    len2 = len(trajectories)
    # Prepare flow data (only for used hits)
    flows = df_used.groupby(['part_id', 'traj_id'])['asso_hit'].apply(len).reset_index()
    flows.columns = ['part_id', 'traj_id', 'hits']
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Position particles on left, trajectories on right
    particle_y = {p: i for i, p in enumerate(sorted(particles_with_flows))}
    traj_y     = {t: i for i, t in enumerate(trajectories)}
    
    # Draw flows as curved lines
    for _, row in flows.iterrows():
        particle = row['part_id']
        traj     = row['traj_id']
        hits     = row['hits']
        
        # Start and end points
        x1, y1 = 0, len1 - particle_y[particle]
        x2, y2 = 1, len1 - traj_y[traj]
        
        # Create curved line
        x_curve = [x1, 0.5, x2]
        y_curve = [y1, (y1 + y2) / 2, y2]
        
        # Line thickness proportional to hits
        linewidth = max(1, hits / flows['hits'].max() * 10)
        
        ax.plot(x_curve, y_curve, linewidth=linewidth, alpha=0.6)
    
    # Add particle labels (only for particles with flows)
    ax.text(-0.25, len1+1,  f'Particle ID:  (used/total hits)',  ha='left', va='center')
    ax.text(0.9,   len1+1,  f'Trajectory ID: (nMeasurements)',  ha='left', va='center')

    for i, p in enumerate(sorted(particles_with_flows)):
        total_hits = particle_totals[p]
        used_hits = particle_used[p]
        color = 'k'
        if total_hits<TRACK_HIT_COUNT_MIN:
            color="grey"
        if mcpart is not None:
            pp=mcpart.iloc[int(p)]
            if abs(pp["vertex_r"])>VERTEX_CUT_R_MAX or abs(pp["vertex.z"])>VERTEX_CUT_Z_MAX or pp["mom"]<TRACK_MOM_MIN:
                color='grey'
            ax.text(-0.03, len1-i,  f'#{p}: {status_to_source[pp.generatorStatus]} ({used_hits}/{total_hits})',
                ha='right', va='center',color=color)
        else:
            ax.text(-0.03, len1-i,  f'#{p}:   ({used_hits}/{total_hits})',
                ha='right', va='center',color=color)
            
    # Add trajectory labels
    for i, t in enumerate(trajectories):
        total_traj_hits = traj_totals[t]
        color="k"
        if total_traj_hits<TRACK_HIT_COUNT_MIN:
            color="grey"
        if params is not None: 
            p = params.iloc[int(t)]
            if abs(p["loc.a"])>VERTEX_CUT_R_MAX or abs(p["loc.b"])>VERTEX_CUT_Z_MAX or p["mom"]<TRACK_MOM_MIN:
                color='grey'
                # print(t, ["loc.a"], p["loc.b"])
        ax.text(1.03, len1-i, f'#{int(t)} ({total_traj_hits})', 
                ha='left', va='center', color=color)
    
    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(-0.5, max(len(particles_with_flows), len(trajectories)) +3)
    ax.set_title('Particle to Trajectory Flow (Particles with Used Hits Only)')
    ax.axis('off')
    
    return plt


def get_part_traj_counts(event,mcpart, ksignal=0, kverbose=0):
    traj_hits=get_traj_hits(event)
    part_hits=get_part_hits(event,traj_hits, ksignal)

    # Prefer explicit particle-id columns; fall back to positional index for compatibility.
    id_col = None
    if "orig_subentry" in mcpart.columns:
        id_col = "orig_subentry"
    elif "subentry" in mcpart.columns:
        id_col = "subentry"

    ## -----------Find good track-----------
    ## Get counts per (traj_id, part_id) and total per traj_id
    traj_counts = get_traj_purity(traj_hits)
    ## get generator status
    if id_col is not None:
        status_map = (
            mcpart.drop_duplicates(subset=[id_col], keep="first")
            .set_index(id_col)["generatorStatus"]
        )
        traj_counts["part_status"] = (
            traj_counts["most_common_source"].map(status_map).fillna(-1).astype(int)
        )
    else:
        traj_counts["part_status"]=traj_counts["most_common_source"].apply(lambda x: mcpart.iloc[x].generatorStatus)

    ## -----------Find good particles-----------
    ## track hit cut
    part_counts = part_hits.groupby("part_id").size()
    part_counts = part_counts[part_counts>=TRACK_HIT_COUNT_MIN_MIN]
    if id_col is not None:
        mcpart_hits = mcpart[mcpart[id_col].isin(part_counts.index)].copy()
        mcpart_hits["hit_counts"] = mcpart_hits[id_col].map(part_counts).astype(int)
    else:
        mcpart_hits = mcpart.iloc[part_counts.index].copy()
        mcpart_hits["hit_counts"] = part_counts.values
    # traj_counts["part_status"]=traj_counts["most_common_source"].apply(lambda x: mcpart.iloc[x].generatorStatus)
    ## only do event-by-event quality check when required. Otherwise return the dataframe for further analysis
    if kverbose:
        good_part_id   = (part_hits.groupby("part_id").size()>=TRACK_HIT_COUNT_MIN)
        part_hits_good = part_hits[part_hits["part_id"].isin(good_part_id[good_part_id].index)][["part_id","part_status"]].drop_duplicates()
        if id_col is not None:
            good_mcpart = mcpart[mcpart[id_col].isin(part_hits_good.part_id.unique())].copy()
        else:
            good_mcpart = mcpart.iloc[part_hits_good.part_id.unique()]
        ## vertex and momentum cut
        cond_vertex = (abs(good_mcpart.vertex_r)<VERTEX_CUT_R_MAX)&(abs(good_mcpart["vertex.z"])<VERTEX_CUT_Z_MAX)
        cond_mom    = (good_mcpart.mom>TRACK_MOM_MIN)
        good_mcpart = good_mcpart[cond_vertex&cond_mom]
        ## signal or background
        cond_sig    = (good_mcpart.generatorStatus==1)|(good_mcpart.generatorStatus==2)
        good_mcpart_sig  =good_mcpart[cond_sig]
        good_mcpart_other=good_mcpart[(~cond_sig)]
        print("Number of good particles (signal, other):",len(good_mcpart_sig), len(good_mcpart_other))


    if kverbose:
        traj_counts["traj_status"]=0
        traj_counts.loc[(traj_counts.max_fraction<TRACK_HIT_FRACTION_MIN) | (traj_counts.max_count<=TRACK_HIT_COUNT_GHOST_MAX),"traj_status"]=-1 ## ghost status=-1
        traj_counts.loc[(traj_counts.traj_status>-1)&(traj_counts['total_count']>=TRACK_HIT_COUNT_MIN),'traj_status']=1

        ghost_traj_id = traj_counts[traj_counts.traj_status==-1].index.tolist()
        print("list of ghost track id:", ghost_traj_id)
        good_traj=traj_counts[traj_counts.traj_status==1].copy()
        good_traj_sig=good_traj[(good_traj.part_status==1)|(good_traj.part_status==2)]
        ntraj_sig   = len(good_traj_sig)
        ntraj_other = len(good_traj)-ntraj_sig
        print("Number of good tracks from sig/others:",ntraj_sig, ntraj_other)
    ## hit purity for good tracks
    ## Percentage of track hits from a single source 
        purity_hit=good_traj.max_fraction.to_list()
        print("Hit purity:" , purity_hit)
    ##-----------Track to particle efficiency----------------
    # fraction of good signal particles that are linked to some good track
        if id_col is not None:
            good_part_sig_id = good_mcpart_sig[id_col].tolist()
        else:
            good_part_sig_id = good_mcpart_sig.index.get_level_values('subentry').tolist()
        good_traj_id = good_traj.most_common_source.unique()
        good_traj_sig_id = good_traj_sig.most_common_source.unique()
        common = list(set(good_part_sig_id) & set(good_traj_id))
        print(set(good_part_sig_id),  set(good_traj_sig_id), common)
        ## 50% of the total particle hits go to that traj <---- skip this for now for a looser check
        # good_traj['fract'] = good_traj['max_count'] / good_traj['most_common_source'].map(part_hits['part_id'].value_counts())
        # good_traj=good_traj.fract>=0.5
    return mcpart_hits, traj_counts



__all__ = [
    "read_podio",
    "show_getter_podio",
    "show_collections_podio",
    "get_collection_member_podio",
    "PodioCollectionWrapper",
    "check_type",
    "build_rawhit_lookup",
    "get_rawhit_index",
    "is_valid_podio_object",
    "build_obj_lookup",
    "get_obj_indices",
    "get_traj_hits",
    "get_part_hits",
    "get_traj_purity",
    "plot_part_traj_flow",
    "get_part_traj_counts",
]
