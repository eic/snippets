'''
    Utility functions for epic tracking analysis. 
    See also epic_analysis_podio.py
    Shujie Li, Aug 2025
'''

# %load_ext memory_profiler
import numpy as np
import pandas as pd
import seaborn as sns

import awkward as ak
import uproot as ur
# import podio
# from podio import root_io

import time
from fnmatch import fnmatch
import types
from particle import Particle

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import matplotlib as mpl

def configure_analysis_environment(
    apply_pandas=True,
    apply_matplotlib=True,
    apply_sns=True
):
    """Apply pandas/matplotlib defaults for interactive analysis."""
    if apply_pandas:
        pd.options.display.max_rows = 200
        pd.options.display.min_rows = 20
        pd.options.display.max_columns = 100
    if apply_matplotlib:
        plt.rcParams['figure.figsize'] = [8.0, 6.0]
        plt.rcParams['ytick.direction'] = 'in'
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['xaxis.labellocation'] = 'right'
        plt.rcParams['yaxis.labellocation'] = 'top'
        small_size = 10
        medium_size = 12
        bigger_size = 20
        plt.rc('font', size=small_size)
        plt.rc('axes', titlesize=medium_size)
        plt.rc('axes', labelsize=medium_size)
        plt.rc('xtick', labelsize=medium_size)
        plt.rc('ytick', labelsize=medium_size)
        plt.rc('legend', fontsize=small_size)
        plt.rc('figure', titlesize=bigger_size)
    if apply_sns:
        sns.set_theme(
            style='whitegrid',
            context='notebook',
            palette='bright',
            font_scale=1.0,
            rc={'figure.figsize': (6, 4)},
        )

# Constants
deg2rad = np.pi/180.0

## event source 
status_to_source = {
    1: "DIS", 2: "DIS",
    2001: "SR", 2002: "SR", 
    3001: "Bremstrahlung", 3002: "Bremstrahlung",
    4001: "Coulomb", 4002: "Coulomb",
    5001: "Touschek", 5002: "Touschek",
    6001: "Proton beam gas", 6002: "Proton beam gas"
}

## track quality cuts
TRACK_HIT_COUNT_MIN_MIN   = 4 ## absolute min to form a track with CKF
TRACK_HIT_COUNT_MIN       = 4
TRACK_MOM_MIN             = 0.2
TRACK_PT_MIN              = 0.2
TRACK_HIT_FRACTION_MIN    = 0.5
TRACK_HIT_COUNT_GHOST_MAX = 2
VERTEX_CUT_R_MAX          = 1
VERTEX_CUT_Z_MAX          = 100#mm

# Detector geometry definitions (unchanged)
barrel_range = [(30,42),(46,60),(115,130),(250,290),(400,450),(540,600),(620,655),(700,760)]
barrel_name = ["L0","L1","L2","L3","L4","inner MPGD","TOF","outer MPGD"]
name_sim_barrel = ["VertexBarrelHits","VertexBarrelHits","VertexBarrelHits","SiBarrelHits","SiBarrelHits","MPGDBarrelHits","TOFBarrelHits","OuterMPGDBarrelHits"]
name_rec_barrel = ["SiBarrelVertexRecHits","SiBarrelVertexRecHits","SiBarrelVertexRecHits","SiBarrelTrackerRecHits","SiBarrelTrackerRecHits","MPGDBarrelRecHits","TOFBarrelRecHits","OuterMPGDBarrelRecHits"]

disk_range = [(-1210.0, -1190.0), (-1110.0, -1090.0),(-1055.0, -1000.0), (-860.0, -840.0),
 (-660.0, -640.0), (-460.0, -440.0), (-260.0, -240.0), (240.0, 260.0),
 (440.0, 460.0), (690.0, 710.0), (940.0, 960.0), (1150.0, 1250.0),
 (1480.0, 1500.0), (1600.0, 1620.0), (1840.0, 1860.0), (1865.0, 1885.0)]
disk_name = ["E-MPGD Disk2","E-MPGD Disk 1","E-Si Disk 4","E-Si Disk 3","E-Si Disk 2","E-Si Disk 1","E-Si Disk 0",
                "H-Si Disk 0","H-Si Disk 1","H-Si Disk 2","H-Si Disk 3","H-Si Disk 4","H-MPGD Disk 1","H-MPGD Disk 2", "H-TOF Disk1","H-TOF Disk2"]
name_rec_disk = ["BackwardMPGDEndcapRecHits","BackwardMPGDEndcapRecHits",
               "SiEndcapTrackerRecHits","SiEndcapTrackerRecHits","SiEndcapTrackerRecHits","SiEndcapTrackerRecHits","SiEndcapTrackerRecHits","SiEndcapTrackerRecHits","SiEndcapTrackerRecHits","SiEndcapTrackerRecHits","SiEndcapTrackerRecHits","SiEndcapTrackerRecHits",
               "ForwardMPGDEndcapRecHits","ForwardMPGDEndcapRecHits",
               "TOFEndcapRecHits","TOFEndcapRecHits"]
name_sim_disk = ["BackwardMPGDEndcapHits","BackwardMPGDEndcapHits",
               "TrackerEndcapHits","TrackerEndcapHits","TrackerEndcapHits","TrackerEndcapHits","TrackerEndcapHits","TrackerEndcapHits","TrackerEndcapHits","TrackerEndcapHits","TrackerEndcapHits","TrackerEndcapHits","ForwardMPGDEndcapHits","ForwardMPGDEndcapHits",
               "TOFEndcapHits","TOFEndcapHits"]

# ACTS geometry ID masks
geo_mask_dict = {
    "approach": 0x0000000ff0000000,
    "boundary": 0x00ff000000000000,
    "extra": 0x00000000000000ff,
    "layer": 0x0000fff000000000,
    "sensitive": 0x000000000fffff00,
    "volume": 0xff00000000000000
}
geo_mask_values = types.MappingProxyType(geo_mask_dict)

# Global variables for caching
COL_TABLE = {}
CACHED_DATA = {}

def ak_flat(ak_array):
    return ak.to_numpy(ak.flatten(ak_array,axis=0))

def ak_df(ak_array):
    return ak.to_dataframe(ak_array)

def ak_hist(ak_array, **kwargs):
    return plt.hist(ak_flat(ak_array), **kwargs)

def ak_filter(br, cond, field=None):
    filtered = br[cond]
    return filtered[field] if field else filtered


def ak_sns(ak_array, **kwargs):
    """Histogram helper for awkward arrays using seaborn."""
    if isinstance(ak_array, (tuple, list)) and len(ak_array) == 2:
        x_data = ak_flat(ak_array[0])
        y_data = ak_flat(ak_array[1])
        kwargs.pop('element', None)
        kwargs.pop('fill', None)
        return sns.histplot(x=x_data, y=y_data, **kwargs)
    return sns.histplot(ak_flat(ak_array), element="step", fill=False, **kwargs)

def get_pdg_info(PDG):
    """Get particle info from PDG code"""
    try:
        return Particle.from_pdgid(PDG)
    except Exception:
        if PDG == 9902210:
            return Particle.from_pdgid(2212)
        print(f"ERROR (get_pdg_info): unknown PDG ID {PDG}")
        return Particle.empty()

def get_geoID(geoID, name="layer"):
    """Extract geometry ID components"""
    kMask = geo_mask_values[name]
    shift = 0
    mask_temp = kMask
    while (mask_temp & 1) == 0:
        mask_temp >>= 1
        shift += 1
    return (geoID & kMask) >> shift


def theta2eta(xx, inverse=0):
    """Convert theta to eta or vice versa"""
    if type(xx)==list:
        xx = np.array(xx)
    if inverse==1:
        return np.arctan((np.e)**(-xx))*2
    else:
        return -np.log(np.tan(xx/2.))

def select_string(strings, patterns):
    """Select strings matching patterns with wildcards"""
    if not isinstance(patterns, list):
        raise ValueError("The 'patterns' argument must be a list.")
    
    patterns = [pattern.lower() for pattern in patterns]
    return [s for s in strings if any(fnmatch(s.lower(), pattern) for pattern in patterns)]

def read_ur(fname, tname, s3_dir="", entry_start=0, entry_stop=None, return_range=False):
    """Read ROOT file with uproot
    fname: path to file
    tname: tree name
    s3_dir: if provided, try the BNL and JLab simulation campaign locations"""
    if len(s3_dir) > 0:
        # The same campaign path is stored below different roots at BNL and JLab.
        servers = (
            'root://epicxrd1.sdcc.bnl.gov:1095//eic/',
            'root://dtn-eic.jlab.org//volatile/eic/',
            "root://epicxrd1.sdcc.bnl.gov:1095/",
        )
        remote_path = s3_dir.rstrip('/') + '/' + fname.lstrip('/')
        open_errors = []
        last_error = None
        for server in servers:
            candidate = server + remote_path
            try:
                root_file = ur.open(candidate, timeout=5)
                fname = candidate
                break
            except OSError as error:
                last_error = error
                open_errors.append(f"{candidate}: {error}")
        else:
            attempted = '\n  '.join(open_errors)
            raise FileNotFoundError(
                f"read_ur: remote file was not found on either server:\n  {attempted}"
            ) from last_error
    else:
        root_file = ur.open(fname)

    tree = root_file[tname]
    if entry_stop is None or entry_stop == -1:
        entry_stop = tree.num_entries
    entry_stop = min(entry_stop, tree.num_entries)
    if entry_start < 0 or entry_start >= entry_stop:
        raise ValueError(f"read_ur: invalid entry range {entry_start}:{entry_stop}")
    print(
        f"read_ur: read {fname}:{tname}. "
        f"{tree.num_entries} events total; using [{entry_start}, {entry_stop})"
    )
    tree._entry_start = entry_start
    tree._entry_stop = entry_stop
    if return_range:
        return tree, entry_start, entry_stop
    return tree

def get_col_table(fname, s3_dir="", verb=0):
    """Get collection table from metadata"""
    global COL_TABLE
    meta = read_ur(fname, "podio_metadata", s3_dir)
    if "events___idTable" in meta.keys(): ## < eic-shell 25.09
        col_name = np.array(meta["m_names"].array()[0])
        col_id = np.array(meta["m_collectionIDs"].array()[0])
    else:
        col_id = get_branch_df(meta,"events___CollectionTypeInfo/events___CollectionTypeInfo.collectionID")["values"].tolist()
        col_name = get_branch_df(meta,"events___CollectionTypeInfo/events___CollectionTypeInfo.name")["values"].tolist()

    COL_TABLE = {}
    for ii, nn in zip(col_id, col_name):
        if verb:
            print(ii, nn)
        COL_TABLE[ii] = nn
    return COL_TABLE


# ============= BRANCH READING with ak or df =============

def get_branch_ak(tree, bname="", entry_start=0, entry_stop=-1, 
                           fields_subset=None, chunk_size=1000, verb=0):
    """Optimized branch reading with awkward arrays"""
    if bname not in tree.keys():
        raise KeyError(f"get_branch_ak: can't find branch {bname}")
    if verb:    
        print(f"Reading branch: {bname}")
    start_time = time.time()
    
    # Determine actual entry range
    if entry_start == 0 and hasattr(tree, "_entry_start"):
        entry_start = tree._entry_start
    if entry_stop == -1:
        entry_stop = tree._entry_stop if hasattr(tree, "_entry_stop") else tree.num_entries
    
    total_entries = entry_stop - entry_start
    if verb:
        print(f"Reading {total_entries} entries")
    
    # For large datasets, read in chunks
    if total_entries > chunk_size:
        if verb:
            print(f"Using chunked reading with chunk_size={chunk_size}")
        all_data = []
        
        for chunk_start in range(entry_start, entry_stop, chunk_size):
            chunk_end = min(chunk_start + chunk_size, entry_stop)
            if verb:
                print(f"  Reading chunk: {chunk_start} to {chunk_end}")
            
            chunk_data = tree[bname].array(
                library="ak",
                entry_start=chunk_start,
                entry_stop=chunk_end
            )
            all_data.append(chunk_data)
        
        # Concatenate all chunks
        ak_data = ak.concatenate(all_data)
    else:
        # Read all at once for smaller datasets
        ak_data = tree[bname].array(
            library="ak",
            entry_start=entry_start,
            entry_stop=entry_stop
        )
    
    read_time = time.time()
    if verb:
        try:
            size_bytes = ak.nbytes(ak_data)
        except Exception:
            size_bytes = None
        if size_bytes is not None:
            print(f"Awkward read: {read_time - start_time:.2f}s ({size_bytes/1e6:.1f} MB)")
        else:
            print(f"Awkward read: {read_time - start_time:.2f}s")
    
    # Rename fields by dropping branch prefix
    if hasattr(ak_data, 'fields'):
        if not ak_data.fields:
            ## return a single array
            return ak_data
        renamed_fields = {}
        for field in ak_data.fields:
            if "[" in field: ## drop nested array such as CentralCKFTrackParameters.covariance.covariance[21]
                continue
            if field.startswith(f'{bname}.'):
                new_name = field.replace(f'{bname}.', '')
                renamed_fields[new_name] = ak_data[field]
            else:
                renamed_fields[field] = ak_data[field]
        
        ak_data = ak.zip(renamed_fields)
        if verb:
            print(f"Renamed {len(renamed_fields)} fields")
    
    # Extract subset of fields if specified
    if fields_subset and hasattr(ak_data, 'fields'):
        subset_data = {}
        for field in fields_subset:
            if field in ak_data.fields:
                subset_data[field] = ak_data[field]
        ak_data = ak.zip(subset_data)
        if verb:
            print(f"Extracted subset: {fields_subset}")
    
    total_time = time.time()
    if verb:
            print(f"Total time: {total_time - start_time:.2f}s")
    
    return ak_data


def get_branch_df(tree, bname="", entry_start=0, entry_stop=-1, 
                           fields_subset=None, chunk_size=1000, verb=0):
    """Get branch as DataFrame when needed (for compatibility)"""
    ak_data = get_branch_ak(tree, bname, entry_start, entry_stop, 
                                     fields_subset, chunk_size)
    
    if verb:
        print("Converting to DataFrame...")
    start_time = time.time()
    
    try:
        df = ak.to_dataframe(ak_data)
        df = df.reset_index()  # Flatten the multi-index
        
        convert_time = time.time() - start_time
        if verb:
            print(f"DataFrame conversion: {convert_time:.2f}s")
            print(f"DataFrame shape: {df.shape}")
        
        return df
    
    except Exception as e:
        print(f"DataFrame conversion failed: {e}")
        print("Returning awkward array instead")
        return ak_data

def get_part(tree, entry_start=0, entry_stop=-1, chunk_size=1000, kprimary=1):
    """MC particles reading, return ak with calculated eta, momentum etc"""
    # print("Reading MC particles as akward arrays...")
    
    # Read as awkward array
    ak_data = get_branch_ak(tree, "MCParticles", entry_start, entry_stop, 
                                     chunk_size=chunk_size)
    # Preserve original per-event particle index so downstream code can
    # match podio particle.id().index even after filtering.
    orig_subentry = ak.local_index(ak_data.PDG, axis=1)

    if kprimary:
        print("Select all primary particles with generatorStatus==x001 or x002")
        primary_mask = (ak_data.generatorStatus % 1000 == 1) | (ak_data.generatorStatus % 1000 == 2)
        ak_data = ak_data[primary_mask]
        orig_subentry = orig_subentry[primary_mask]
    # Compute momentum quantities vectorized
    px = ak_data["momentum.x"]
    py = ak_data["momentum.y"]
    pz = ak_data["momentum.z"]
    
    p_mag = np.sqrt(px**2 + py**2 + pz**2)
    safe_p_mag = ak.where(p_mag != 0, p_mag, np.nan)
    theta = np.arccos(pz / safe_p_mag)
    phi = np.arctan2(py, px)
    eta = -np.log(np.tan(theta / 2.0))
    pt = p_mag * np.sin(theta)
    
    # Compute vertex quantities
    vx = ak_data["vertex.x"]
    vy = ak_data["vertex.y"]
    vz = ak_data["vertex.z"]
    vertex_r = np.sqrt(vx**2 + vy**2)
    vertex_dist = np.sqrt(vertex_r**2 + vz**2)
    
    # Compute endpoint quantities  
    ex = ak_data["endpoint.x"]
    ey = ak_data["endpoint.y"]
    endpoint_r = np.sqrt(ex**2 + ey**2)
    
    # Get PDG names (this is slower, so we do it efficiently)
    pdg_codes = ak.to_numpy(ak.flatten(ak_data.PDG))
    unique_pdgs = np.unique(pdg_codes)
    
    # Create PDG name mapping for unique values only
    pdg_name_map = {}
    for pdg in unique_pdgs:
        pdg_name_map[pdg] = get_pdg_info(pdg).name
    
    # Apply mapping vectorized
    pdg_names = ak.unflatten(
        np.array([pdg_name_map[pdg] for pdg in pdg_codes]),
        ak.num(ak_data.PDG)
    )
    
    # Create new awkward array with selected quantities
    my_field=['PDG', 'generatorStatus',  'charge', 'time', 'mass',
       'vertex.x', 'vertex.y', 'vertex.z', 'endpoint.x', 'endpoint.y',
       'endpoint.z', 'momentum.x', 'momentum.y', 'momentum.z']
    enhanced_data = ak.zip({
        # Original fields
        **{field: ak_data[field] for field in my_field},
        # Derived quantities
        'mom': p_mag,
        'theta': theta,
        'phi': phi,
        'eta': eta,
        'pt': pt,
        'vertex_r': vertex_r,
        'vertex_dist': vertex_dist,
        'endpoint_r': endpoint_r,
        'pdg_name': pdg_names,
        'orig_subentry': orig_subentry
    })
    
    return enhanced_data

def get_params(tree, bname="CentralCKFTrackParameters", entry_start=0, entry_stop=-1, chunk_size=1000):
    """Track Parameters reading, return ak with calculated eta, mom, pt"""
    
    # Read as awkward array
    ak_data = get_branch_ak(tree, bname, entry_start, entry_stop, 
                                     chunk_size=chunk_size)
    eta = theta2eta(ak_data.theta)
    mom = abs(1.0 / ak_data.qOverP)
    pt = abs(mom * np.sin(ak_data.theta))
    ak_data = ak.with_field(ak_data, eta, "eta")
    ak_data = ak.with_field(ak_data, mom, "mom")
    ak_data = ak.with_field(ak_data, pt, "pt")
    return ak_data

def get_branches(trees,bname):
    df=pd.DataFrame()
    for tree in trees:
        if bname=="MCParticles":
            dff = get_part(tree)
        else:
            dff=get_branch_df(tree,bname)
        df = pd.concat([df,dff],ignore_index=True)

        ## add a new counter (event_id)
        event_id = []
        current_id = -1
        prev_entry = None
        for e in df["entry"]:
            if prev_entry is None:
                current_id += 1
            elif e != prev_entry:
                current_id += 1
            event_id.append(current_id)
            prev_entry = e
        df["event_id"] = event_id

    return df


def get_collections(tree, bname='', kflatten=1):
    """Extract collections that a given branch pointed to (one to one/many relation)"""

    if not COL_TABLE:
        raise RuntimeError("COL_TABLE not populated. Call get_col_table() first.")    
    # Use the optimized awkward array reader
    br = get_branch_ak(tree, bname) if kflatten else get_branch_df(tree, bname, chunk_size=1000)
    
    # Convert to DataFrame to check for collectionID
    if hasattr(br, 'fields') and 'collectionID' in br.fields:
        # Get unique collection IDs efficiently
        colID = np.unique(ak.to_numpy(ak.flatten(br.collectionID)))
        collections = {}

        print(f"Loading {len(colID)} collections...")
        for ii in colID:
            if ii in COL_TABLE:
                # Use optimized reader for each collection
                collections[ii] = get_branch_ak(tree, COL_TABLE[ii]) if kflatten else get_branch_df(tree, COL_TABLE[ii], chunk_size=1000)
            else:
                print(f"Warning: Collection ID {ii} not found in COL_TABLE")
        
        return collections
    else:
        # Convert to DataFrame to check columns if it's awkward array
        if hasattr(br, 'fields'):
            br_df = ak.to_dataframe(br).reset_index()
        else:
            br_df = br
            
        if "collectionID" in br_df.columns:
            colID = br_df.collectionID.unique()
            collections = {}
            for ii in colID:
                if ii in COL_TABLE:
                    collections[ii] = get_branch_df(tree, COL_TABLE[ii], chunk_size=1000)
            return collections
        else:
            print("ERROR(get_collections):", bname, "is not a relation.")
            return 0

def get_relation(tree, b_name, v_name):
    """Get relation or vector members with index"""
    print(f"Processing relation: {b_name}.{v_name}")
    
    # Read main branch as DataFrame for index operations
    br = get_branch_df(tree, b_name, chunk_size=1000)
    
    # Check if the required columns exist
    begin_col = v_name + "_begin"
    end_col = v_name + "_end"
    
    if begin_col not in br.columns or end_col not in br.columns:
        print(f"ERROR(get_relation): {begin_col} or {end_col} not found in {b_name}")
        return 0
    
    loc1 = br.columns.get_loc(begin_col)
    loc2 = br.columns.get_loc(end_col)
    in_name = "_" + b_name + "_" + v_name
    
    # Read the relation/vector data
    try:
        app = get_branch_df(tree, in_name, chunk_size=1000)
    except:
        print(f"ERROR(get_relation): Cannot read branch {in_name}")
        return 0
    
    if not isinstance(app, pd.DataFrame):
        print(f"ERROR(get_relation): {in_name} is not a valid DataFrame")
        return 0
    
    # Vector of float values
    if len(app.columns) == 3 and app.columns[2] == "values":
        print("Processing vector values...")
        l_val = []
        l_ind = []
        # Vectorized approach for better performance
        for row in br.itertuples(index=False):
            l_ind.append(row[0])
            
            i1, i2 = row[loc1], row[loc2]
            if i1 == i2:  # empty
                l_val.append(np.array([]))
            else:
                # Use entry to locate the correct event
                event_data = app[app['entry'] == row.entry]
                if not event_data.empty:
                    v_row = np.array(event_data['values'])[i1:i2]
                    l_val.append(v_row)

                else:
                    l_val.append(np.array([]))
        
        return pd.DataFrame({"entry": l_ind, "values": l_val})
    
    # Relations with collectionID and index
    elif len(app.columns) > 1 and 'index' in app.columns and 'collectionID' in app.columns:
        print("Processing relations...")
        l_index = []
        l_id = []
        
        # Group app by entry for faster lookup
        app_grouped = app.groupby('entry')
        
        for row in br.itertuples(index=False):
            i1, i2 = row[loc1], row[loc2]
            if i1 == i2:  # empty
                l_index.append(np.array([]))
                l_id.append(np.array([]))
            else:
                # Get the event data
                if row.entry in app_grouped.groups:
                    v_row = app_grouped.get_group(row.entry)
                    if len(v_row) >= i2:  # Check bounds
                        l1 = np.array(v_row["index"])[i1:i2]
                        l2 = np.array(v_row["collectionID"])[i1:i2]
                        # l2 = np.array(v_row["collectionID"].iloc[0])[i1:i2] if hasattr(v_row["collectionID"].iloc[0], '__iter__') else np.array(v_row["collectionID"])[i1:i2]
                        
                        l_index.append(l1)
                        l_id.append(l2)
                    else:
                        l_index.append(np.array([]))
                        l_id.append(np.array([]))
                else:
                    l_index.append(np.array([]))
                    l_id.append(np.array([]))
        
        return pd.DataFrame({"index": l_index, "collectionID": l_id})
    
    else:
        print("ERROR(get_relation): Invalid vector or relation member structure")
        print(f"Columns found: {app.columns.tolist()}")
        return 0

def get_branch_relation(tree, branch_name="CentralCKFTrajectories", relation_name="measurements_deprecated", relation_variables=["*"]):
    """Get relation or vector members appended to the original branch with optimization"""
    print(f"Processing branch relation: {branch_name}.{relation_name}")
    if not COL_TABLE:
        raise RuntimeError("COL_TABLE not populated. Call get_col_table() first.")    
    # Read main branch
    br = get_branch_df(tree, branch_name, chunk_size=1000)
    df = get_relation(tree, branch_name, relation_name)
    
    if not isinstance(df, pd.DataFrame):
        print("ERROR (get_branch_relation): please provide a valid relation name.")
        return br, None
    
    # Handle vector values case
    if len(df.columns) == 1 and df.columns[0] == "values":
        return br, df["values"]
    
    # Handle relations case
    br[relation_name + "_index"] = df["index"]
    br[relation_name + "_colID"] = df["collectionID"]
    
    # Return early if no relation variables requested
    if len(relation_variables) == 0:
        br = br.explode([relation_name + "_index", relation_name + "_colID"]).reset_index(drop=True)
        return br, None
    
    # Prepare collections
    in_name = "_" + branch_name + "_" + relation_name
    print("Loading collections...")
    collections = get_collections(tree, in_name, 0)
    
    if not collections:
        print(f"ERROR: No collections {in_name} found")
        return br, None
    
    # Process relations efficiently
    print("Processing relation data...")
    l_relations = []
    l_collection_names = []
    loc1 = br.columns.get_loc(relation_name + "_index")
    loc2 = br.columns.get_loc(relation_name + "_colID")
    
    # Pre-convert collections to grouped DataFrames for faster access
    collections_grouped = {}
    sample_columns = None
    
    for col_id, col_data in collections.items():
        if hasattr(col_data, 'fields'):
            col_df = ak.to_dataframe(col_data).reset_index()
        else:
            col_df = col_data
        collections_grouped[col_id] = col_df.groupby('entry')
        if sample_columns is None:
            sample_columns = col_df.columns
    
    # Process each row
    for row in br.itertuples(index=False):
        ind = row[loc1]
        col = row[loc2]
        if len(ind) == 0:  # empty relation
            l_collection_names.append(None)
            # Add NaN row with correct number of columns
            if sample_columns is not None:
                l_relations.append([np.nan] * (len(sample_columns) - 2))  # -2 for entry, subentry
            else:
                l_relations.append([np.nan])
        else:
            for ii, cc in zip(ind, col):
                if cc in COL_TABLE:
                    l_collection_names.append(COL_TABLE[cc])
                    
                    # Get data from pre-grouped collections
                    if cc in collections_grouped and row.entry in collections_grouped[cc].groups:
                        event_data = collections_grouped[cc].get_group(row.entry)
                        if ii < len(event_data):
                            # Get the row data, excluding entry and subentry columns
                            row_data = event_data.iloc[ii]
                            filtered_data = [row_data[col] for col in row_data.index if col not in ['entry', 'subentry']]
                            l_relations.append(filtered_data)
                        else:
                            l_relations.append([np.nan] * (len(sample_columns) - 2))
                    else:
                        l_relations.append([np.nan] * (len(sample_columns) - 2))
                else:
                    l_collection_names.append(f"Unknown_{cc}")
                    l_relations.append([np.nan] * (len(sample_columns) - 2))
    
    # Explode the main branch
    br = br.explode([relation_name + "_index", relation_name + "_colID"]).reset_index(drop=True)
    
    # Create the additional DataFrame
    if sample_columns is not None:
        column_names = [col for col in sample_columns if col not in ['entry', 'subentry']]
    else:
        column_names = ['value']  # fallback
    
    df_add = pd.DataFrame(l_relations, columns=column_names)
    
    # Filter columns based on relation_variables
    if relation_variables != ["*"]:
        available_columns = select_string(column_names, relation_variables)
        df_add = df_add[available_columns]
    
    # Add metadata
    df_add[relation_name + "_colName"] = l_collection_names
    
    # Add entry and subentry columns at the beginning
    df_add.insert(0, 'entry', br['entry'])
    df_add.insert(1, 'subentry', br['subentry'])
    
    print(f"Completed processing. Result shape: {df_add.shape}")
    return br, df_add

def get_traj_relations(tree,bname="CentralCKFTrajectories",l_var=["measurementChi2", "outlierChi2", "trackParameters", "measurements_deprecated", "outliers_deprecated"]): 
    # l_var   = ["measurementChi2", "outlierChi2", "trackParameters", "measurements_deprecated", "outliers_deprecated"]#,"seed"]
    br      = get_branch_df(tree,bname)
    print("get_traj_relations: accessing the following vector members:")
    for vv in l_var:
        print(vv)
        a     = get_relation(tree,bname,vv)
        if not isinstance(a, pd.DataFrame):
            print(f"WARNING(get_traj_relations): {bname}.{vv} returned no data")
            continue
        for cc in a.columns:
            if cc=='values':
                br[vv]=a[cc]
                break
            elif cc=='index':
                br[vv+'_index']=a[cc]
            elif cc=='collectionID':
                br[vv+'_colID']=a[cc]
            else:
                print('WARNING: invalid column ',cc,' in ',bname,'_',vv)
    return br


def get_traj_hits_particle(tree, traj_name = 'CentralCKFTrajectories', measurement_name='measurements_deprecated'):
    """Trace trajectory measurements to reconstructed hits and MC particles.

    ``CentralTrackerMeasurements`` is a PODIO subset collection in current
    EICrecon output. Its ``*_objIdx`` branch points to the physical measurement
    collections, such as ``CentralWithoutTOFTrackerMeasurements`` and the TOF
    cluster measurements. Resolve those references before following each
    measurement's ``hits`` relation.

    A particle index of -1 means that the hit has no valid MC-particle
    reference, as expected for noise hits.
    """
    if not COL_TABLE:
        raise RuntimeError("COL_TABLE not populated. Call get_col_table() first.")

    tree_branches = set(tree.keys())
    measurement_collection = "CentralTrackerMeasurements"
    subset_branch = f"{measurement_collection}_objIdx"
    measurement_hit_frames = []

    if subset_branch in tree_branches:
        # Each subset entry stores the collection ID and object index of the
        # physical Measurement2D object. Discover the source collections from
        # those references instead of looking for a data branch named
        # CentralTrackerMeasurements.
        subset_references = get_branch_df(tree, subset_branch)
        source_collection_ids = pd.unique(subset_references["collectionID"])

        for raw_collection_id in source_collection_ids:
            if pd.isna(raw_collection_id):
                continue
            collection_id = int(raw_collection_id)
            if collection_id in (0, 0xFFFFFFFF):
                continue

            source_name = COL_TABLE.get(collection_id)
            if source_name is None:
                print(
                    f"WARNING: measurement collection ID {collection_id} "
                    "is not in COL_TABLE"
                )
                continue
            if source_name not in tree_branches or f"_{source_name}_hits" not in tree_branches:
                print(f"WARNING: {source_name} does not contain a readable hits relation")
                continue

            _, source_hits = get_branch_relation(
                tree,
                branch_name=source_name,
                relation_name="hits",
            )
            if source_hits is None:
                print(
                    f"WARNING: could not resolve the hits referenced by {source_name}. "
                    "The referenced reconstructed-hit collection may be missing "
                    "from the input file."
                )
                continue

            source_hits = source_hits.copy()
            source_hits["measurement_colID"] = collection_id
            source_hits["measurement_colName"] = source_name
            measurement_hit_frames.append(source_hits)
    else:
        # Backward compatibility for older files where
        # CentralTrackerMeasurements was stored as a normal collection.
        _, source_hits = get_branch_relation(
            tree,
            branch_name=measurement_collection,
            relation_name="hits",
        )
        if source_hits is not None:
            measurement_ids = [
                collection_id
                for collection_id, name in COL_TABLE.items()
                if name == measurement_collection
            ]
            source_hits = source_hits.copy()
            source_hits["measurement_colID"] = (
                int(measurement_ids[0]) if measurement_ids else -1
            )
            source_hits["measurement_colName"] = measurement_collection
            measurement_hit_frames.append(source_hits)

    if measurement_hit_frames:
        df_hits = pd.concat(measurement_hit_frames, ignore_index=True)
    else:
        df_hits = pd.DataFrame(
            columns=[
                "entry",
                "subentry",
                "hits_colName",
                "cellID",
                "measurement_colID",
                "measurement_colName",
            ]
        )

    # Match reconstructed tracker hits to simulated hits by event and cell ID.
    # This is exact for the current one-to-one tracker digitization. The shared
    # TOF reconstructed-hit names are aliases for the corresponding TOF sim
    # hit collections.
    name_sim_tracker = name_sim_barrel + name_sim_disk + ["B0TrackerHits"]
    name_rec_tracker = name_rec_barrel + name_rec_disk + ["B0TrackerRecHits"]
    rec_to_sim_hit = dict(zip(name_rec_tracker, name_sim_tracker))
    rec_to_sim_hit.update({
        "TOFBarrelSharedRecHits": "TOFBarrelHits",
        "TOFEndcapSharedRecHits": "TOFEndcapHits",
    })

    df_hits["particle_index"] = -1
    sim_particle_lookups = {}

    if not df_hits.empty and "cellID" in df_hits.columns:
        for rec_hit_name in pd.unique(df_hits["hits_colName"].dropna()):
            sim_hit_name = rec_to_sim_hit.get(rec_hit_name)
            if sim_hit_name is None:
                print(f"WARNING: {rec_hit_name} is not a recognized tracker hit collection")
                continue

            if sim_hit_name not in sim_particle_lookups:
                particle_relation_name = f"_{sim_hit_name}_particle"
                if (
                    sim_hit_name not in tree_branches
                    or particle_relation_name not in tree_branches
                ):
                    print(f"WARNING: {particle_relation_name} is not a branch")
                    sim_particle_lookups[sim_hit_name] = pd.DataFrame()
                    continue

                sim_hits = get_branch_df(tree, sim_hit_name)
                particle_relations = get_branch_df(tree, particle_relation_name)

                particle_indices = pd.to_numeric(
                    particle_relations["index"], errors="coerce"
                ).fillna(-1).astype(np.int64)
                particle_collection_names = particle_relations["collectionID"].map(
                    lambda collection_id: COL_TABLE.get(int(collection_id))
                    if pd.notna(collection_id)
                    else None
                )
                valid_particle_reference = (
                    particle_collection_names.eq("MCParticles")
                    & particle_indices.ge(0)
                    & particle_indices.ne(0xFFFFFFFF)
                )

                particle_relations = particle_relations[
                    ["entry", "subentry"]
                ].copy()
                particle_relations["particle_index"] = np.where(
                    valid_particle_reference,
                    particle_indices,
                    -1,
                )

                sim_hit_lookup = sim_hits[
                    ["entry", "subentry", "cellID"]
                ].merge(
                    particle_relations,
                    on=["entry", "subentry"],
                    how="left",
                    validate="one_to_one",
                )
                sim_hit_lookup["particle_index"] = (
                    sim_hit_lookup["particle_index"].fillna(-1).astype(np.int64)
                )
                sim_particle_lookups[sim_hit_name] = sim_hit_lookup.drop_duplicates(
                    ["entry", "cellID"], keep="first"
                )

            sim_hit_lookup = sim_particle_lookups[sim_hit_name]
            if sim_hit_lookup.empty:
                continue

            rec_hit_mask = df_hits["hits_colName"].eq(rec_hit_name)
            rec_hits_to_match = df_hits.loc[
                rec_hit_mask, ["entry", "cellID"]
            ].copy()
            rec_hits_to_match["_df_hits_index"] = rec_hits_to_match.index
            matched_particles = rec_hits_to_match.merge(
                sim_hit_lookup[["entry", "cellID", "particle_index"]],
                on=["entry", "cellID"],
                how="left",
            )
            df_hits.loc[
                matched_particles["_df_hits_index"].to_numpy(),
                "particle_index",
            ] = matched_particles["particle_index"].fillna(-1).to_numpy()

    # A trajectory relation stores the source measurement collection ID and
    # source object index. Match both values; measurement indices are only
    # unique within a collection and event.
    trajectory_relations, _ = get_branch_relation(
        tree,
        branch_name=traj_name,
        relation_name=measurement_name,
        relation_variables=[],
    )

    # get_branch_relation preserves a trajectory with an empty relation as one
    # row whose relation index and collection ID are NaN. Such a row represents
    # no hit (for example, a trajectory with nOutliers == 0), so do not label it
    # as a noise hit in the returned hit table. Keep non-empty but invalid PODIO
    # references; those are real relation entries and remain particle_index=-1.
    relation_index = f"{measurement_name}_index"
    relation_collection = f"{measurement_name}_colID"
    trajectory_relations = trajectory_relations[
        trajectory_relations[relation_index].notna()
        & trajectory_relations[relation_collection].notna()
    ].reset_index(drop=True)

    if df_hits.empty:
        traj_hits = trajectory_relations.copy()
        traj_hits["particle_index"] = -1
    else:
        hit_lookup = df_hits.drop_duplicates(
            ["entry", "measurement_colID", "subentry"],
            keep="first",
        ).rename(columns={"subentry": "_measurement_index"})

        traj_hits = trajectory_relations.merge(
            hit_lookup,
            left_on=[
                "entry",
                f"{measurement_name}_colID",
                f"{measurement_name}_index",
            ],
            right_on=["entry", "measurement_colID", "_measurement_index"],
            how="left",
            sort=False,
        )
        traj_hits.drop(columns=["_measurement_index"], inplace=True)
        traj_hits["particle_index"] = (
            traj_hits["particle_index"].fillna(-1).astype(np.int64)
        )

    if {"position.x", "position.y"}.issubset(traj_hits.columns):
        traj_hits["position.r"] = np.sqrt(
            traj_hits["position.x"]**2 + traj_hits["position.y"]**2
        )
        traj_hits["position.phi"] = np.arctan2(
            traj_hits["position.y"], traj_hits["position.x"]
        )
    else:
        traj_hits["position.r"] = np.nan
        traj_hits["position.phi"] = np.nan

    return traj_hits

from lmfit.models import GaussianModel

def gaussian(x, amplitude, mean, std):
    return amplitude * np.exp(-0.5 * ((x - mean) / std) ** 2) / (std * np.sqrt(2 * np.pi))

def hist_gaus(
    data, ax,
    bins=100, klog=False, header=None,
    clip=3.0, max_iters=5, min_points=50,
    verbose=False,
):
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    if data.size < min_points:
        if verbose:
            print('hist_gaus: not enough finite points')
        return np.nan, np.nan, np.nan

    center = np.median(data)
    mad = np.median(np.abs(data - center))
    scale = 1.4826 * mad if mad > 0 else np.std(data)
    if not np.isfinite(scale) or scale <= 0:
        if verbose:
            print('hist_gaus: invalid initial scale')
        return np.nan, np.nan, np.nan

    for _ in range(max_iters):
        mask = np.abs(data - center) <= clip * scale
        clipped = data[mask]
        if clipped.size < min_points:
            break
        new_center = np.median(clipped)
        mad = np.median(np.abs(clipped - new_center))
        new_scale = 1.4826 * mad if mad > 0 else np.std(clipped)
        if not np.isfinite(new_scale) or new_scale <= 0:
            break
        if np.isclose(new_center, center) and np.isclose(new_scale, scale):
            center, scale = new_center, new_scale
            break
        center, scale = new_center, new_scale

    lo, hi = center - clip * scale, center + clip * scale
    if lo == hi:
        if verbose:
            print('hist_gaus: degenerate range')
        return np.nan, np.nan, np.nan

    counts, edges = np.histogram(data, bins=bins, range=(lo, hi))
    mid = 0.5 * (edges[:-1] + edges[1:])

    if ax is not None:
        ax.hist(data, bins=bins, range=(lo, hi), histtype='stepfilled', alpha=0.3)

    model = GaussianModel()
    params = model.make_params(center=center, sigma=scale, amplitude=np.max(counts))
    try:
        result = model.fit(counts, params, x=mid)
    except Exception as exc:
        if verbose:
            print('hist_gaus: fit failed', exc)
        return np.nan, np.nan, np.nan

    sigma = float(result.params['sigma'])
    sigma_err = result.params['sigma'].stderr
    if sigma_err is None or not np.isfinite(sigma) or sigma <= 0:
        if verbose:
            print('hist_gaus: invalid fit result')
        return np.nan, np.nan, np.nan

    mean = float(result.params['center'])
    sigma_err = float(sigma_err)
    ampl = float(result.params['amplitude'])
    peak = ampl / (sigma * np.sqrt(2 * np.pi))

    if ax is not None:
        ax.plot(mid, gaussian(mid, ampl, mean, sigma), 'r-')
        if header:
            ax.set_title(header)
        ax.set_xlabel('value')
        ax.set_ylabel('entries')
        if klog:
            ax.set_yscale('log')
        else:
            ymax = max(np.max(counts), peak)
            ax.set_ylim(0, ymax / 0.7)

    return mean, sigma, sigma_err


__all__ = [
    "configure_analysis_environment",
    "ak_flat",
    "ak_df",
    "ak_hist",
    "ak_filter",
    "get_pdg_info",
    "get_geoID",
    "theta2eta",
    "select_string",
    "read_ur",
    "get_col_table",
    "get_branch_ak",
    "get_branch_df",
    "get_part",
    "get_params",
    "get_branches",
    "get_collections",
    "get_relation",
    "get_branch_relation",
    "get_traj_relations",
    "get_traj_hits_particle",
    "deg2rad",
    "status_to_source",
    "geo_mask_values",
    "hist_gaus"
]
