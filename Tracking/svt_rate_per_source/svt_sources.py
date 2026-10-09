"""SVT geometry selections, PODIO reading, and rec-hit source assignment.

All coordinates are mm. These windows reproduce the parent workflow's SVT
selections; they are analysis cuts, not a geometry or sensor-layout service.
"""

from collections import Counter, defaultdict

import awkward as ak
import numpy as np


STATUS_TO_SOURCE = {
    1: "DIS", 2: "DIS",
    2001: "SR", 2002: "SR",
    3001: "Bremstrahlung", 3002: "Bremstrahlung",
    4001: "Coulomb", 4002: "Coulomb",
    5001: "Touschek", 5002: "Touschek",
    6001: "Proton beam gas", 6002: "Proton beam gas",
}
BACKGROUND_COMPONENTS = {"Bremstrahlung", "Coulomb", "Touschek", "Proton beam gas"}
SOURCE_ORDER = ["ALL", "DIS", "SR", "other background",
                "Bremstrahlung", "Coulomb", "Touschek", "Proton beam gas",
                "other", "unmatched", "invalid_relation"]

# (surface name, rec-hit collection, barrel?, radius or z window in mm).
SURFACES = [
    ("L0", "SiBarrelVertexRecHits", True, (30, 42)),
    ("L1", "SiBarrelVertexRecHits", True, (46, 60)),
    ("L2", "SiBarrelVertexRecHits", True, (115, 130)),
    ("L3", "SiBarrelTrackerRecHits", True, (250, 290)),
    ("L4", "SiBarrelTrackerRecHits", True, (400, 450)),
    ("E-Si Disk 4", "SiEndcapTrackerRecHits", False, (-1055, -1000)),
    ("E-Si Disk 3", "SiEndcapTrackerRecHits", False, (-860, -840)),
    ("E-Si Disk 2", "SiEndcapTrackerRecHits", False, (-660, -640)),
    ("E-Si Disk 1", "SiEndcapTrackerRecHits", False, (-460, -440)),
    ("E-Si Disk 0", "SiEndcapTrackerRecHits", False, (-260, -240)),
    ("H-Si Disk 0", "SiEndcapTrackerRecHits", False, (240, 260)),
    ("H-Si Disk 1", "SiEndcapTrackerRecHits", False, (440, 460)),
    ("H-Si Disk 2", "SiEndcapTrackerRecHits", False, (690, 710)),
    ("H-Si Disk 3", "SiEndcapTrackerRecHits", False, (940, 960)),
    ("H-Si Disk 4", "SiEndcapTrackerRecHits", False, (1150, 1250)),
]

# Explicit relation paths keep the EDM convention auditable.
RELATIONS = [
    {
        "rec": "SiEndcapTrackerRecHits",
        "rec_raw": "_SiEndcapTrackerRecHits_rawHit",
        "association": "SiEndcapTrackerRawHitAssociations",
        "association_raw": "_SiEndcapTrackerRawHitAssociations_rawHit",
        "association_sim": "_SiEndcapTrackerRawHitAssociations_simHit",
        "sim_particle": "_TrackerEndcapHits_particle",
    },
    {
        "rec": "SiBarrelVertexRecHits",
        "rec_raw": "_SiBarrelVertexRecHits_rawHit",
        "association": "SiBarrelVertexRawHitAssociations",
        "association_raw": "_SiBarrelVertexRawHitAssociations_rawHit",
        "association_sim": "_SiBarrelVertexRawHitAssociations_simHit",
        "sim_particle": "_VertexBarrelHits_particle",
    },
    {
        "rec": "SiBarrelTrackerRecHits",
        "rec_raw": "_SiBarrelTrackerRecHits_rawHit",
        "association": "SiBarrelRawHitAssociations",
        "association_raw": "_SiBarrelRawHitAssociations_rawHit",
        "association_sim": "_SiBarrelRawHitAssociations_simHit",
        "sim_particle": "_SiBarrelHits_particle",
    },
]


def read_branch(tree, name, start, stop, fields=None):
    """Read one event range and shorten split PODIO names (e.g. position.x).

    The outer workflow handles chunking. All required fields must be present;
    missing collections or incompatible schemas fail rather than yield zeros.
    """
    try:
        data = tree[name].array(library="ak", entry_start=start, entry_stop=stop)
    except KeyError as error:
        raise KeyError(f"Missing required PODIO collection/relation: {name}") from error
    renamed = {}
    for field in data.fields:
        short = field[len(name) + 1:] if field.startswith(name + ".") else field
        if "[" not in short and (fields is None or short in fields):
            renamed[short] = data[field]
    if fields is not None:
        missing = set(fields) - renamed.keys()
        if missing:
            raise ValueError(f"{name}: missing fields {sorted(missing)}")
    if not renamed:
        raise ValueError(f"{name}: expected a split PODIO record, got {data.type}")
    return ak.zip(renamed)


def resolve_status(index, statuses, begin, end, parents):
    """Follow the first parent of status-0 secondaries, with a 100-link limit.

    Returns (resolved status, traced secondary?, unresolved ancestry?). None
    denotes an invalid initial MC index; an unresolved secondary retains 0.
    Multiple parents are deliberately not combined, matching the parent code.
    """
    if index < 0 or index >= len(statuses):
        return None, False, False
    if statuses[index] != 0:
        return int(statuses[index]), False, False
    seen = set()
    for _ in range(100):
        if index in seen:
            break
        seen.add(index)
        first, last = int(begin[index]), int(end[index])
        if first < 0 or first >= last or last > len(parents):
            break
        index = int(parents[first])
        if index < 0 or index >= len(statuses):
            break
        if statuses[index] != 0:
            return int(statuses[index]), True, False
    return 0, True, True


def assign_sources(rec_raw, weights, association_raw, association_sim,
                   sim_particle, mc_particles, mc_parents):
    """Assign one exclusive source to each rec hit in one event.

    ObjectID indices are event-local within the named collections in RELATIONS.
    Choose the largest association weight; ties choose the first row. These
    assumptions follow the original analysis. Collection IDs are not remapped.
    Invalid selected relations remain visible instead of trying another source.
    """
    raw_indices = ak.to_numpy(association_raw["index"])
    sim_indices = ak.to_numpy(association_sim["index"])
    particle_indices = ak.to_numpy(sim_particle["index"])
    weights = ak.to_numpy(weights)
    statuses = ak.to_numpy(mc_particles["generatorStatus"])
    begin = ak.to_numpy(mc_particles["parents_begin"])
    end = ak.to_numpy(mc_particles["parents_end"])
    parents = ak.to_numpy(mc_parents["index"])
    if not (len(raw_indices) == len(sim_indices) == len(weights)):
        raise ValueError("RawHitAssociation weights and relation lengths disagree")
    if not (len(statuses) == len(begin) == len(end)):
        raise ValueError("MCParticle status and parent-range lengths disagree")

    raw_to_rows = defaultdict(list)
    for row, raw_index in enumerate(raw_indices):
        if raw_index >= 0:
            raw_to_rows[int(raw_index)].append(row)
    # Cache ancestry results for particles shared by many associations/hits.
    resolved = {}

    def association_status(row):
        sim_index = int(sim_indices[row])
        if sim_index < 0 or sim_index >= len(particle_indices):
            return None, False, False
        particle_index = int(particle_indices[sim_index])
        if particle_index not in resolved:
            resolved[particle_index] = resolve_status(particle_index, statuses, begin, end, parents)
        return resolved[particle_index]

    labels = []
    diagnostics = Counter()
    for raw_index in ak.to_numpy(rec_raw["index"]):
        rows = raw_to_rows.get(int(raw_index), [])
        if raw_index < 0:
            labels.append("invalid_relation")
            diagnostics["invalid_relation_count"] += 1
            continue
        if not rows:
            labels.append("unmatched")
            diagnostics["unmatched_count"] += 1
            continue
        if len(rows) > 1:
            diagnostics["ambiguous_association_count"] += 1
            candidate_sources = {
                STATUS_TO_SOURCE.get(status, "other")
                for status, _, _ in (association_status(row) for row in rows)
                if status is not None
            }
            diagnostics["ambiguous_source_count"] += int(len(candidate_sources) > 1)
        if not np.all(np.isfinite(weights[rows])) or np.any(weights[rows] < 0):
            labels.append("invalid_relation")
            diagnostics["invalid_weight_count"] += 1
            diagnostics["invalid_relation_count"] += 1
            continue
        best = max(rows, key=lambda row: float(weights[row]))
        status, traced, unresolved = association_status(best)
        if status is None:
            labels.append("invalid_relation")
            diagnostics["invalid_relation_count"] += 1
            continue
        if traced:
            diagnostics["secondary_unresolved_count" if unresolved else "secondary_traced_count"] += 1
        labels.append(STATUS_TO_SOURCE.get(status, "other"))
    diagnostics["rec_hits"] += len(labels)
    return np.asarray(labels, dtype=str), diagnostics
