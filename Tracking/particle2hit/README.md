# MC particles to tracking hits

Shujie Li with ChatGPT, 10.2026

Read an EDM4hep simulation ROOT file into linked pandas tables of MC particles
and their directly associated **simulated tracking hits**. The reader discovers
tracker-hit collections and CellID layouts from the file; ROOT, PODIO and DD4hep
Python bindings are not required.

* This script is designed for tracker acceptance study. A similar script to extract
 Reconstructed hits on each surface and group by source particles are available
  at `snippets/Tracking/svt_rate_per_source`.

## Getting started

Requires Python 3.9+, NumPy, Awkward Array, Uproot and pandas. The notebook also
uses Matplotlib and IPython/Jupyter. Tested in the EIC 26.08 environment.

Open [mc_tracking_hits_example.ipynb](mc_tracking_hits_example.ipynb) from this
directory, set `ROOT_FILE` and the MC cuts in section 1, and run all cells with
a kernel containing these dependencies. The notebook demonstrates layer
selection, complete hit attributes, per-particle hit counts and diagnostic plots.

## Python usage

```python
from mc_tracking_hits import read_mc_tracking_hits

tables = read_mc_tracking_hits(
    "simulation.root", pdgs=[211, -211],
    eta_range=(0.9, 1.0), p_range=(4.5, 5.5), stop=2000,
)
particles, hits = tables["particles"], tables["hits"]

# VertexBarrel L0 in the example geometry; verify IDs for your geometry.
l0_hits = hits.loc[
    (hits.hit_collection == "VertexBarrelHits")
    & (hits.system == 31) & (hits.layer == 0)
]
```

| Result | Contents |
| --- | --- |
| `particles` | Selected MC particles, including zero-hit particles, with kinematics and total hit counts |
| `hits` | Associated hits with all stored scalar attributes, original references, exact uint64 CellIDs and decoded fields |
| `surfaces` | Observed collection/system/layer groups with hit and distinct-particle counts |
| `metadata` | Source identity, selections, collection IDs, CellID layouts and units |

Join particles and hits on `file_id`, `event_entry`, `mc_index`. Select a layer
by **collection, system and layer**, since layer numbers are reused across
detectors. These groups are detector layers, not individual ACTS surfaces;
decoded module/sensor fields retain finer segmentation. Missing encoded fields
are represented by NA. Positions/path lengths are mm, momentum/energy GeV,
time ns and angles radians; decoded CellID fields are integer identifiers.

By default, `primary_only=True` selects `generatorStatus == 1`; multiple MC
particles per event are supported. Set it to `False` to include other stored MC
particles. Eta and momentum cuts use MC kinematics and half-open `[low, high)`
ranges. Only direct hit-to-particle links are followed, without secondary
ancestry. Hit presence does not imply track reconstruction or a geometric
intersection test.

## Large files

`read_mc_tracking_hits()` collects one file/prefix in memory. For batch analysis,
reduce each block instead of retaining all raw hit tables:

```python
from mc_tracking_hits import iter_mc_tracking_hits

for tables in iter_mc_tracking_hits("simulation.root", pdgs=[211], step_size=5000):
    hits = tables["hits"]
    # Select layers and accumulate counts or histograms here.
```

Use `collections=[...]` to restrict I/O. If the file lacks CellID-layout metadata,
pass `encodings` from the matching geometry; `load_readout_encodings(xml_files)`
reads readout descriptions from explicit detector XML files without expanding
includes. Input ROOT files are opened read-only.
