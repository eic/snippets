# Recover generator particles from npsim output
Shujie Li with Chatgpt, 10.2026

[`extract_generator_hepmc.py`](extract_generator_hepmc.py) converts the generator
record saved in MCparticles either from EDM4hep or EDM4eic file into HepMC3 ASCII, suitable for use
as input to another npsim run. It is a standalone script: no local imports,
PODIO Python bindings, or detector geometry are needed.

## Requirements and use

Python 3 with `numpy`, `awkward`, `uproot`, and `pyhepmc`:

```bash
python -m pip install numpy awkward uproot pyhepmc
python extract_generator_hepmc.py sample.edm4hep.root source.hepmc3 --max-events 2
```

Validated with `pyhepmc 2.14.0` and `uproot 5.5.1`.

Omit `--max-events` to process every event. The output directory must already
exist. Existing output files are never overwritten. During conversion the
script writes `source.hepmc3.partial`, then renames it on success. If conversion
fails, the partial file stays available for inspection; choose a new output
name or move it aside before retrying.

## What is preserved

- Keep MCParticles with **simulatorStatus bit 30 clear**
  (`BITCreatedInSimulation`); discard particles created by detector simulation.
  This keeps all saved generator statuses, including incoming and intermediate
  particles.
- Preserve particle order, PDG IDs, generator statuses, masses, event numbers,
  and parent/daughter connections within the retained record.
- Recompute energy as `E = sqrt(px² + py² + pz² + m²)`. Momentum, energy, and
  mass use GeV; positions use mm. Convert EDM4hep time in ns to HepMC `ct` in
  mm using `ct = 299.792458 × time`.

The input must contain the standard `events` tree, `EventHeader`, `MCParticles`,
and its parent/daughter relation branches (listed in `BRANCHES` in the script).
Each relation index is interpreted within that event's `MCParticles` collection.
Events are read individually with the array cache disabled, so memory depends
on event size rather than the number of events in the file.

## Reconstruction limits

This reconstructs the saved particle record, not the original HepMC file byte
for byte. Generator weights, run metadata, and unsaved attributes are not
recovered. Children with the same parent set share a vertex at the parents'
common endpoint; the nearest child's stored production position supplies the
time. Parentless particles must start at `(x, y, z, ct) = (0, 0, 0, 0)`.
Inconsistent relations, incompatible parent endpoints, and unsupported source
positions raise errors rather than silently changing the record.
