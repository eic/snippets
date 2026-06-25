# Jet Macro For Validation

**Last update:** 06.25.2026, DMA

### Shell Used

Utilized eic-shell version `26.05.0-stable` for testing.

### Files Used

Performed two tests: one before we switched the jet reconstruction output to
`edm4eic::Jet` (26.05.0) and one after (26.06.0).

#### 26.05.0
------------

Utilized 26.05.0 Pythia8 NC DIS 10x100 (q2 = 100 - 1K) dataset for before test:

```bash
rucio did content list epic:/RECO/26.05.0/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/10x100/q2_100to1000
```

Complete filelist saved to `filelists/files26050.py8ncdis10x100q100t1000.list`. Example
plots generated using file `*_run000.0000.*` and saved to `examples/26050_000_0000`.

#### 26.06.0
------------

Utilized 26.06.0 Pythia8 NC DIS 10x100 (q2 = 100 - 1K) dataset for after test:

```bash
rucio did content list epic:/RECO/26.06.0/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/10x100/q2_100to1000
```

Complete filelist saved to `filelists/files26060.py8ncdis10x100q100t1000.list`. Example
Plots generated using file `*_run000.0000.*` and saved to `examples/26060_000_000`.
