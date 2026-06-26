# Jet Macro For Validation

**Last update:** 06.25.2026, DMA

### Shell Used
--------------

Utilized eic-shell version `26.05.0-stable` for testing.

### Files Used
--------------

Performed two tests: one before we switched the jet reconstruction output to
`edm4eic::Jet` (26.05.0) and one after (26.06.0).

#### 26.05.0

Utilized 26.05.0 Pythia8 NC DIS 10x100 (q2 = 100 - 1K) dataset for before test:

```bash
rucio did content list epic:/RECO/26.05.0/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/10x100/q2_100to1000
```

Complete filelist saved to `filelists/files26050.py8ncdis10x100q100t1000.list`. Example
plots generated using file `*_run000.0000.*` and saved to `examples/26050_000_0000`.

#### 26.06.0

Utilized 26.06.0 Pythia8 NC DIS 10x100 (q2 = 100 - 1K) dataset for after test:

```bash
rucio did content list epic:/RECO/26.06.0/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/10x100/q2_100to1000
```

Complete filelist saved to `filelists/files26060.py8ncdis10x100q100t1000.list`. Example
Plots generated using file `*_run000.0000.*` and saved to `examples/26060_000_000`.

### Plot File Names and Titles
------------------------------

These tables list the plot names and their titles which might be be best for Hydra.
Some (such as the resolution fit summaries) likely won't render well on the web. The
order of entries reflects the order which might be good to display on the web.

**Note:** the area plots won't appear if you're running on a campaign prior to 26.06.0.

#### Reconstructed Jets

| File Name | Title |
| --------- | ----- |
| numberRecoJets.png | Number of Charged Jets per Event (Reconstruction Level) |
| numberConstituentPerRecoJet.png | Number of Constituents per Charged Jet (Reconstruction Level) |
| recoJetEnergy.png | Charged Jet Energy ($|\eta| < 2.5$, Reconstruction Level) |
| recoJetEta.png | Charged Jet $\eta$ ($E_{jet}$ > 5 GeV, Reconstruction Level) |
| recoJetEnergyVsArea.png | Charged Jet Energy vs. Area (Reconstruction Level) |
| recoJetEnergyVsEta.png | Charged Jet Energy vs. $\eta$ (With $e^{-}$, Reconstruction Level) |
| recoJetEnergyVsEtaNoElectron.png | Charged Jet Energy vs. $\eta$ (No $e^{-}$, Reconstruction Level) |
| recoJetEnergyVsArea.png | Charged Jet Energy vs. Area (Reconstruction Level) |
| recoJetPhiVsEta.png | Charged Jet $\varphi$ vs. $\eta$ (With $e^{-}$, $E_{jet}$ > 5 GeV, Reconstruction Level) |
| recoJetPhiVsEtaNoElectron.png | Charged Jet $\varphi$ vs. $\eta$ (No $e^{-}$, $E_{jet}$ > 5 GeV, Reconstruction Level) |
| recoJetConstituent.png | Charged Constituent $\eta$ (Reconstruction Level) |
| recoJetConstituentMomentum.png | Charged Constituent Momenta (Reconstruction Level) |
| recoJetConstituentMomentumVsEta.png | Charged Constituent Momenta vs. $\eta$ (With $e^{-}$, Reconstruction Level) |
| recoJetConstituentMomentumVsEtaNoElectron.png | Charged Constituent Momenta vs. $\eta$ (No $e^{-}$, Reconstruction Level) |
| recoJetConstituentPairwiseDR.png | Pairwise Charged Constituent $\Delta R$ (Reconstruction Level) |
| recoJetConstituentPhiVsEta.png | Jet Constituent $\varphi$ vs. $\eta$ (With $e^{-}$, Reconstruction Level) |
| recoJetConstituentPhiVsEtaNoElectron.png | Jet Constituent $\varphi$ vs. $\eta$ (No $e^{-}$, Reconstruction Level) |

#### Generated Jets

| File Name | Title |
| --------- | ----- |
| numberGenJets.png | Number of Charged Jets per Event (Generator Level) |
| numberConstituentPerGenJet.png | Number of Constituents per Charged Jet (Generator Level) |
| genJetEnergy.png | Charged Jet Energy ($|\eta| < 2.5$, Generator Level) |
| genJetEta.png | Charged Jet $\eta$ ($E_{jet}$ > 5 GeV, Generator Level) |
| genJetAea.png | Charged Jet Area ($E_{jet}$ > 5 GeV, Generator Level) |
| genJetEnergyVsEta.png | Charged Jet Energy vs. $\eta$ (With $e^{-}$, Generator Level) |
| genJetEnergyVsEtaNoElectron.png | Charged Jet Energy vs. $\eta$ (No $e^{-}$, Generator Level) |
| genJetPhiVsEta.png | Charged Jet $\varphi$ vs. $\eta$ (With $e^{-}$, $E_{jet}$ > 5 GeV, Generator Level) |
| genJetPhiVsEtaNoElectron.png | Charged Jet $\varphi$ vs. $\eta$ (No $e^{-}$, $E_{jet}$ > 5 GeV, Generator Level) |
| genJetEnergyVsArea.png | Charged Jet Energy vs. Area (Generator Level) |
| genJetConstituent.png | Charged Constituent $\eta$ (Generator Level) |
| genJetConstituentMomentum.png | Charged Constituent Momenta (Generator Level) |
| genJetConstituentMomentumVsEta.png | Charged Constituent Momenta vs. $\eta$ (With $e^{-}$, Generator Level) |
| genJetConstituentMomentumVsEtaNoElectron.png | Charged Constituent Momenta vs. $\eta$ (No $e^{-}$, Generator Level) |
| genJetConstituentPairwiseDR.png | Pairwise Charged Constituent $\Delta R$ (Generator Level) |
| genJetConstituentPhiVsEta.png | Jet Constituent $\varphi$ vs. $\eta$ (With $e^{-}$, Generator Level) |
| genJetConstituentPhiVsEtaNoElectron.png | Jet Constituent $\varphi$ vs. $\eta$ (No $e^{-}$, Generator Level) |

#### Generated vs. Reconstructed Jets

| File Name | Title |
| --------- | ----- |
| genRecoJetDeltaR.png | $\Delta R$ Between Generated and Reconstructed Charged Jet Matches |
| matchedJetScaleResolutionSummary.png | Charged Jet Energy Resolution and Scale vs. $\eta$ |
| matchedJetResolutionVsEnergy.png | Charged Jet Energy Resolution vs. Generated Jet Energy (All $\eta$) |
| matchedJetResolutionVsEta.png | Charged Jet Energy Resolution vs. Generated Jet $\eta$ |
| matchedJetResolutionVsEnergyMidEta.png | Central Charged Jet Energy Resolution vs. Generated Jet Energy ($|\eta|$ < 1) |
| matchedJetResolutionVsEnergyNegEta.png | Backward Charged Jet Energy Resolution vs. Generated Jet Energy ($\eta \in (-2,5, -1.0)) |
| matchedJetResolutionVsEnergyPosEta.png | Forward Charged Jet Energy Resolution vs. Generated Jet Energy ($\eta \in (1.0, 2.5)) |
| matchedRecoVsGenJetEnergy.png | Reconstructed vs. Generated Charged Jet Energy |
| matchedRecoVsGenJetEta.png | Reconstructed vs. Generated Charged Jet $\eta$ |
| matchedRecoVsGenJetPhi.png | Reconstructed vs. Generated Charged jet $\varphi$ |
| matchedRecoVsGenJetArea.png | Reconstructed vs. Generated Charged Jet Area |
