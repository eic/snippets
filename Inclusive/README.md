# snippets/Inclusive

Snippets related to inclusive processes at the EIC.

## `ePIC_IncKinUtils.h`

This is a header file which contains functions for calculating inclusive kinematic quantities by hand. In order to use these functions within any ePIC analysis, the following line should be added to the relevant scripts:

```
#include "snippets/Inclusive/ePIC_IncKinUtils.h"
```

The header file here duplicates the effects of the `InclusiveKinematics` branches within the output trees created by `EICrecon`, for those who wish to calculate Q2, x and y by hand. Functions are included for the following:
- Caclulating Q2, x and y separately.
- Calculating all three quantities simultaneously, and storing them in holding variables passed to the function.
- All of the above for all methods used for the calculation of the `InclusiveKinematics` branches (Electron, Jaquet-Blondel, Double Angle, Sigma and e-Sigma).

The functions take as input the 4 vectors for particles identified within an ePIC event. These 4-vectors can take the form of any type which includes the relevant operators: `.X()`, `.Y()`, `.Z()`, `.E()`, `.Pt()`, `.M2()` and `.Dot()`. ROOT's legacy `TLorentzVector` class, and its more recent replacement, `ROOT::Math::LorentzVector`, both contain all of the operators required.
