# snippets/Exclusive

Snippets related to exclusive, diffractive, and tagged processes at the EIC.

## `ePIC_ExcKinUtils.h`

This is a header file which contains functions for calculating exclusive kinematic quantities. In order to use these functions within any ePIC analysis, the following line should be added to the relevant scripts:

```
#include "snippets/Exclusive/ePIC_ExcKinUtils.h"
```

The header file here contains functions for calculating the following quantities:
- The energy of a particle, from its scalar mass and vector momentum.
- Missing mass (squared), momentum, pT and energy, for a 2-body final state (a + b -> c + d)
- Missing mass (squared), momentum, pT and energy, for a 3-body final state (a + b -> c + d + f)
- The Mandelstam variable t, using multiple methods as defined in the tRECO convention note  **STILL IN PROGRESS**


The functions take as input the 4 vectors for particles identified within an ePIC event. These 4-vectors can take the form of any type which includes the relevant operators: `.E()`, `.P()`, `.Pt()` and `.M2()`. ROOT's legacy `TLorentzVector` class, and its more recent replacement, `ROOT::Math::LorentzVector`, both contain all of the operators required.
