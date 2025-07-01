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
- The Mandelstam variable t, using multiple methods as defined in the tRECO convention note


The functions take as input the 4 vectors for particles identified within an ePIC event. These 4-vectors can take the form of any type which includes the relevant operators: `.E()`, `.P()`, `.Pt()` and `.M2()`. ROOT's legacy `TLorentzVector` class, and its more recent replacement, `ROOT::Math::LorentzVector`, both contain all of the operators required.

## `tEXBABE.h`

A header file to calculate teXBABE. This is a generic formulation of a calculation used in the analysis of DEMP events. Vectors and inputs have been named generically to match the tRECO convention for ease of use and readability.

```
#include "snippets/Exclusive/teXBABE.h"
```

The header file contains three methods for calculating teXBABE -

- An explicit method which takes 6 inputs
  - Detected scattered electron vector
  - Detected X vector
  - Detected baryon (BA) vector
  - Mass of the expected baryon (BA)
  - Electron beam vector
  - Hadron beam vector
- A simplified method with 4 inputs
  - PMiss vector
  - Detected BA vector
  - Mass of the expected baryon (BA)
  - Hadron beam vector
- Shortest method with 2 inputs
  - Corrected BA vector
  - Hadron beam vector

Some comments and additional explanations are included in the file.

## Notes

A test of equivalence between ePIC_ExcKinUtils.h and teXBABE.h should be conducted.