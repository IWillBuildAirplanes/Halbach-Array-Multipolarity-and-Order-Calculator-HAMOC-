# HAMOC
Halbach Array Multipolarity and Order Calculator
MatLab software for designing Halbach arrays out of permanent magnet cubes and for determining the multipole elements of given arrays up to a certain specificed order.
This code is an edit of the previous code by Volker Ziemann called Permanent Magnet Design (PMD). Reading the manuscript I wrote may be helpful in deciphering these files. The only files that were added on top of PMD were
## AntiStrengthsOf.m
Gives strength of anti-array as a function of angle.
## MAIN_DegCoeffs.m
Outputs the coefficients for a Halbach array under a specified amount of degradation k in [0, 1].
## MAIN_DegCoeffsDiQuad.m
Outputs the coefficients for a Halbach array with nominal dipole and quadrupole elements under a specified amount of degradation k in [0, 1].
## MAIN_FigureCycle.m
Generates a figure of the element argument (direction) for a Halbach array.
## MAIN_Figures.m
Generates figures used in the paper.
## MAIN_FiguresDiQuad.m
Generates figures used in the paper for the dual dipole-quadrupole array.
## MAIN_VSTime.m
Generates figure of element degradation over time.
## MAIN_VSTimeDiQuad.m
Generates figure of element degradation over time for the dual dipole-quadrupole array.
## MPoleElemOf.m
Retrieves the coefficients from a group of sheets created from PMD.
## MakeArray
Makes a Halbach array with specified characteristics.
