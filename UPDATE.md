# 1.2.1 - 1.2.3

May 15-25, 2023

Add episol to call eprism3d/gmxtop2solute/gensolvent/ts4sdump

episol and eprism3d can directly read top files, no need to translate with gmxtop2solute

Fixed a bug of residue number counting in eprism3d


# 1.1.322 - 1.1.326

1.1.326: Miscellaneous bug fix

1.1.325 : Dec 10, 2022 : Minor bug fix

1.1.324: Nov 8, 2022 : Change: the cutoff of PLHNC is now set in closure-factor command (-cmd cf=...)

1.1.322 & 1.1.323 : Nov 3-4, 2022 : Change the distribution name to EPRISM3D


# 1.1.321

Nov 1, 2022

Add: closure: PSE, allows fractional order of PSE

e.g.: “cf=3 closure=pse” = PSE3; “cf=3.5 closure=PSE3” = 0.5 PSE3 + 0.5 PSE4

Add: free energy functionals for: PY, D2, MS, BPGG


# 1.1.320

Oct 18, 2022

Add: -allow-original-err, display original (unscaled) SCF error in RISM-closure interations

Add from from 0.234: -cr, detect and remove unrealistic water in cavities or pockets

Fix: rewrote MS and BPGG closures


# 1.0.316

June 3, 2022

The license changed to LGPL3 from GPL3.
