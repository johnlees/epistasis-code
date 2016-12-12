# epistasis-code
pairwise bacterial/human variant association

Uses the chi^2 and logistic regressions as implemented in seer to look for epistasis between two matched populations

## QC
Will check MAF and missing rate but that's it.

You should:
* Make sure the sample line order in all three input files (bacteria, human,
  covariates) matches.
* Make sure bacteria is coded 0 or 1, human is coded 0/0, 0/1 or 1/1.
  (missing . or ./.)
