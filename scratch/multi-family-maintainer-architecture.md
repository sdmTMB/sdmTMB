# Multi-family maintainer architecture

Multi-family support uses one row-wise family table and at most two shared
latent linear predictors. It is not a separate model backend: each observation
selects a family, and that family selects one or both latent predictors.

```text
user family input
  -> family-spec compiler
  -> analysis-row and response preparation
  -> TMB family metadata
  -> C++ family resolver
  -> likelihood, simulation, prediction, and index reports
```

The canonical metadata is `family_spec`; do not route internals through
`object$family`. That public field preserves the complete user input, including
the named family list in multi-family fits.

`R/family-spec.R` compiles, validates, and serializes that metadata for TMB.
`R/family-response.R` prepares row-aligned responses. `R/family-prediction.R`
maps prediction rows to families, formats component outputs, and combines
realized simulation draws. TMB remains authoritative for combined prediction
reports.

To add a supported family safely:

1. Update the central registry in `R/family-spec.R`.
2. Keep its R and C++ enum codes aligned.
3. Declare any family-level auxiliary parameter in the registry.
4. Add ordinary-fit and multi-family equivalence tests.
5. Do not add row-wise family logic to downstream methods. Use `family_spec`,
   or add an early capability guard when the method is unsupported.
