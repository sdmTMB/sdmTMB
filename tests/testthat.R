options(Matrix.warnDeprecatedCoerce = 2)
options(lifecycle_verbosity = "quiet")

library(testthat)
library(sdmTMB)

test_check("sdmTMB")
