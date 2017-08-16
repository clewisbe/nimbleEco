Sys.setenv("R_TESTS" = "")  #Per Fritz suggestion

library(testthat)
library(nimbleEcology)

#For some reason devtools::test() works, but test_check() says there are no tests.
test_check("nimbleEcology")

