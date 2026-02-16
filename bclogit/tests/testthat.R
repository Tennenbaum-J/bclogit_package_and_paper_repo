library(testthat)
library(bclogit)

envs <- Sys.getenv()
check_envs <- envs[grep("^_R_CHECK", names(envs))]
print(check_envs)

if (Sys.getenv("_R_CHECK_CRAN_INCOMING_") != "true") {
  test_check("bclogit")
}
