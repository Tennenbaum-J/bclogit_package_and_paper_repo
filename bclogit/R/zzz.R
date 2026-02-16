#' @importFrom methods new
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel CxxFlags
.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste("Welcome to bclogit v",
            utils::packageVersion("bclogit"),
            sep = ""
        )
    )
}
