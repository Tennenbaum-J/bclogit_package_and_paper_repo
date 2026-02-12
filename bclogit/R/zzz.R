.onLoad <- function(libname, pkgname) {
    # Initialize globals in the package namespace
    assign("bclogit_globals", new.env(), envir = parent.env(environment()))
}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste("Welcome to bclogit v",
            utils::packageVersion("bclogit"),
            sep = ""
        )
    )
}
