# Environment to store compiled Stan models
.bclogit_cache <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
    # Initialize cache if needed, though new.env above does it.
    # We can pre-compile here if we wanted, but lazy is better.
}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste("Welcome to bclogit v",
            utils::packageVersion("bclogit"),
            sep = ""
        )
    )
}
