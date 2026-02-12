#' Extract model formula
#'
#' @param x A `bclogit` object.
#' @param ... Additional arguments.
#' @return The formula used in the model.
#' @export
formula.bclogit <- function(x, ...) {
    if (!is.null(x$terms)) {
        formula(x$terms)
    } else {
        NULL
    }
}
