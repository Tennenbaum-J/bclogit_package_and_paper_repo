#' bclogit
#'
#' Fits a conditional logistic regression model to incidence data. Allows for use of the concordant pairs in the fitting.
#'
#' @name 		bclogit-package
#' @title Bayesian Conditional Logistic Regression with Concordant Pairs
#' @author 		Jacob Tennenbaum \email{Jacob.Tennenbaum51@@qmail.cuny.edu}
#' @author 		Adam Kapelner \email{kapelner@@qc.cuny.edu}
#' @references 	Jacob Tennenbaum and Adam Kapelner (2026). "Improved Conditional Logistic Regression using Information in Concordant Pairs with Software." arXiv preprint arXiv:2602.08212.
#' @keywords 	Conditional Logistic Regression
#' @import      checkmate
#' @import      Rcpp
#' @importFrom  Rcpp evalCpp
#' @useDynLib   bclogit, .registration=TRUE
##### Run "library(roxygen2); roxygenise("bclogit", clean = TRUE)" to regenerate all Rd files and NAMESPACE and DESCRIPTION file
"_PACKAGE"
