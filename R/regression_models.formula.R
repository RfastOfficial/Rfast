#' Fit Linear Regression Model
#'
#' This function fits a linear regression model using either numeric data or a formula interface.
#' 
#' @param x matrix of predictor variables when using the default method.
#' @param y vector of response variable when using the default method.
#' @param form a formula specifying the model when using the formula method.
#' @param data an optional data frame in which to evaluate the variables in the formula.
#' @param w optional vector of weights.
#' 
#' @return A list with components:
#'   \item{be}{coefficients of the linear regression.}
#'   \item{residuals}{residuals of the linear regression.}
#' 
#' @examples
#' # Example with default method
#' x <- matrix(rnorm(100), ncol = 2)
#' y <- rnorm(50)
#' result_default <- lmfit(x, y)
#' 
#' # Example with formula method
#' data <- data.frame(y = rnorm(50), x1 = rnorm(50), x2 = rnorm(50))
#' result_formula <- lmfit(y ~ x1 + x2, data = data)
#' 
#' @export
#' @rdname lmfit
lmfit <- function(x, ...) {
  UseMethod("lmfit")
}

#' @export
#' @rdname lmfit
lmfit.default <- function(x, y, ...) {
  return(Rfast::lmfit(x = x, y = y, ...))
}

#' @export
#' @rdname lmfit
lmfit.formula <- function(form, data, w = NULL, ...) {
  form <- as.formula(form)
  vars <- all.vars(form)
  
  x_vars <- vars[-1]
  y_vars <- vars[1]
  
  if (missing(data)) {
    mm <- model.frame(form)
    X <- as.matrix(mm[, x_vars])
    y <- mm[, y_vars]
  } else {
    mm <- model.frame(form, data = data)
    X <- as.matrix(mm[, -1])
    y <- model.response(mm)
  }
  
  if (is.null(w)) {
    lm_fit <- Rfast::lmfit(x = X, y = y, ...)
  } else {
    lm_fit <- Rfast::lmfit(x = X, y = y, w = w, ...)
  }
  
  return(lm_fit)
}
