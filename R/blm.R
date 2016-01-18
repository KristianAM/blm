#' Construct model matrix
#'
#' This functions constructs a model matrix with or without provided data
#'
#'@param formula an object of type formula
#'@param Data an object of type dataframe
#'@return a model matrix with data either from df or from the global environment
construct_model_matrix <- function(formula, Data = NULL){
  if(is.null(Data)) return(model.matrix(model.frame(formula)))
  else return(model.matrix(formula, model.frame(formula, Data)))
}
#'Construct responseless model matrix
#'
#'This function constructs a model matrix for data with no response
#'
#'@param formula an object of type formula
#'@param Data an object of type dataframe
#'@return a responseless model matrix
construct_responseless_model_matrix <- function(formula, Data){
  responseless <- delete.response(terms(formula))
  return(model.matrix(responseless, model.frame(responseless, data = Data)))
}

#' get posterior of blm object
#'
#' This function fetches the mean(s) and covariance matrix of the posterior distribution
#'
#' @param object an object of class blm
#' @return a list of mean(s) and covariance matrix
  posterior <-function(object) object$posterior

#'bayesian linear regression function
#'
#'main constructor function for creating an object of class blm
#'
#'@param formula an object of type formula
#'@param Data an object of type dataframe
#'@param prior an object of type blm
#'@param beta the precision
#'@return an object of class blm
#'
#'@export
blm <- function(formula, prior = NULL, beta = 1.3, Data = NULL){
  call <- match.call()
  if(is.null(prior)) {alpha = 0.5}
  else {alpha <- prior$variance}
  if(is.null(Data)){ model_matrix <-construct_model_matrix(formula)
                     response <- model.response(model.frame(formula))
                     model <- model.frame(formula)}
  else{
                     model_matrix <- construct_model_matrix(formula, Data = Data)
                     response <- model.response(model.frame(formula, Data))
                     model <- model.frame(formula, data = Data)}
  Sxy <- diag(alpha, nrow = ncol(model_matrix)) + beta * t(model_matrix) %*% model_matrix
  Sxy <- solve(Sxy)
  Mxy <- beta * Sxy %*% t(model_matrix) %*% response
  posterior <- list()
  posterior$variance <- Sxy
  posterior$mean <- Mxy

  coef <- t(posterior$mean)
  names(coef) <- colnames(model_matrix)

  object <- list()
  object$call <- call
  object$coefficients <- coef
  object$model <- model
  object$posterior <- posterior
  class(object) <- 'blm'
  object$fitted.values <- predict(object)
  object$residuals <- response - object$fitted.values
  object
}

#' predict new values from blm
#'
#'A function for predicting values of y from new data
#'
#'@param object an object of class blm
#'@param newdata an object of type dataframe
#'@return predicted y values from newdata
#'
#'@export
predict.blm <- function(object, newdata = NULL){
  if (!inherits(object, "blm")) {warning("calling predict.blm(<fake-blm-object>) ...")}
  if (is.null(newdata)){ Data <- object$model}
  else{ Data <- newdata}
  model_matrix <- construct_responseless_model_matrix(formula(formula(object$call)), Data = Data)
  posterior <- posterior(object)
  y <- t(posterior$mean) %*% t(model_matrix)
  y
}

#' Update a blm
#'
#' A function for refitting data with a new prior.
#'
#' @param object an object of class blm
#' @param Data an object of type dataframe
#' @return an object of class blm
#'
#' @export
update.blm <- function(object, Data = NULL){
  if (!inherits(object, "blm")) {warning("calling predict.blm(<fake-blm-object>) ...")}
  model <- blm(formula, prior = posterior(blm), beta = 1.3, Data = Data)
  model
}

#' coefficients
#' A function for extracting the coefficients of a blm
#'
#' @param object an object of class blm
#' @return a named vector of coefficients
#'
#' @export
coefficients.blm <-function(object){
  if (!inherits(object, "blm")) {warning("calling predict.blm(<fake-blm-object>) ...")}
  object$coefficients
}

#' fitted values
#' A function for retrieving the fitted values of the response
#'@param object an object of class blm
#'@param ... further arguments passed to or from other methods
#'@return a vector of fitted values
#'
#'@export
fitted.blm <- function(object, ...){
  if (!inherits(object, "blm")) {warning("calling predict.blm(<fake-blm-object>) ...")}
  object$fitted.values
}

#' residual function
#' A function for calculating residuals
#' @param object an object of class blm
#' @param ... further arguments passed to or from other methods
#' @return a vector of residuals
#'
#' @export
residuals.blm <- function(object, ...){
  object$residuals
}

#' Deviance function
#' A function for calculating the deviance between fitted and observed values
#'
#' @param object an object of class blm
#' @param ... further arguments passed to or from other methods
#' @return the sum of squares of the model
#'
#' @export
deviance.blm <- function(object, ...){
  sum(residuals(object)^2)
}

#' confidence interval function
#' A function for calculating confidence intervals of coefficients
#'
#'@param object an object of class blm
#'@param parm a specification of which parameters are to be given confidence intervals
#'@param level a double/float
#'@param ... further arguments passed to or from other methods
#'@return confidence intervals of coefficients
#'
#'@export
confint.blm <- function(object, parm, level = 0.95, ...){
  if (!inherits(object, "blm")) {warning("calling confint.blm(<fake-blm-object>) ...")}
  cf <- coefficients(object)
  pnames <- names(cf)
  if (missing(parm))
      parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- format(a,3)
  lowererror <- qnorm(a[1]) * (diag(posterior(object)$variance)/sqrt(length(object$model[,1])))
  uppererror <- qnorm(a[2]) * (diag(posterior(object)$variance)/sqrt(length(object$model[,1])))
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  error <- matrix(c(lowererror, uppererror), nrow = nrow(ci))
  ci[] <- cf[parm] + error
  ci
}
#' plot function for blm
#' A function for plotting residuals versus fitted y values
#'
#' @param x an object of class blm
#' @return a plot of residuals versus fitted values
#' @param ... further arguments passed to or from other methods
#'
#' @export
plot.blm <- function(x, ...){
  if (!inherits(object, "blm")) {warning("calling plot.blm(<fake-blm-object>) ...")}
  plot(x = x$fitted.values, y = x$residuals,
       ylim = c(min(x$residuals), max(x$residuals)), xlim = c(min(x$fitted.values), max(x$fitted.values)),
       main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
  abline(h = 0, lty = 3, col = "gray")
}

#' print function for blm
#' A function for printing the call and coefficients of a blm
#'
#' @param x an object of class blm
#' @param ... further arguments passed to or from other methods
#' @return the function call and coefficients of the object
#'
#' @export
print.blm <- function(x, ...){
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(coefficients(object))) {
    cat("Coefficients:\n")
    print.default(format(coefficients(object),digits = max(3L, getOption("digits") - 3L)), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}
