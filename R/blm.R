#' Construct model matrix
#'
#' This functions constructs a model matrix with or without provided data
#'
#'@param formula an object of type formula
#'@param ... can be an object of type dataframe
#'@return a model matrix with data either from df or from the global environment
construct_model_matrix <- function(formula, ...){
  if(!(hasArg(df))) return(model.matrix(model.frame(formula)))
  else return(model.matrix(formula, model.frame(formula, df)))
}
#'Construct responseless model matrix
#'
#'This function constructs a model matrix for data with no response
#'
#'@param formula an object of type formula
#'@param df an object of type dataframe
#'@return a model matrix with data from df
construct_responseless_model_matrix <- function(formula, df){
  responseless <- delete.response(terms(formula))
  return(model.matrix(responseless, model.frame(responseless, data = df)))
}

#' get posterior of blm object
#'
#' This function fetches the mean(s) and covariance matrix posterior distribution
#'
#' @param blm an object of class blm
#' @return a list of mean(s) and covariance matrix
posterior <-function(blm){
  post <- list()
  post$mean <- blm$mean
  post$variance <- blm$variance
  post
}

#'bayesian linear regression function
#'
#'main constructor function for creating an object of class blm
#'
#'@param formula an object of type formula
#'@param df an object of type dataframe
#'@param prior an object of type blm
#'@param beta the precision
#'@return an object of class blm
#'
#'@export
blm <- function(formula, df, beta = 1.3, prior = NULL){
  call <- match.call()
  if(is.null(prior)) {alpha = 0.5}
  else {alpha <- posterior(prior)$var}
  model_matrix <- construct_model_matrix(formula, df = df)
  Sxy <- diag(alpha, nrow = ncol(model_matrix)) + beta * t(model_matrix) %*% model_matrix
  Sxy <- solve(Sxy)
  Mxy <- beta * Sxy %*% t(model_matrix) %*% model.response(model.frame(formula, df))
  object <- list()
  object$call <- call
  object$variance <- Sxy
  object$mean <- Mxy
  class(object) <- 'blm'
  object
}

#' predict new values from blm
#'
#'A function for predicting values of y from new data
#'
#'@param blm an object of class blm
#'@param x an object of type dataframe
#'@return predicted y values from x
#'
#'@export
predict.blm <- function(blm, x){
  model_matrix <- construct_responseless_model_matrix(formula(formula(blm$call)), df = x)
  y <- t(blm$mean) %*% t(model_matrix)
  y
}
