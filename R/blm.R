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
#' @param blm an object of class blm
#' @return a list of mean(s) and covariance matrix
posterior <-function(blm) blm$posterior

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
blm <- function(formula, beta = 1.3, prior = NULL, Data = NULL){
  call <- match.call()
  if(is.null(prior)) {alpha = 0.5}
  else {alpha <- posterior(prior)$var}
  if(is.null(Data)){ model_matrix <-construct_model_matrix(formula)
                     response <- model.response(model.frame(formula))}
  else{
                     model_matrix <- construct_model_matrix(formula, Data = Data)
                     response <- model.response(model.frame(formula, Data))}
  Sxy <- diag(alpha, nrow = ncol(model_matrix)) + beta * t(model_matrix) %*% model_matrix
  Sxy <- solve(Sxy)
  Mxy <- beta * Sxy %*% t(model_matrix) %*% response
  posterior <- list()
  object <- list()
  object$call <- call
  posterior$variance <- Sxy
  posterior$mean <- Mxy
  object$posterior <- posterior
  class(object) <- 'blm'
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
  if (!inherits(object, "blm"))
    warning("calling predict.blm(<fake-blm-object>) ...")
  if (is.null(newdata)){ Data <- model.matrix(object)}
  else{ Data <- newdata}
  model_matrix <- construct_responseless_model_matrix(formula(formula(object$call)), Data = Data)
  posterior <- posterior(object)
  y <- t(posterior$mean) %*% t(model_matrix)
  y
}



