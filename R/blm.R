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

#'bayesian linear regression function
#'
#'main constructor function for creating an object of class blm
#'
#'@param formula an object of type formula
#'@param df an object of type dataframe
#'@param alpha a number
#'@param beta a number
#'@return an object of class blm
#'
#'@export
blm <- function(formula, df, alpha = 0.5, beta = 1.3){
  call <- as.list(sys.call())
  model_matrix <- construct_model_matrix(formula, df = df)
  Sxyneg <- diag(alpha, nrow = ncol(model_matrix)) + beta * t(model_matrix) %*% model_matrix
  Sxy <- Sxyneg * Sxyneg^2
  Mxy <- beta * Sxy %*% t(model_matrix) %*% model.response(model.frame(formula, df))
  posterior <- list()
  posterior$call <- call[[2]]
  posterior$Sxy <- Sxy
  posterior$Mxy <- Mxy
  class(posterior) <- 'blm'
  posterior
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
  model_matrix <- construct_responseless_model_matrix(formula(blm$call), df = x)
  y <- t(blm$Mxy) %*% t(model_matrix)
  y
}
