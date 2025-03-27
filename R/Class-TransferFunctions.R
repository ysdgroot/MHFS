#' Construction of the S-shaped and V-shaped transfer functions
#'
#' @return R6 Class object of TransferFunction
#' @export
R6::R6Class("TransferFunction",
            public = list(
              #' @field name character. Name of the transfer function
              name = NULL,
              #' @field fun function. function to transform into values between 0 and 1
              fun = NULL,
              #' @field type character, "S" or "V", if it is S-shaped or V-shaped transfer class
              type = NULL,
              #' @param name character. Name of the transfer function
              #' @param fun function. function to transform into values between 0 and 1
              #' @param type character, "S" or "V", if it is S-shaped or V-shaped transfer class
              #'
              #' @return R6 object
              initialize = function(name,
                                    fun,
                                    type = "S"){
                self$name <- name
                self$fun <- fun

                if(!(type %in% c("S", "V"))){
                  stop(sprintf("Type should be 'S' or 'V' not %s", type))
                }
                self$type <- type
              },
              #' Transfer the values
              #'
              #' @param x vector or matrix to calculate the function on
              #'
              #' @return results of function `fun`
              transfer = function(x){
                  return(self$fun(x))
                },

              #' Update the next positions based on the transfer functions
              #'
              #' @param x binary vector of current position
              #' @param velocity vector of the next position which is also called the velocity.
              #' based on the velocity the next binary positions are produced
              #' @param use_var_importance (Experimental) logical
              #' If variable importance can be done using transfer functions
              #' @return transfer the final results
              changePosition = function(x, velocity, use_var_importance = FALSE){
                rand <- runif(length(x))
                if("matrix" %in% class(x)){
                  rand <- matrix(rand,
                                 nrow = dim(x)[1],
                                 ncol = dim(x)[2])
                }

                transferedy <- self$transfer(velocity)
                if(self$type == "S"){
                  if (use_var_importance){
                    return(list(1 * (rand < transferedy),
                                "TransferValues" = transferedy))
                  }
                  return( 1 * (rand < transferedy))
                } else if(self$type == "V"){
                  if (use_var_importance){
                    return(list(1 * ((rand < transferedy) * abs(x - 1) +
                                       (rand >= transferedy) * x),
                                "TransferValues" = 1 * ((rand < transferedy) * transferedy  +
                                                          (rand >= transferedy) * (1-transferedy)) *
                                  ( 1 - 2 * x)
                                ))
                  }
                  return(1 * ((rand < transferedy) * abs(x - 1) +
                           (rand >= transferedy) * x))
                } else{
                  stop("Type is not correct")
                }
              },
#' @description
#' Get the name of the Transfer Function
#' @return character with the name
              get_name = function(){
                return(self$name)
              },

#' @description
#' Get the Type of transfer function (S-shaped or V-shaped)
#' @return character with the name
              get_type = function(){
                return(self$type)
              },
#' @description
#' Get the Transfer Function as function
#' @return Function
               get_fun = function(){
                return(self$fun)
              }
            )) -> TransferFunction
