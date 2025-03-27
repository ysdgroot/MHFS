#' @title Particle
#'
#' @description
#' Abstract class for a particle
R6::R6Class("Particle",
            public = list(
              #' @field position n-dimensional value which states the position
              position = NULL,
              #' @field track History of positions of this particle
              track = NULL,
              #' @field personal_best Best position or the past
              personal_best = NULL,
              #' @param position vector with only 0 or 1
              initialize = function(position){

                if(length(setdiff(unique(position), c(0, 1))) != 0){
                  stop("Position should only have 0 or 1 values")
                }

                self$position <- position
                self$track <- list()
              },
              #' @description
              #' save the result of the current position
              #'
              #' @param result a single numeric value
              #'
              #' @return returns nothing
              save_result = function(result){

                add_track <- list("Position" = self$position,
                                  "Result" = result)

                self$track <- append(self$track,
                                        list(add_track))

                if(is.null(self$personal_best) ||
                   self$personal_best[["Result"]] < result){
                  self$personal_best <- add_track
                }

              },
              #' @description
              #' return the best result of the particle
              #' @return list with 2 elements with names
              #' `Position` and `Result`
              get_personal_best = function(){
                return(self$personal_best)
              },
              #' @description
              #' return the current position of the particle
              #' @return position, vector of a length
              get_position = function(){
                return(self$position)
              },
              #' @description
              #' return the current position of the particle
              #' @param position new position of the particle
              #' @return position, vector of a length
              set_position = function(position){
                if(length(position) != length(self$position)){
                  stop("Different length compared to current position")
                }
                if(!self$check_position(position)){
                  stop("New position should only have 0 or 1")
                }
                self$position <- position
              },
              #' @description
              #' returns list of the history of this particle
              #'
              #' @return list of lists with 2 elements with names
              #' `Position` and `Result`
              get_track = function(){
                return(self$track)
              },

              #' @description
              #' Get string with the current position#'
              #'
              #' @return string of length `length(position)`
              get_char_position = function(){
                return(paste(self$position, collapse = ""))
              },
              #' @description
              #' check if the new position only has 0 and 1
              #' @param position vector to be checked
              #' @return `logical(1)` TRUE if the position only contains 1 and/or 0
              check_position = function(position){
                return(length(setdiff(unique(position), c(0, 1))) == 0)
              }
            )) -> Particle

#' @title ParticleVelocity
#'
#' @description
#' Abstract class for a particle with velocity.
#' This inherits the functionality of a simple particle
R6::R6Class("ParticleVelocity",
            inherit = Particle,
            public = list(
              #' @field velocity n-dimensional velocity of the particle
              velocity = NULL,
              #' @description
              #' Initialize the Particle with Velocity
              #'
              #' @param position n-dimensional position of the particle
              #' @param velocity n-dimensional velocity of the particle
              #'
              #' @return NULL
              initialize = function(position,
                                    velocity){
                super$initialize(position)

                if(length(position) != length(velocity)){
                  stop("Position and Velocity should have the same length")
                }

                self$velocity <- velocity
              },

              #' @description
              #' Get the current Velocity
              #'
              #' @return vector of fixed length
              get_velocity = function(){
                return(self$velocity)
              },
              #' @description
              #' Set the current velocity into a new one
              #' @param velocity new velocity of the particle
              #' @return NULL
              set_velocity = function(velocity){
                if(length(velocity) != length(self$velocity)){
                  stop("Different length compared to current velocity")
                }
                self$velocity <- velocity
              },
              #' @description
              #' Set the position and velocity into a new one
              #' @param position new position of the particle
              #' @param velocity new velocity of the particle
              #' @return NULL
              set_position_velocity = function(position,
                                               velocity){
                super$set_position(position)
                self$set_velocity(velocity)
              }
            )) -> ParticleVelocity

