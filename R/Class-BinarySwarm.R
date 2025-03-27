#' @title BinarySwarm
#'
#' @description
#' Abstract base class for a Binary Swarm
#'
#' @import R6
#'
R6::R6Class("BinarySwarm",
            public = list(
              #' @field population_size integer, size of population
              population_size = NULL,
              #' @field n_bits integer, number of bits to represent
              n_bits = NULL,
              #' @field transferFun object of class `TransferFunction` to update from continuous to binary
              transferFun = NULL,
              #' @field global_best list with 2 elements with names `Position` and `OptimResult`
              global_best = NULL,
              #' @field population list with `Particle` objects
              population = NULL,
              #' @field has_run logical, if the process has been run
              has_run = NULL,
              #' @field iteration integer, position of the iteration during the process
              iteration = NULL,
              #' @field max_iteration integer, maximum number of iterations to perform
              max_iteration = NULL,
              #' @field use_var_importance logical, if the variable importance should be used
              use_var_importance = NULL,
              #' @field variable_importance (Experimental) list of length `n_bits` to store the feature importance
              variable_importance = NULL,
              #' @description
              #' Creation of object of a Binary Swarm-based algorithm
              #'
              #' @param population_size integer, size of the population
              #' @param n_bits integer, number of bits to be represented
              #' @param transferFun object of class `TransferFunction`
              #' @param particleGenerator object of class `BPG`
              #' @param seed seed for the generation of the population
              #' @param use_var_importance (Experimental) logical
              #' If variable importance can be done using transfer functions
              initialize = function(population_size,
                                   n_bits,
                                   transferFun,
                                   particleGenerator,
                                   use_var_importance = FALSE,
                                   seed = NULL){
                # checks
                if(!("BPG" %in% class(particleGenerator))){
                  stop("particleGenerator should be of class 'BPG' ")
                }

                if(!("TransferFunction" %in% class(transferFun))){
                  stop("transferFun should be of class 'TransferFunction'")
                }

                if(n_bits != as.integer(n_bits)){
                  stop("n_bits should be an integer")
                }
                if(population_size != as.integer(population_size)){
                  stop("population_size should be an integer")
                }
                # end checks

                self$population_size <- population_size
                self$n_bits <- n_bits
                self$transferFun <- transferFun
                self$use_var_importance <- use_var_importance

                self$variable_importance <- list()

                # set the global best
                self$global_best <- list("Position" = 0,
                                            "Result" = -Inf)
                self$has_run <- FALSE
                self$iteration <- 1
                self$max_iteration <- 30

                # particleGenerator returns list of elements
                self$population <- particleGenerator$get(population_size,
                                                            n_bits,
                                                            seed = seed)
              },
              #' @description
              #' Get the current population of the swarm
              #' @return list with Particles with the population
              get_population = function(){
                return(self$population)
              },
              #' @description
              #' Get the global best after a run
              get_global_best = function(){
                if (!self$has_run){warning("Process hasn't run yet")}
                return(self$global_best)
              },
              #' @description
              #' Returns the iteration during the process
              #' Needed if the algorithm needs the iteration number
              get_iteration = function(){return(self$iteration)},
              #' @description
              #' Returns logical if it uses to calculate the variable importance
              get_use_variable_importance = function(){
                return(self$use_var_importance)
              },
              #' @description
              #' Returns the variable importance
              #' This is just testing the results
              get_variable_importance = function(){
                return(self$variable_importance)
              },
              #' @description
              #' (Experimental)
              #' Add values to the variable importance
              #'
              #' @param x list with the values for the variable importance
              add_variable_importance = function(x){
                # add the results to the variable importance
                position <- length(self$variable_importance) + 1
                self$variable_importance[[position]] <- matrix(unlist(x),
                                                                  ncol = self$n_bits)
              },
              #' @description
              #' Run the Swarm-based algorithm
              #'
              #' @param fun function to maximize
              #' @param args_fun extra arguments for the function
              #' @param max_iter integer,
              #' maximum number of iterations to take place
              #' @param max_stable integer,
              #' maximum number of iterations the results needs to be stable to early stop the result
              #' @param seed a single value, interpreted as an integer, or NULL
              #' @param is_maximize logical, if the solution is a maximization (TRUE);
              #' if FALSE it is a minimization
              #' @param show_process logical, if some in between output should be shown
              #' @param pos_fun Output of `fun` can give multiple results as a list.
              #' This should give the position of the list which it should maximize/minimize
              #'
              #' @return
              #' list with 2 values `AllResults` and `BestResult`
              #' `BestResult` is list of 2 elements with the position and the result of the position
              run_process = function(fun,
                                     args_fun = list(),
                                     max_iter = 30,
                                     max_stable = 5,
                                     pos_fun = 1,
                                     is_maximize = TRUE,
                                     show_process = TRUE,
                                     seed = NULL){

                set.seed(seed)

                self$has_run <- TRUE
                self$max_iteration <- max_iter

                all_results <- list()
                counter_stable <- 0

                for(i in 1:max_iter){
                  if (show_process) {
                    cat(sprintf("\r Iteration: %d", i))
                  }

                  self$iteration <- i

                  # get all the results of the particles
                  results <- self$get_results(fun,
                                                 args_fun,
                                                 pos_fun = pos_fun,
                                                 is_maximize = is_maximize)

                  has_global_best <- self$set_global_best(results)

                  # update All results
                  all_results <- append(all_results,
                                       list(results))

                  # keep a counter if the numbers are stable
                  counter_stable <- ifelse(has_global_best, 1, counter_stable + 1)

                  # it is stable enough, so get out the first loop
                  if(counter_stable >= max_stable){break}

                  # update the positions of all particles
                  self$update_all_positions()

                }
                if (show_process) {
                  # put the other everything on a new line
                  cat("\n")
                }

                # convert the results into a proper format (data.table)
                total_list <- data.table()
                for (iter in 1:length(all_results)) {

                  result_iter <- rbindlist(all_results[[iter]]$AllResults)
                  result_iter[, Iteration := iter]

                  total_list <- rbindlist(list(total_list,
                                               result_iter))
                }

                return(list(AllResults = total_list,
                            BestResult = self$global_best))
              },
              #' @description
              #' Interface function.
              #' This specific function of this object should not be used
              update_all_positions = function(){
                stop("Should be implemented in the subclass")
              },
              #' @description
              #' set the global best based on the results given
              #' @param results list of results.
              #' List of the result should have 3 elements
              #' `Positions`, `Results` & `AllResults`
              #' based on output of function `get_results()`
              #' @return TRUE if the global best is updated
              #' FALSE when the global best did not change
              set_global_best = function(results){
                changed_global <- FALSE
                # update Global best
                for(i in seq_along(results[["Result"]])){
                  result <- results[["Result"]][[i]]
                  position <- results[["Positions"]][[i]]

                  if (result > self$global_best[["Result"]]) {
                    self$global_best <- list("Position" = position,
                                                "Result" = result)
                    changed_global <- TRUE
                  }
                }
                return(changed_global)
              },
              #' @description
              #' get the results of the function to maximize
              #'
              #' @param fun function to maximize,
              #' with a binary sequence as first input.
              #' Function to maximize can give a list of results.
              #' Only the first element will be used to maximize.
              #' This can be useful in case other results should to be monitored but should not optimized.
              #' @param argsFun list with arguments passed to `fun`
              #' @param pos_fun position to be used of the output
              #' @param is_maximize logical, if the result should be maximized.
              #' if FALSE, the negative will be given.
              #'
              #' @return list with 3 elements with names:
              #' `Result` list with all the output results
              #' `Positions` list with all the positions checked
              #' `AllResults` list with the results which should be optimized
              #' The results and positions are in the same order of each other
              get_results = function(fun,
                                     argsFun,
                                     pos_fun = 1,
                                     is_maximize = TRUE){
                results <- list() # results to optimize
                positions <- list() # positions
                all_result <- list() # all results that are given

                for(part in self$population){
                  # first get the position
                  position <- part$get_position()
                  results_part <- do.call(fun,
                                          args = append(list(position), argsFun))
                  # first element is to maximize
                  result <- results_part[[pos_fun]]

                  # then update result so that it would maximize
                  if (!is_maximize) {
                    result <- -result
                  }

                  #save result in the particle
                  part$save_result(result)

                  # store the results
                  positions <- append(positions, list(position))
                  results <- append(results, list(result))
                  all_result <- append(all_result, list(results_part))
                }
                return(list("Result" = results,
                            "Positions" = positions,
                            "AllResults" = all_result))
              })) -> BinarySwarm
