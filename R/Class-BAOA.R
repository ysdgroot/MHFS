# Particle ----------------------------------------------------------------

#' R6 Object Representing a Particle for the BAOA algorithm
#' @export
R6::R6Class("ParticleBAOA",
            inherit = Particle) -> ParticleBAOA


# Swarm for the particle --------------------------------------------------

#' R6 Object Representing the Swarm for the BAOA algorithm
#'
#' @description
#' The Swarm based on the Binary Arithmetic Optimizer Algorithm with
#' @export
R6::R6Class("SwarmBAOA",
            inherit = BinarySwarm,
            public = list(
              #' @field minMoa minimum value for MOA
              minMoa = NULL,
              #' @field maxMoa maximum value for MOA
              maxMoa = NULL,
              #' @field beta beta for the BAOA
              beta = NULL,
              #' @field k beta for the BAOA
              k = NULL,
              #' @field delta beta for the BAOA
              delta = NULL,
              #' @description
              #' initialize the Binary Arithmetic Optimizer Algorithm (BAOA)
              #'
              #' @param population_size integer, size of the population
              #' @param n_bits integer, number of bits to be represented
              #' @param transferFun object of class `TransferFunction`
              #' @param particleGenerator object of class `BPG`
              #' @param seed seed for the generation of the population
              #' @param beta Parameter for the BAOA Algoritm
              #' @param k Parameter for the BAOA Algoritm
              #' @param minMoa Parameter for the BAOA Algoritm
              #' @param maxMoa Parameter for the BAOA Algoritm
              #' @param delta Parameter for the BAOA Algoritm
              #'
              #' @returns NULL
              initialize = function(population_size,
                                    n_bits,
                                    transferFun,
                                    particleGenerator,
                                    beta,
                                    k,
                                    minMoa,
                                    maxMoa,
                                    delta = 1e-8,
                                    seed = seed){

                if (!("BPG" %in% class(particleGenerator))){
                  stop("Particle should be of class 'BPG'")
                }

                super$initialize(population_size,
                                 n_bits,
                                 transferFun,
                                 particleGenerator,
                                 seed = seed)

                self$beta <- beta
                self$k <- k
                self$delta <- delta
                self$minMoa <- minMoa
                self$maxMoa <- maxMoa
              },
              #' @description
              #' Math optimizer accelerated (MOA)
              #' Function which is used the determine the probability of the exploration or exploitation phase
              #' @param t integer, current step of iterations
              #' @param Tmax the maximum number of iterations
              #' @returns numeric value
              MOA = function(t, Tmax){
                return(self$minMoa + t * ((self$maxMoa - self$minMoa)/Tmax))
              },
              #' @description
              #' Math optimizer probability (MOP)
              #' Additional value to determine the search
              #' @param t integer, current step of iterations
              #' @param Tmax the maximum number of iterations
              #' @returns numeric value
              MOP = function(t, Tmax){
                return(1 - t ** (1/self$beta)/(Tmax**(1/self$beta)))
              },
              #' @description
              #'  Calculate the new positions
              #' @param t integer, current step of iterations
              #' @param Tmax the maximum number of iterations
              #' @returns matrix with all the new positions.
              #' rows represents a particle.
              #' columns represents the positions in that dimension
              new_positions = function(t,
                                       Tmax){

                # get MOA(t)
                moa_t <- self$MOA(t, Tmax)

                # set the dimensions of the matrix
                n <- self$popSize
                m <- self$nBits

                # exploration phase when r > MOA(t)
                rand1 <- matrix(runif(n * m),
                                nrow = n,
                                ncol = m)

                # vector of values
                isExploration <- rand1 > moa_t

                # generate random values for the type of calculation
                rand2 <- matrix(runif(n * m),
                                nrow = n,
                                ncol = m)

                # the MOP of at iteration t
                mop_t <- self$MOP(t, Tmax)

                # the original division is
                # (MOP + delta) * ((ub_j - lb_j) * k + lb_j)
                # but because of binary levels ub_j = 1 and lb_j = 0
                division <- (mop_t + self$delta) * self$k

                # original
                # MOP * ((ub_j - lb_j) * k + lb_j)
                multi_add_minus <- mop_t * self$k

                # calculated everything at once
                # matrix addition and multiplication is iterating row and then column
                # not in the same 'shape' as the vector
                matrix_global_best <- t(matrix(rep(self$globalBest$Position[[1]],
                                                   self$popSize),
                                               ncol = self$popSize))

                result <- isExploration *
                  matrix_global_best *
                  ((rand2 < 0.5) * 1/division +
                     (rand2 >= 0.5) * multi_add_minus) +

                  !isExploration *
                  ( matrix_global_best +
                      multi_add_minus *
                      ((rand2 < 0.5) * (-1) +
                         (rand2 >= 0.5)))

                return(result)
              },

              #' @description
              #' Update all the positions based on the parameters in the object
              #' @returns NULL.
              #' All positions needs to be changed inside the R6 object
              update_all_positions = function(){

                t <- self$iteration
                Tmax <- self$max_iteration

                # get MOA(t)
                moa_t <- self$MOA(t, Tmax)

                # exploration phase when r > MOA(t)
                rand1 <- matrix(runif(self$population_size * self$n_bits),
                                nrow = self$population_size,
                                ncol = self$n_bits)

                # vector of values
                isExploration <- rand1 > moa_t

                # generate random values for the type of calculation
                rand2 <- matrix(runif(self$population_size * self$n_bits),
                                nrow = self$population_size,
                                ncol = self$n_bits)

                # the MOP of at iteration t
                mop_t <- self$MOP(t, Tmax)

                # the original division is
                # (MOP + delta) * ((ub_j - lb_j) * k + lb_j)
                # but because of binary levels ub_j = 1 and lb_j = 0
                division <- (mop_t + self$delta) * self$k

                # original
                # MOP * ((ub_j - lb_j) * k + lb_j)
                multi_add_minus <- mop_t * self$k

                # calculated everything at once
                # matrix addition and multiplication is iterating row and then column
                # not in the same 'shape' as the vector
                matrix_global_best <- t(matrix(rep(self$global_best[["Position"]],
                                                   self$population_size),
                                               ncol = self$population_size))

                result <- isExploration *
                  matrix_global_best *
                  ((rand2 < 0.5) * 1/division +
                     (rand2 >= 0.5) * multi_add_minus) +

                  !isExploration *
                  ( matrix_global_best +
                      multi_add_minus *
                      ((rand2 < 0.5) * (-1) +
                         (rand2 >= 0.5)))

                i <- 1
                # update position of each particle
                for(particle in self$population){
                  result_i <- result[i, ]
                  position <- particle$get_position()
                  new_position <- self$transferFun$changePosition(position,
                                                                     result_i)

                  particle$set_position(new_position)
                  i <- i + 1
                }
              }
            )) -> SwarmBAOA
