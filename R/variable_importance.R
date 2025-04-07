#' Get the variable importance of a variable based on the Meta Heuristic run.
#'
#' @param results data.table with the results of the meta heuristic run.
#' The results should be in line with the transformation of `result_2_dt`
#' @param column_name column name of `results` with the results of the Meta heuristic run
#' @param col_binary column name of `results` where the binary coding is of the results
#' @param var_names vector of labels for the binary coding of the positions.
#' Length of `var_names` should be the same as the binary coding
#'
#' @import data.table
#' @import matrixStats
#' @importFrom stats rnorm
#' @importFrom purrr map2
#' @returns list, with names if `var_names` is given, with order as the importance of the variable.
#' @export
variable_importance <- function(results,
                                column_name,
                                col_binary = "Position",
                                var_names = NULL){

  c_results <- copy(results)

  c_results[, ConvertedPosition := lapply(strsplit(get(col_binary) ,
                                                   split = ""),
                                          as.numeric)]
  c_results[, ValuePosition := map2(ConvertedPosition,
                                    get(column_name),
                                    \(x, y) x * y)]

  n_features <- length(c_results$ValuePosition[[1]])

  # construction of a matrix
  # unlist is necessary due to the construction of the results
  mat_test <- matrix(unlist(c_results$ValuePosition),
                     ncol = n_features,
                     byrow = TRUE)

  # get the average result when used
  avg_when_used <- colSums(mat_test) / colSums(mat_test != 0)
  # case that the feature is never used
  avg_when_used[is.nan(avg_when_used)] <- 0

  # get the average result when variable is not used
  avg_not_used <- c()

  for (iter in 1:n_features) {

    # average when the column is not used
    avg_not_used_iter <- mean(rowMaxs(matrix(mat_test[which(mat_test[,iter] == 0), ],
                                             ncol = n_features)))
    avg_not_used <- c(avg_not_used,
                      avg_not_used_iter)
  }

  result <- avg_when_used - avg_not_used

  if (!is.null(var_names)) {
    names(result) <- var_names
  }

  return(result)
}

