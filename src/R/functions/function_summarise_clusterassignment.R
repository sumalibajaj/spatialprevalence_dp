# Cluster summary
summary_cluster_assign <- function(results, maxIters, data){
  # Keep only columns of cluster assignments
  assignments <- results
  
  # Fake assignment matrix to see if code is doing the right thing
  # assignments <- matrix(c(10, 10, 21,
  #                         23, 5, 23), nrow = 3, ncol = 2)
  # maxIters <- 2
  
  # Create an empty matrix
  matrix_size <- nrow(assignments)
  matrix <- matrix(0, nrow = matrix_size, ncol = matrix_size)
  
  # 
  # Iterate over assignments in the second half
  # for (col_i in ((maxIters %/% 2) + 1):maxIters) {
  for (col_i in (maxIters-9):maxIters) {
    assignment <- assignments[, col_i]
    
    # Iterate over characters and indices in the assignment
    for (ii in seq_along(assignment)) {
      for (jj in seq_along(assignment)) {
        # Check if indices are equal
        if (ii == jj) {
          matrix[ii, jj] <- matrix[ii, jj] + 1
        } else if (ii != jj && assignment[ii] == assignment[jj]) {
          matrix[ii, jj] <- matrix[ii, jj] + 1
        }
      }
    }
  }
  
  # Similarity matrix between 0 and 1
  # matrix <- matrix/(maxIters-(maxIters %/% 2))
  matrix <- matrix/10
  
  # number of clusters at each iteration
  len_unique <- function(x) length(unique(x))
  num_clusters <- apply(results, MARGIN = 2, FUN = len_unique)
  mode_num_clusters <- names(sort(-table(num_clusters)))[1] %>% as.numeric() # MODE!
  mode_num_clusters
  
  matrix_kernel <- as.kernelMatrix(matrix)
  # Perform spectral clustering on this similarity or affinity matrix
  final_clusters <- specc(matrix_kernel, centers = mode_num_clusters, iterations = 10000)
  final_clusters <- as.numeric(final_clusters)
  
  length(unique(final_clusters)) == mode_num_clusters # check final no. of clusters is as expected

  # # assign final cluster assignments to prevalence data
  data_cluster <- cbind(data, final_clusters)
  
  return(data_cluster)
}
