GenerateNewActiveSets <- function(s, delta) {
  # 1
  k = length(s)
  # 2
  Z <- list()

  # 3
  if (length(delta) == 0) return (Z)
  length_list <- lapply(delta, length)
  subset_index <- which(length_list == k)
  delta_k <- list()
  for (i in subset_index) {
    delta_k[[length(delta_k) + 1]] <- unlist(delta[[i]])
  }
  delta_k[[length(delta_k) + 1]] <- unlist(s)


  # 5
  supp <- table(unlist(unique(delta_k)))
  # 6
  omega <- setdiff(strtoi(names(supp[supp >= k])), s)

  # 7
  if (all(supp, supp[names(supp[s])] >= k)) {
    # 8
    for (a in omega) {
      # 9
      r = sort(c(s, a))
      # 10
      if (all(combn(r, k, simplify = FALSE) %in% delta_k)) {
        # 11
        Z[[length(Z) + 1]] <- r
      }
    }
  }

  # 12
  return(Z)
}
