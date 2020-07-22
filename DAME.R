### NOTE: Edit drop to include adaptive weights from machine learning
decide_drop <-
  function(data, active_cov_sets, weights) {
    # Helper function for DAME, decides what cov set to drop

    max_weight <- 0
    s_h <- list()
    p <- length(weights)
    J <- c(1:p)
    for (s in active_cov_sets) {
      theta_s <- as.integer(!(J %in% s))
      temp <- theta_s %*% weights
      if (temp > max_weight) {
        max_weight <- temp
        s_h <- s
      }
    }

    return(s_h)
  }

make_MGs <-
  function(data, index, matched_units, covs, cov_names) {

    # Takes all the units that were matched on these p covariates and separates
    # them into matched groups based off their unique values of those covariates
    # Returns a list with two items
    # MGs: a vector, each entry of which is a list corresponding to a different MG
    #   and contains the indices of the corresponding units
    # matched_on: a vector the same length as MGs, each entry of which is a list
    #  detailing the covariates and their values that the units in the
    #  corresponding MG matched on

    unique_MGs <- unique(index)
    n_MGs <- length(unique_MGs)

    MGs <- vector('list', length = n_MGs)
    matched_on <- vector('list', length = n_MGs)

    for (i in 1:n_MGs) {
      members <- matched_units[which(index == unique_MGs[i])]
      MGs[[i]] <- members
      matched_on[[i]] <- data[members[1], covs, drop = FALSE]
      rownames(matched_on[[i]]) <- NULL
      names(matched_on[[i]]) <- cov_names[covs]
    }

    return(list(MGs = MGs,
                matched_on= matched_on))
  }

DAME <- function(data, w, repeats = FALSE, num_iter) {
  h <- 1
  data_all <- data
  unmatched_indices <- c(1:length(data[[1]]))
  MGs <- vector('list', length = num_iter)
  matched_on <- vector('list', length = num_iter)
  data_matched <- vector('list', length = num_iter)
  p <- length(w)
  active_cov_sets <- list()

  for (i in 1:p) {
    active_cov_sets <- append(active_cov_sets, c(i))
  }
  processed_cov_sets <- list()
  J = c(1:p)

  # First try grouping on all covs
  match <- GroupedMR(data_all, unmatched_indices, J, repeats)
  new_matched_indices <- match[[1]]
  match_index <- match[[2]][[1]]
  index <- match[[2]][[2]]

  data_all$matched[new_matched_indices] <- TRUE

  made_matches <- sum(match_index) > 0

  if (made_matches) {
    new_MGs <- make_MGs(data_all, index, new_matched_indices, J, colnames(data[,J]))
    MGs[[h]] <- new_MGs$MGs
    matched_on[[h]] <- new_MGs$matched_on
    data_matched[[h]] <- new_matched_indices
  }

  unmatched_indices <- setdiff(unmatched_indices, new_matched_indices)
  h <- h + 1

  # Main loop
  while (length(unmatched_indices) > 0 & h <= num_iter) {

    print(data_all)

    # Find best covariate-set to drop from active covariate sets
    curr_cov_set <- decide_drop(data, active_cov_sets, w)

    # Check early stopping condition theta_s * w < 5% sum(w)
    ### EDIT
    if (sum(w) * 0.05 > as.integer(!(J %in% curr_cov_set)) %*% w) {
        break
      }
    cov <- setdiff(J, curr_cov_set)


    match <- GroupedMR(data_all, unmatched_indices, cov, repeats)
    new_matched_indices <- match[[1]]
    match_index <- match[[2]][[1]]
    index <- match[[2]][[2]]

    data_all$matched[new_matched_indices] <- TRUE

    made_matches <- sum(match_index) > 0

    if (made_matches) {
      new_MGs <- make_MGs(data_all, index, new_matched_indices, cov, colnames(data[,cov]))
      MGs[[h]] <- c(MGs[[h - 1]], new_MGs$MGs)
      matched_on[[h]] <- c(matched_on[[h - 1]],  new_MGs$matched_on)
      data_matched[[h]] <- c(data_matched[[h - 1]], new_matched_indices)
    }

    Z_h <- GenerateNewActiveSets(curr_cov_set, processed_cov_sets)
    active_cov_sets <- setdiff(active_cov_sets, curr_cov_set)
    active_cov_sets <- append(active_cov_sets, Z_h)
    processed_cov_sets <- append(processed_cov_sets, curr_cov_set)
    unmatched_indices <- setdiff(unmatched_indices, new_matched_indices)
    h <- h + 1
  }

  return(list(data_matched = data_matched,
              MGs = MGs,
              matched_on = matched_on))
}
