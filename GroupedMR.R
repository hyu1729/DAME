GroupedMR <- function(data, unmatched_indices, covs, repeats = FALSE) {

  # From groups on data by exact matching on covs
  # bit_match copied from FLAME package

  if (repeats) {
    match_out <- bit_match(data, covs)
  }
  else {
    match_out <- bit_match(data[!data$matched,], covs)
  }
  match_index <- match_out[[1]]
  index <- match_out[[2]]

  # Get group values
  unique_matched_vals = unique(index)

  # Check whether unmatched data matches with group values
  # First compute b_u for data_unmatched
  n_levels <- sapply(data[unmatched_indices, covs, drop = FALSE], nlevels)
  data_unmatched_wo_t <- gmp::as.bigz(as.matrix(data[unmatched_indices, covs[order(n_levels)]]))
  n_levels <- sort(n_levels)

  multiplier <- gmp::pow.bigz(n_levels, seq_along(n_levels) - 1)

  b_u <-
    as.vector(gmp::`%*%`(data_unmatched_wo_t, multiplier))

  # Check if b_u matches any group values
  new_matches <- as.vector(as.character(b_u)) %in% as.vector(as.character(unique_matched_vals))

  new_matched_indices <- unmatched_indices[new_matches]

  # Get D_m, same format as output for bit_match
  # data_matched <- list(data_unmatched_matches, b_u[data_unmatched_matches])

  return(list(new_matched_indices = new_matched_indices,
               match_out= match_out))
}


# Helper Functions:

aggregate_table <- function(vals) {
  vals <- as.character(vals)
  tab <- table(vals)
  name <- names(tab)
  return(as.vector(tab[match(vals, name)]))
}

# bit_match takes a dataframe, a set of covariates to match on, the
# treatment indicator column and the matched indicator column. it returns the
# array indicating whether each unit is matched (the first return value), and a
# list of indices for the matched units (the second return value). dataframe must
# have a matched indicator column

bit_match <- function(data, covs) {
  # gmp::as.bigz is for handling lots and lots of covariates so we don't
  # have trouble with overflow issues

  n_levels <- sapply(data[, covs, drop = FALSE], nlevels)
  data_wo_t <- gmp::as.bigz(as.matrix(data[, covs[order(n_levels)]]))
  n_levels <- sort(n_levels)

  # Compute b_u
  multiplier <- gmp::pow.bigz(n_levels, seq_along(n_levels) - 1)

  b_u <-
    as.vector(gmp::`%*%`(data_wo_t, multiplier))

  # Compute b_u+
  b_u_plus <-
    as.vector(gmp::add.bigz(data$treated, gmp::mul.bigz(b_u, as.bigz(n_levels))))

  # Compute c_u
  c_u <- aggregate_table(b_u)

  # Compute c_u+
  c_u_plus <- aggregate_table(b_u_plus)

  match_index <- (c_u != c_u_plus) & (c_u >= 2)

  index <- b_u[match_index]

  return(list(match_index = match_index,
              index = index))
}
