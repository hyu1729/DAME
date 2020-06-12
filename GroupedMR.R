GroupedMR <- function(data, data_unmatched, covs) {
  # From groups on data by exact matching on covs
  # bit_match copied from FLAME package
  match_out <- bit_match(data, covs)
  match_index <- match_out[[1]]
  index <- match_out[[2]]

  # Get group values
  unique_matched_vals = unique(index)

  # Check whether unmatched data matches with group values
  # First compute b_u for data_unmatched
  n_levels <- sapply(data_unmatched[, covs, drop = FALSE], nlevels) - 1
  data_unmatched_wo_t <- gmp::as.bigz(as.matrix(data_unmatched[, covs[order(n_levels)]]))
  n_levels <- sort(n_levels)

  multiplier <- gmp::pow.bigz(n_levels, seq_along(n_levels) - 1)

  b_u <-
    gmp::`%*%`(data_unmatched_wo_t, multiplier) %>%
    as.vector()

  print(b_u)
  print(unique_matched_vals)
  # Check if b_u matches any group values
  data_unmatched_matches <- as.vector(as.character(b_u)) %in% as.vector(as.character(unique_matched_vals))

  # Get D_m, same format as output for bit_match
  data_matched <- list(data_unmatched_matches, b_u[data_unmatched_matches])

  return(list(data_matched, match_out))
}


# Dependencies:

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
  # subtract one for the matched indicator column
  n_levels <- sapply(data[, covs, drop = FALSE], nlevels) - 1
  data_wo_t <- gmp::as.bigz(as.matrix(data[, covs[order(n_levels)]]))
  n_levels <- sort(n_levels)

  # Compute b_u
  multiplier <- gmp::pow.bigz(n_levels, seq_along(n_levels) - 1)

  b_u <-
    gmp::`%*%`(data_wo_t, multiplier) %>%
    as.vector()


  # Compute b_u+
  b_u_plus <-
    gmp::mul.bigz(b_u, as.bigz(n_levels)) %>%
    gmp::add.bigz(data$treated) %>%
    as.vector()

  # Compute c_u
  c_u <- aggregate_table(b_u)

  # Compute c_u+
  c_u_plus <- aggregate_table(b_u_plus)

  match_index <- (c_u != c_u_plus) & (c_u >= 2)

  index <- b_u[match_index]

  return(list(match_index = match_index,
              index = index))
}
