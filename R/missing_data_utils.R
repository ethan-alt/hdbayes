#' Utility functions for handling missing covariates using Gaussian copula
#'
#' @name missing_data_utils
#' @keywords internal
NULL


#' Detect missing data patterns in a dataset
#'
#' @param data a data.frame that may contain missing values
#' @param formula a two-sided formula specifying the model
#' @return a list containing:
#'   \itemize{
#'     \item has_missing: logical indicating if there are missing values in covariates
#'     \item response_name: name of the response variable
#'     \item covariate_names: names of covariates in the formula
#'     \item miss_idx: matrix with columns (row_index, col_index) for each missing value
#'     \item n_miss: total number of missing covariate values
#'     \item complete_cases: logical vector indicating complete cases
#'     \item missing_covs: names of covariates that have missing values
#'   }
#' @noRd
detect.missing.pattern = function(data, formula) {
  # Get variable names from formula (use all.vars to get raw variable names)
  all_varnames    = all.vars(formula)
  response_name   = all_varnames[1]
  covariate_names = all_varnames[-1]

  # Check response - cannot have missing values
  if (any(is.na(data[[response_name]]))) {
    stop("Response variable '", response_name, "' cannot have missing values")
  }

  # Identify which covariates have missing values
  cov_data = data[, covariate_names, drop = FALSE]
  missing_mask = is.na(cov_data)

  # Find covariates with at least one missing value
  missing_covs = covariate_names[colSums(missing_mask) > 0]

  # Create matrix of missing indices: (row, col) pairs
  if (length(missing_covs) > 0) {
    miss_locations = which(missing_mask, arr.ind = TRUE)
    # Convert column indices to covariate indices within missing_covs
    col_in_missing = match(covariate_names[miss_locations[, 2]], missing_covs)
    miss_idx = cbind(miss_locations[, 1], col_in_missing)
    colnames(miss_idx) = c("row", "col")
  } else {
    miss_idx = matrix(integer(0), ncol = 2)
    colnames(miss_idx) = c("row", "col")
  }

  # Identify complete cases (no missing covariates)
  complete_cases = rowSums(missing_mask) == 0

  list(
    has_missing     = length(missing_covs) > 0,
    response_name   = response_name,
    covariate_names = covariate_names,
    miss_idx        = miss_idx,
    n_miss          = nrow(miss_idx),
    complete_cases  = complete_cases,
    missing_covs    = missing_covs,
    n_complete      = sum(complete_cases)
  )
}


#' Detect covariate types (continuous, binary, ordinal, count)
#'
#' @param x a vector of covariate values (may contain NA)
#' @return a list with:
#'   \itemize{
#'     \item type: character string ("continuous", "binary", "ordinal", "count")
#'     \item type_id: integer (1=continuous, 2=binary, 3=ordinal, 4=count)
#'   }
#' @noRd
detect.covariate.type = function(x) {
  x_obs = x[!is.na(x)]

  # Check if binary (only 0 and 1)
  unique_vals = sort(unique(x_obs))
  if (length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))) {
    return(list(type = "binary", type_id = 2L))
  }

  # Check if factor/ordered (treat as ordinal)
  if (is.factor(x) || is.ordered(x)) {
    return(list(type = "ordinal", type_id = 3L))
  }

  # Check if integer with small range (likely count or ordinal)
  if (is.integer(x_obs) || all(x_obs == floor(x_obs))) {
    # If all non-negative integers with limited range, could be count
    if (all(x_obs >= 0) && max(x_obs) - min(x_obs) < 20) {
      n_unique = length(unique_vals)
      if (n_unique <= 10) {
        return(list(type = "ordinal", type_id = 3L))
      }
      return(list(type = "count", type_id = 4L))
    }
  }

  # Default to continuous
  return(list(type = "continuous", type_id = 1L))
}


#' Estimate marginal distribution parameters for a covariate
#'
#' @param x a vector of covariate values (may contain NA)
#' @param type_id integer indicating type (1=continuous, 2=binary, 3=ordinal, 4=count)
#' @return a list of marginal parameters appropriate for the type
#' @noRd
estimate.marginal.params = function(x, type_id) {
  x_obs = x[!is.na(x)]

  if (type_id == 1L) {
    # Continuous: estimate mean and sd
    return(list(
      mean = mean(x_obs),
      sd   = sd(x_obs)
    ))
  } else if (type_id == 2L) {
    # Binary: estimate probability of 1
    return(list(
      prob = mean(x_obs)
    ))
  } else if (type_id == 3L) {
    # Ordinal: estimate cutpoints for latent normal
    # Use cumulative probabilities to determine thresholds
    unique_vals = sort(unique(x_obs))
    n_levels = length(unique_vals)
    cum_probs = cumsum(table(x_obs)) / length(x_obs)
    # Convert to latent normal scale via qnorm
    # Don't include the last cutpoint (infinity)
    cutpoints = stats::qnorm(cum_probs[-n_levels])
    return(list(
      levels    = unique_vals,
      cutpoints = cutpoints,
      n_levels  = n_levels
    ))
  } else if (type_id == 4L) {
    # Count: estimate Poisson rate (lambda)
    return(list(
      lambda = mean(x_obs)
    ))
  }
}


#' Estimate copula correlation matrix from complete cases
#'
#' Transforms covariates to normal scale using probability integral transform
#' and estimates the correlation matrix.
#'
#' @param data a data.frame
#' @param covariate_names character vector of covariate names to include
#' @param marginal_types list of type_ids for each covariate
#' @param marginal_params list of marginal parameters for each covariate
#' @return a correlation matrix
#' @noRd
estimate.copula.corr = function(data, covariate_names, marginal_types, marginal_params) {
  n_covs = length(covariate_names)

  # Get complete cases only
  cov_data = data[, covariate_names, drop = FALSE]
  complete_idx = complete.cases(cov_data)

  if (sum(complete_idx) < 3) {
    stop("Need at least 3 complete cases to estimate copula correlation")
  }

  cov_data_complete = cov_data[complete_idx, , drop = FALSE]
  n_complete = nrow(cov_data_complete)

  # Transform each covariate to normal scale
  Z = matrix(NA_real_, nrow = n_complete, ncol = n_covs)

  for (j in seq_len(n_covs)) {
    x = cov_data_complete[[j]]
    type_id = marginal_types[[j]]
    params = marginal_params[[j]]

    Z[, j] = transform.to.normal(x, type_id, params)
  }

  # Compute correlation matrix (handle edge cases)
  if (n_covs == 1) {
    corr_mat = matrix(1, 1, 1)
  } else {
    corr_mat = stats::cor(Z, use = "complete.obs")
    # Ensure positive definiteness (can fail with small samples)
    corr_mat = make.positive.definite(corr_mat)
  }

  return(corr_mat)
}


#' Transform covariate to standard normal scale
#'
#' @param x vector of observed values
#' @param type_id integer indicating type
#' @param params list of marginal parameters
#' @return vector of transformed values on standard normal scale
#' @noRd
transform.to.normal = function(x, type_id, params) {
  if (type_id == 1L) {
    # Continuous: standardize
    z = (x - params$mean) / params$sd
  } else if (type_id == 2L) {
    # Binary: use probit (inverse normal CDF)
    # Add jitter to avoid infinite values
    p = params$prob
    # Map 0 -> qnorm(p/2), 1 -> qnorm(1 - (1-p)/2) roughly
    z = ifelse(x == 1,
               stats::qnorm(p + (1 - p) * 0.5),
               stats::qnorm(p * 0.5))
  } else if (type_id == 3L) {
    # Ordinal: map to midpoint of latent interval
    n_levels = params$n_levels
    cutpoints = c(-Inf, params$cutpoints, Inf)
    level_idx = match(x, params$levels)
    # Midpoint of interval
    z = (cutpoints[level_idx] + cutpoints[level_idx + 1]) / 2
    # Handle infinities
    z[is.infinite(z) & z < 0] = cutpoints[2] - 1
    z[is.infinite(z) & z > 0] = cutpoints[n_levels] + 1
  } else if (type_id == 4L) {
    # Count: use inverse Poisson CDF
    lambda = params$lambda
    p = stats::ppois(x, lambda)
    # Avoid exact 0 or 1
    p = pmin(pmax(p, 1e-6), 1 - 1e-6)
    z = stats::qnorm(p)
  }
  return(z)
}


#' Ensure a matrix is positive definite
#'
#' @param mat a square symmetric matrix
#' @param tol tolerance for eigenvalues
#' @return a positive definite matrix close to the input
#' @noRd
make.positive.definite = function(mat, tol = 1e-6) {
  eig = eigen(mat, symmetric = TRUE)
  eig$values = pmax(eig$values, tol)
  pd_mat = eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  # Re-normalize to correlation matrix
  d = sqrt(diag(pd_mat))
  pd_mat = pd_mat / outer(d, d)
  return(pd_mat)
}


#' Prepare copula data for all datasets in data.list
#'
#' @param formula a two-sided formula
#' @param data.list a list of data.frames (first is current, rest are historical)
#' @param copula.corr optional user-provided correlation matrix
#' @return a list containing all information needed for Stan
#' @noRd
prepare.copula.data = function(formula, data.list, copula.corr = NULL) {
  n_datasets = length(data.list)

  # Detect missing patterns for each dataset
  miss_patterns = lapply(data.list, detect.missing.pattern, formula = formula)

  # Identify all covariates that have missing values across any dataset
  all_missing_covs = unique(unlist(lapply(miss_patterns, function(x) x$missing_covs)))

  if (length(all_missing_covs) == 0) {
    # No missing data - return NULL to signal use of standard methods
    return(NULL)
  }

  # Get all covariate names (should be same across datasets)
  covariate_names = miss_patterns[[1]]$covariate_names
  response_name = miss_patterns[[1]]$response_name

  # Determine types and estimate marginals from pooled complete data
  # Pool all datasets for parameter estimation
  pooled_data = do.call(rbind, data.list)

  marginal_types = list()
  marginal_params = list()

  for (cov in all_missing_covs) {
    type_info = detect.covariate.type(pooled_data[[cov]])
    marginal_types[[cov]] = type_info$type_id
    marginal_params[[cov]] = estimate.marginal.params(pooled_data[[cov]], type_info$type_id)
  }

  # Estimate copula correlation if not provided
  if (is.null(copula.corr)) {
    copula.corr = estimate.copula.corr(
      pooled_data, all_missing_covs, marginal_types, marginal_params
    )
  } else {
    # Validate user-provided correlation matrix
    q = length(all_missing_covs)
    if (!is.matrix(copula.corr) || nrow(copula.corr) != q || ncol(copula.corr) != q) {
      stop("copula.corr must be a ", q, " x ", q, " correlation matrix")
    }
  }

  list(
    has_missing      = TRUE,
    n_datasets       = n_datasets,
    miss_patterns    = miss_patterns,
    all_missing_covs = all_missing_covs,
    covariate_names  = covariate_names,
    response_name    = response_name,
    marginal_types   = marginal_types,
    marginal_params  = marginal_params,
    copula_corr      = copula.corr
  )
}


#' Stack data for Stan with missing covariate handling
#'
#' Similar to stack.data but handles missing values specially
#'
#' @param formula a two-sided formula
#' @param data.list a list of data.frames
#' @param copula_info output from prepare.copula.data
#' @return a list with stacked data suitable for Stan
#' @noRd
stack.data.mc = function(formula, data.list, copula_info) {
  n_datasets = length(data.list)
  all_missing_covs = copula_info$all_missing_covs
  q = length(all_missing_covs)

  # Get response for all datasets
  y_list = lapply(data.list, function(d) d[[copula_info$response_name]])
  y = unlist(y_list)
  N = length(y)

  # Build design matrix, replacing missing values with 0 (placeholder)
  # We'll track which are missing via miss_idx
  X_list = list()
  miss_idx_list = list()
  row_offset = 0

  for (k in seq_len(n_datasets)) {
    d = data.list[[k]]
    # Create a copy with missing values filled with column means (for model.matrix)
    d_filled = d
    for (cov in all_missing_covs) {
      if (any(is.na(d[[cov]]))) {
        d_filled[[cov]][is.na(d[[cov]])] = mean(d[[cov]], na.rm = TRUE)
      }
    }
    X_k = stats::model.matrix(formula, d_filled)

    # Track missing indices (in terms of design matrix columns)
    miss_info = copula_info$miss_patterns[[k]]
    if (miss_info$n_miss > 0) {
      # Map covariate names to design matrix column indices
      # This is tricky because factors expand to multiple columns
      dm_colnames = colnames(X_k)
      for (m in seq_len(nrow(miss_info$miss_idx))) {
        row_in_data = miss_info$miss_idx[m, "row"]
        cov_idx = miss_info$miss_idx[m, "col"]
        cov_name = miss_info$missing_covs[cov_idx]

        # Find which design matrix columns correspond to this covariate
        # Need to handle transformations like scale(age) -> age
        # Simple case: covariate name matches column name exactly
        dm_cols = which(dm_colnames == cov_name)
        if (length(dm_cols) == 0) {
          # Check for transformations: columns containing the covariate name
          # e.g., "scale(age)" contains "age"
          dm_cols = grep(cov_name, dm_colnames, fixed = TRUE)
        }
        if (length(dm_cols) == 0) {
          # Check for factor expansion (columns starting with covariate name)
          dm_cols = grep(paste0("^", cov_name), dm_colnames)
        }

        for (dm_col in dm_cols) {
          miss_idx_list[[length(miss_idx_list) + 1]] = c(
            row_offset + row_in_data,  # global row index
            dm_col,                     # design matrix column
            which(all_missing_covs == cov_name)  # index in missing covs
          )
        }
      }
    }

    X_list[[k]] = X_k
    row_offset = row_offset + nrow(X_k)
  }

  X = do.call(rbind, X_list)
  rownames(X) = NULL

  # Convert miss_idx_list to matrix
  if (length(miss_idx_list) > 0) {
    miss_idx = do.call(rbind, miss_idx_list)
    colnames(miss_idx) = c("row", "dm_col", "cov_idx")
  } else {
    miss_idx = matrix(integer(0), ncol = 3)
    colnames(miss_idx) = c("row", "dm_col", "cov_idx")
  }

  # Dataset boundaries
  num_obs = sapply(data.list, nrow)
  end_index = cumsum(num_obs)
  start_index = c(1, end_index[-n_datasets] + 1)

  list(
    X           = X,
    y           = y,
    N           = N,
    p           = ncol(X),
    K           = n_datasets,
    start_index = start_index,
    end_index   = end_index,
    miss_idx    = miss_idx,
    n_miss      = nrow(miss_idx),
    q           = q
  )
}


#' Prepare marginal parameters for Stan
#'
#' Convert R marginal parameters to Stan-ready format
#'
#' @param marginal_types list of type_ids
#' @param marginal_params list of marginal parameter lists
#' @param all_missing_covs character vector of covariate names
#' @return a list ready for Stan data block
#' @noRd
prepare.marginal.stan = function(marginal_types, marginal_params, all_missing_covs) {
  q = length(all_missing_covs)

  type_ids = integer(q)
  means = numeric(q)
  sds = numeric(q)
  # For ordinal: max of 10 levels supported
  max_levels = 10
  n_levels = integer(q)
  cutpoints = matrix(-Inf, nrow = q, ncol = max_levels - 1)

  for (j in seq_len(q)) {
    cov = all_missing_covs[j]
    type_id = marginal_types[[cov]]
    params = marginal_params[[cov]]

    type_ids[j] = type_id

    if (type_id == 1L) {
      # Continuous
      means[j] = params$mean
      sds[j] = params$sd
    } else if (type_id == 2L) {
      # Binary: store prob as "mean", sd = 1 (not used)
      means[j] = params$prob
      sds[j] = 1.0
    } else if (type_id == 3L) {
      # Ordinal
      means[j] = 0.0
      sds[j] = 1.0
      n_levels[j] = params$n_levels
      n_cuts = length(params$cutpoints)
      cutpoints[j, 1:n_cuts] = params$cutpoints
    } else if (type_id == 4L) {
      # Count
      means[j] = params$lambda
      sds[j] = 1.0
    }
  }

  list(
    marginal_type = type_ids,
    marginal_mean = means,
    marginal_sd   = sds,
    n_levels      = n_levels,
    cutpoints     = cutpoints,
    max_levels    = max_levels
  )
}
