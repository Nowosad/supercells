# shared helpers for sc_slic workflows

# prepare and validate slic arguments
.sc_slic_prep_args = function(x, step, compactness, dist_fun, avg_fun, clean, minarea, iter,
                              k, centers, metadata, chunks, iter_diagnostics, verbose) {
  # Validate core arguments and types
  .sc_slic_validate_args(step, compactness, k, centers, chunks, dist_fun, avg_fun, iter, metadata, minarea)
  # Normalize input to SpatRaster
  x = .sc_util_prep_raster(x)
  # Resolve step from k when needed
  step = .sc_slic_prep_step(x, step, k)
  # Adjust numeric chunks to match step size
  if (is.numeric(chunks)) {
    if (chunks < step) {
      stop("The 'chunks' argument must be >= 'step' when numeric", call. = FALSE)
    }
    if (chunks %% step != 0) {
      new_chunks = ceiling(chunks / step) * step
      warning(sprintf("The 'chunks' argument was rounded up to %s to match step %s", new_chunks, step), call. = FALSE)
      chunks = new_chunks
    }
  }
  # Prepare initial centers (either custom or placeholder)
  input_centers = .sc_slic_prep_centers(x, centers)
  # Resolve distance/averaging functions and names
  funs = .sc_slic_prep_funs(dist_fun, avg_fun)
  # Validate and normalize minarea
  minarea = .sc_slic_prep_minarea(minarea, step)
  # Compute chunk extents based on size/limits
  chunk_ext = .sc_chunk_extents(dim(x), limit = chunks)
  # Disable iter diagnostics when chunking is active
  if (iter_diagnostics && nrow(chunk_ext) > 1) {
    warning("Iteration diagnostics are only available when chunks = FALSE (single chunk). Iteration diagnostics were disabled.", call. = FALSE)
    iter_diagnostics = FALSE
  }
  # Package prep results for downstream functions
  return(list(x = x, step = step, input_centers = input_centers, funs = funs,
              minarea = minarea, chunk_ext = chunk_ext,
              iter_diagnostics = iter_diagnostics, metadata = metadata,
              compactness = compactness, clean = clean, iter = iter,
              verbose = verbose))
}

# validate slic arguments and types
.sc_slic_validate_args = function(step, compactness, k, centers, chunks, dist_fun, avg_fun, iter, metadata, minarea) {
  if (!missing(metadata)) {
    if (!is.logical(metadata) || length(metadata) != 1 || is.na(metadata)) {
      stop("The 'metadata' argument must be TRUE or FALSE", call. = FALSE)
    }
  }
  if (!is.numeric(iter) || length(iter) != 1 || is.na(iter) || iter < 0) {
    stop("The 'iter' argument must be a non-negative numeric value", call. = FALSE)
  }
  if (!is.character(avg_fun) && !is.function(avg_fun)) {
    stop("The 'avg_fun' argument must be a string or a function", call. = FALSE)
  }
  if (is.character(avg_fun) && !(avg_fun %in% c("mean", "median"))) {
    stop("The 'avg_fun' argument must be 'mean', 'median', or a custom function", call. = FALSE)
  }
  if (!is.character(dist_fun) && !is.function(dist_fun)) {
    stop("The 'dist_fun' argument must be a string or a function", call. = FALSE)
  }
  if (!is.logical(chunks) && !is.numeric(chunks)) {
    stop("The 'chunks' argument must be FALSE, TRUE, or a numeric value", call. = FALSE)
  }
  if (is.numeric(chunks) && (length(chunks) != 1 || is.na(chunks) || chunks <= 0)) {
    stop("The 'chunks' argument must be a positive number", call. = FALSE)
  }
  if (!is.null(centers)) {
    if (!inherits(centers, "sf")) {
      stop("The 'centers' argument must be an sf object", call. = FALSE)
    }
    if (is.null(step)) {
      stop("The 'step' argument is required when 'centers' is provided", call. = FALSE)
    }
    if (!is.null(k)) {
      stop("Use either 'k' or 'centers', not both", call. = FALSE)
    }
    if (chunks > 0) {
      stop("Chunks cannot be used for custom cluster centers!", call. = FALSE)
    }
  }
  if (is.null(step) && is.null(k) && is.null(centers)) {
    stop("You need to specify either 'k' or 'step'", call. = FALSE)
  }
  if (!is.null(step) && !is.null(k)) {
    stop("You can specify either 'k' or 'step', not both", call. = FALSE)
  }
  if (missing(compactness)) {
    stop("The 'compactness' argument is required", call. = FALSE)
  }
  if (!missing(minarea) && (!is.numeric(minarea) || length(minarea) != 1 || is.na(minarea) || minarea < 0)) {
    stop("The 'minarea' argument must be a non-negative numeric value", call. = FALSE)
  }
  return(invisible(TRUE))
}

# derive step from k when needed
.sc_slic_prep_step = function(x, step, k) {
  if (!is.null(step)) {
    return(step)
  }
  mat = dim(x)[1:2]
  superpixelsize = round((mat[1] * mat[2]) / k + 0.5)
  return(round(sqrt(superpixelsize) + 0.5))
}

# normalize custom centers or create placeholder
.sc_slic_prep_centers = function(x, centers) {
  if (is.null(centers)) {
    return(matrix(c(0L, 0L), ncol = 2))
  }
  return(.sc_util_centers_to_dims(x, centers))
}

# resolve distance and averaging functions
.sc_slic_prep_funs = function(dist_fun, avg_fun) {
  avg_prep = .sc_util_prep_avg_fun(avg_fun)
  dist_prep = .sc_util_prep_dist_fun(dist_fun)
  return(list(dist_name = dist_prep$dist_name, dist_fun = dist_prep$dist_fun,
              avg_fun_name = avg_prep$avg_fun_name, avg_fun_fun = avg_prep$avg_fun_fun))
}

# validate and normalize minarea
.sc_slic_prep_minarea = function(minarea, step) {
  if (missing(minarea)) {
    return(0)
  }
  if (minarea > step^2) {
    warning("The provided minarea value is larger than than the average supercell (step^2). The connectivity cleaning is likely to fail.", call. = FALSE)
  }
  return(minarea)
}

# apply slic over chunks
.sc_slic_apply_chunks = function(chunk_ext, fun, args) {
  return(do.call(apply, c(list(chunk_ext, MARGIN = 1, FUN = fun), args)))
}

# add iter diagnostics attribute when enabled
.sc_slic_add_iter_attr = function(chunks, iter_diagnostics) {
  if (isTRUE(iter_diagnostics) && length(chunks) > 0) {
    return(attr(chunks[[1]], "iter_diagnostics"))
  }
  NULL
}

# finalize slic output with ids, metadata, and attributes
.sc_slic_post = function(chunks, prep, iter_attr) {

  slic_sf = .sc_chunk_update_ids(chunks)

  if (!isTRUE(prep$metadata)) {
    remove_cols = intersect(names(slic_sf), c("supercells", "x", "y"))
    if (length(remove_cols) > 0) {
      slic_sf = slic_sf[, -which(names(slic_sf) %in% remove_cols)]
    }
  }

  attr(slic_sf, "step") = prep$step
  attr(slic_sf, "compactness") = prep$compactness
  attr(slic_sf, "method") = "slic"
  class(slic_sf) = c(class(slic_sf), "supercells")
  if (!is.null(iter_attr)) {
    attr(slic_sf, "iter_diagnostics") = iter_attr
  }
  return(slic_sf)
}

# dispatch slic run for single or chunked input
.sc_slic_segment = function(prep, single_runner, chunk_runner) {
  args = list(
    x = prep$x,
    step = prep$step,
    compactness = prep$compactness,
    dist_name = prep$funs$dist_name,
    dist_fun = prep$funs$dist_fun,
    avg_fun_fun = prep$funs$avg_fun_fun,
    avg_fun_name = prep$funs$avg_fun_name,
    clean = prep$clean,
    iter = prep$iter,
    minarea = prep$minarea,
    input_centers = prep$input_centers,
    verbose = prep$verbose,
    iter_diagnostics = prep$iter_diagnostics,
    metadata = prep$metadata
  )
  if (nrow(prep$chunk_ext) == 1) {
    list(chunks = list(do.call(single_runner, args)))
  } else {
    chunks = .sc_slic_apply_chunks(prep$chunk_ext, chunk_runner, args)
    list(chunks = chunks)
  }
}
