# shared helpers for sc_slic workflows

# prepare and validate slic arguments
.sc_slic_prep_args = function(x, step, step_unit, compactness, dist_fun, avg_fun, clean, minarea, iter,
                              k, centers, outcomes, chunks, iter_diagnostics, verbose) {
  # Validate core arguments and types
  .sc_slic_validate_args(step, step_unit, compactness, k, centers, chunks, dist_fun, avg_fun, iter, minarea)
  outcomes = .sc_slic_prep_outcomes(outcomes)
  # Normalize input to SpatRaster
  x = .sc_util_prep_raster(x)
  if (terra::is.lonlat(x)) {
    warning("The input raster uses a geographic (lon/lat) CRS; consider projecting it before using SLIC", call. = FALSE)
  }
  # Resolve step from k when needed
  step = .sc_slic_prep_step(x, step, k, step_unit)
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
  chunk_ext = .sc_chunk_extents(dim(x), limit = chunks, step = step)
  # Disable iter diagnostics when chunking is active
  if (iter_diagnostics && nrow(chunk_ext) > 1) {
    warning("Iteration diagnostics are only available when chunks = FALSE (single chunk). Iteration diagnostics were disabled.", call. = FALSE)
    iter_diagnostics = FALSE
  }
  verbose_cpp = if (is.numeric(verbose) && length(verbose) == 1 && !is.na(verbose) && verbose >= 2) verbose else 0
  adaptive_compactness = is.character(compactness) &&
    length(compactness) == 1 && !is.na(compactness) && compactness == "auto"
  if (adaptive_compactness) {
    compactness = 0
  }
  # Package prep results for downstream functions
  return(list(x = x, step = step, step_unit = step_unit,
              dist_fun_input = dist_fun,
              input_centers = input_centers, funs = funs,
              minarea = minarea, chunk_ext = chunk_ext,
              iter_diagnostics = iter_diagnostics, outcomes = outcomes,
              compactness = compactness, adaptive_compactness = adaptive_compactness,
              clean = clean, iter = iter,
              verbose = verbose, verbose_cpp = verbose_cpp))
}

# validate slic arguments and types
.sc_slic_validate_args = function(step, step_unit, compactness, k, centers, chunks, dist_fun, avg_fun, iter, minarea) {
  if (!missing(step_unit)) {
    if (!is.character(step_unit) || length(step_unit) != 1 || is.na(step_unit) ||
        !(step_unit %in% c("cells", "map"))) {
      stop("The 'step_unit' argument must be 'cells' or 'map'", call. = FALSE)
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
  if (!is.numeric(compactness) && !(is.character(compactness) && length(compactness) == 1 && !is.na(compactness) && compactness == "auto")) {
    stop("The 'compactness' argument must be numeric or 'auto'", call. = FALSE)
  }
  if (!missing(minarea) && (!is.numeric(minarea) || length(minarea) != 1 || is.na(minarea) || minarea < 0)) {
    stop("The 'minarea' argument must be a non-negative numeric value", call. = FALSE)
  }
  return(invisible(TRUE))
}

# derive step from k when needed
.sc_slic_prep_step = function(x, step, k, step_unit = "cells") {
  if (!is.null(step)) {
    if (step_unit == "map") {
      res = terra::res(x)
      if (!isTRUE(all.equal(res[[1]], res[[2]]))) {
        stop("Map-unit step requires square cells; res(x) has different x/y resolution", call. = FALSE)
      }
      step = round(step / res[[1]])
      if (step < 1) {
        step = 1
      }
    }
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

# resolve outcomes for returned fields
.sc_slic_prep_outcomes = function(outcomes) {
  allowed = c("supercells", "coordinates", "values")
  if (is.null(outcomes) || !is.character(outcomes) || anyNA(outcomes)) {
    stop("The 'outcomes' argument must be a character vector", call. = FALSE)
  }
  if (length(outcomes) == 0) {
    return(character(0))
  }
  outcomes = unique(outcomes)
  bad = setdiff(outcomes, allowed)
  if (length(bad) > 0) {
    stop("The 'outcomes' argument must be one or more of: ",
         paste(allowed, collapse = ", "), call. = FALSE)
  }
  outcomes
}

# select requested outcomes columns from sf output
.sc_slic_select_outcomes = function(x, outcomes) {
  geom_col = attr(x, "sf_column")
  value_cols = setdiff(names(x), c("supercells", "x", "y", geom_col))
  keep = character(0)
  if ("supercells" %in% outcomes && "supercells" %in% names(x)) {
    keep = c(keep, "supercells")
  }
  if ("coordinates" %in% outcomes) {
    if ("x" %in% names(x)) keep = c(keep, "x")
    if ("y" %in% names(x)) keep = c(keep, "y")
  }
  if ("values" %in% outcomes) {
    keep = c(keep, value_cols)
  }
  keep = unique(c(keep, geom_col))
  x[, keep, drop = FALSE]
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

  slic_sf = .sc_slic_select_outcomes(slic_sf, prep$outcomes)

  attr(slic_sf, "step") = prep$step
  attr(slic_sf, "step_unit") = prep$step_unit
  attr(slic_sf, "compactness") = prep$compactness
  attr(slic_sf, "dist_fun") = prep$dist_fun_input
  attr(slic_sf, "method") = if (isTRUE(prep$adaptive_compactness)) "slic0" else "slic"
  class(slic_sf) = unique(c("supercells", class(slic_sf)))
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
    adaptive_compactness = prep$adaptive_compactness,
    dist_name = prep$funs$dist_name,
    dist_fun = prep$funs$dist_fun,
    avg_fun_fun = prep$funs$avg_fun_fun,
    avg_fun_name = prep$funs$avg_fun_name,
    clean = prep$clean,
    iter = prep$iter,
    minarea = prep$minarea,
    input_centers = prep$input_centers,
    verbose = prep$verbose_cpp,
    iter_diagnostics = prep$iter_diagnostics
  )
  if (nrow(prep$chunk_ext) == 1) {
    list(chunks = list(do.call(single_runner, args)))
  } else {
    chunks = .sc_slic_apply_chunks(prep$chunk_ext, chunk_runner, args)
    list(chunks = chunks)
  }
}
