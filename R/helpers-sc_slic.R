# Shared helpers for sc_slic* workflows

.sc_slic_prep_args = function(x, step, compactness, k, centers, dist_fun, avg_fun,
                         minarea, chunks, iter, transform, metadata, iter_diagnostics) {
  # Validate core arguments and types
  .sc_slic_validate_args(step, compactness, k, centers, chunks, dist_fun, avg_fun, iter, transform, metadata, minarea)
  # Normalize input to SpatRaster
  x = .sc_prep_raster(x)
  # Resolve step from k when needed
  step = .sc_slic_prep_step(x, step, k)
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
  chunk_ext = prep_chunks_ext(dim(x), limit = chunks)
  # Disable iter diagnostics when chunking is active
  if (iter_diagnostics && nrow(chunk_ext) > 1) {
    warning("Iteration diagnostics are only available when chunks = FALSE (single chunk). Iteration diagnostics were disabled.", call. = FALSE)
    iter_diagnostics = FALSE
  }
  # Package prep results for downstream functions
  return(list(x = x, step = step, input_centers = input_centers, funs = funs,
              minarea = minarea, chunk_ext = chunk_ext, iter_diagnostics = iter_diagnostics))
}

.sc_slic_validate_args = function(step, compactness, k, centers, chunks, dist_fun, avg_fun, iter, transform, metadata, minarea) {
  if (!missing(metadata)) {
    if (!is.logical(metadata) || length(metadata) != 1 || is.na(metadata)) {
      stop("The 'metadata' argument must be TRUE or FALSE", call. = FALSE)
    }
  }
  if (!is.numeric(iter) || length(iter) != 1 || is.na(iter) || iter < 0) {
    stop("The 'iter' argument must be a non-negative numeric value", call. = FALSE)
  }
  if (!is.null(transform) && !identical(transform, "to_LAB")) {
    stop("The 'transform' argument must be NULL or 'to_LAB'", call. = FALSE)
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

.sc_slic_prep_step = function(x, step, k) {
  if (!is.null(step)) {
    return(step)
  }
  mat = dim(x)[1:2]
  superpixelsize = round((mat[1] * mat[2]) / k + 0.5)
  return(round(sqrt(superpixelsize) + 0.5))
}

.sc_slic_prep_centers = function(x, centers) {
  if (is.null(centers)) {
    return(matrix(c(0L, 0L), ncol = 2))
  }
  return(centers_to_dims(x, centers))
}

.sc_slic_prep_funs = function(dist_fun, avg_fun) {
  avg_prep = .sc_prep_avg_fun(avg_fun)
  dist_prep = .sc_prep_dist_fun(dist_fun)
  return(list(dist_name = dist_prep$dist_name, dist_fun = dist_prep$dist_fun,
              avg_fun_name = avg_prep$avg_fun_name, avg_fun_fun = avg_prep$avg_fun_fun))
}

.sc_slic_prep_minarea = function(minarea, step) {
  if (missing(minarea)) {
    return(0)
  }
  if (minarea > step^2) {
    warning("The provided minarea value is larger than than the average supercell (step^2). The connectivity cleaning is likely to fail.", call. = FALSE)
  }
  return(minarea)
}

.sc_slic_apply_chunks = function(prep, fun, compactness, clean, iter, transform, verbose, future, metadata = NULL) {
  args = list(
    x = prep$x,
    step = prep$step,
    compactness = compactness,
    dist_name = prep$funs$dist_name,
    dist_fun = prep$funs$dist_fun,
    avg_fun_fun = prep$funs$avg_fun_fun,
    avg_fun_name = prep$funs$avg_fun_name,
    clean = clean,
    iter = iter,
    minarea = prep$minarea,
    transform = transform,
    input_centers = prep$input_centers,
    verbose = verbose,
    iter_diagnostics = prep$iter_diagnostics
  )
  if (!is.null(metadata)) {
    args$metadata = metadata
  }

  if (future) {
    if (in_memory(prep$x)) {
      names_x = names(prep$x)
      prep$x = terra::writeRaster(prep$x, tempfile(fileext = ".tif"))
      names(prep$x) = names_x
    }
    if (!in_memory(prep$x)) {
      prep$x = terra::sources(prep$x)[[1]]
    }
    args$x = prep$x
    oopts = options(future.globals.maxSize = +Inf)
    on.exit(options(oopts))
    return(do.call(future.apply::future_apply, c(list(prep$chunk_ext, MARGIN = 1, FUN = fun), args, list(future.seed = TRUE))))
  } else {
    return(do.call(apply, c(list(prep$chunk_ext, MARGIN = 1, FUN = fun), args)))
  }
}

.sc_slic_post = function(slic_sf, metadata, step, compactness, iter_attr) {

  slic_sf = update_supercells_ids(slic_sf)

  if (!isTRUE(metadata)) {
    remove_cols = intersect(names(slic_sf), c("supercells", "x", "y"))
    if (length(remove_cols) > 0) {
      slic_sf = slic_sf[, -which(names(slic_sf) %in% remove_cols)]
    }
  }

  attr(slic_sf, "step") = step
  attr(slic_sf, "compactness") = compactness
  attr(slic_sf, "method") = "slic"
  class(slic_sf) = c(class(slic_sf), "supercells")
  if (!is.null(iter_attr)) {
    attr(slic_sf, "iter_diagnostics") = iter_attr
  }
  return(slic_sf)
}
