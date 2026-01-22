# Shared helpers for sc_slic* workflows.

.sc_slic_prep = function(x, k, step, centers, compactness, dist_fun, avg_fun,
                         minarea, chunks, iter, transform, metadata, iter_diagnostics) {
  x = .sc_slic_input(x)
  .sc_slic_validate(k, step, centers, compactness, chunks, dist_fun, avg_fun, iter, transform, metadata)
  step = .sc_slic_step(x, k, step)
  input_centers = .sc_slic_centers(x, centers)
  funs = .sc_slic_funs(dist_fun, avg_fun)
  minarea = .sc_slic_minarea(minarea, step)
  chunk_ext = prep_chunks_ext(dim(x), limit = chunks)
  if (iter_diagnostics && nrow(chunk_ext) > 1) {
    warning("Iteration diagnostics are only available when chunks = FALSE (single chunk). Iteration diagnostics were disabled.", call. = FALSE)
    iter_diagnostics = FALSE
  }
  list(x = x, step = step, input_centers = input_centers, funs = funs,
       minarea = minarea, chunk_ext = chunk_ext, iter_diagnostics = iter_diagnostics)
}

.sc_slic_apply_chunks = function(prep, fun, compactness, clean, iter, transform, verbose, future) {
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
    do.call(future.apply::future_apply, c(list(prep$chunk_ext, MARGIN = 1, FUN = fun), args, list(future.seed = TRUE)))
  } else {
    do.call(apply, c(list(prep$chunk_ext, MARGIN = 1, FUN = fun), args))
  }
}

.sc_slic_input = function(x) {
  if (inherits(x, "SpatRaster")) {
    return(x)
  }
  if (inherits(x, "stars")) {
    return(terra::rast(x))
  }
  stop("The SpatRaster class is expected as an input", call. = FALSE)
}

.sc_slic_validate = function(k, step, centers, compactness, chunks, dist_fun, avg_fun, iter, transform, metadata) {
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
}

.sc_slic_step = function(x, k, step) {
  if (!is.null(step)) {
    return(step)
  }
  mat = dim(x)[1:2]
  superpixelsize = round((mat[1] * mat[2]) / k + 0.5)
  round(sqrt(superpixelsize) + 0.5)
}

.sc_slic_centers = function(x, centers) {
  if (is.null(centers)) {
    return(matrix(c(0L, 0L), ncol = 2))
  }
  centers_to_dims(x, centers)
}

.sc_slic_funs = function(dist_fun, avg_fun) {
  if (is.character(avg_fun)) {
    avg_fun_name = avg_fun
    avg_fun_fun = function() ""
  } else {
    avg_fun_name = ""
    avg_fun_fun = avg_fun
  }
  if (is.character(dist_fun)) {
    if (!(dist_fun %in% c("euclidean", "jsd", "dtw", "dtw2d", philentropy::getDistMethods()))) {
      stop("The provided distance function ('dist_fun') does not exist!", call. = FALSE)
    }
    dist_name = dist_fun
    dist_fun = function() ""
  } else {
    dist_name = ""
  }
  list(dist_name = dist_name, dist_fun = dist_fun,
       avg_fun_name = avg_fun_name, avg_fun_fun = avg_fun_fun)
}

.sc_slic_minarea = function(minarea, step) {
  if (missing(minarea)) {
    return(0)
  }
  if (!is.numeric(minarea) || length(minarea) != 1 || is.na(minarea) || minarea < 0) {
    stop("The 'minarea' argument must be a non-negative numeric value", call. = FALSE)
  }
  if (minarea > step^2) {
    warning("The provided minarea value is larger than than the average supercell (step^2). The connectivity cleaning is likely to fail.", call. = FALSE)
  }
  minarea
}

.sc_slic_post = function(slic_sf, metadata, step, compactness, iter_attr) {
  if (!is.list(slic_sf) || inherits(slic_sf, "sf")) {
    slic_sf = list(slic_sf)
  }
  slic_sf = update_supercells_ids(slic_sf)
  if (isFALSE(metadata)) {
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
  slic_sf
}
