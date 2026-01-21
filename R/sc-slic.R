#' Create supercells
#'
#' Creates supercells from single- or multi-band rasters using a SLIC-style algorithm.
#' The function supports either a target number of supercells (`k`) or a fixed grid
#' spacing (`step`), as well as optional custom centers and chunked processing.
#'
#' @param x An object of class SpatRaster (terra) or class stars (stars).
#' @param k The number of supercells desired (alternative to `step`).
#' @param step The distance (number of cells) between initial centers (alternative to `k`).
#' @param centers Optional sf object of custom centers. Requires `step`.
#' @param compactness A compactness value.
#' @param dist_fun A distance function name or a custom function. Supported names:
#' "euclidean", "jsd", "dtw", "dtw2d", or any method from `philentropy::getDistMethods()`.
#' A custom function must accept two numeric vectors and return a single numeric value.
#' @param avg_fun An averaging function name or custom function used to summarize
#' values within each supercell. Supported names: "mean" and "median". A custom
#' function must accept a numeric vector and return a single numeric value.
#' @param clean Should connectivity of the supercells be enforced?
#' @param iter Number of iterations.
#' @param transform Optional transformation applied before segmentation. Currently
#' supports "to_LAB" for RGB inputs.
#' @param minarea Minimal size of a supercell (in cells).
#' @param metadata Logical. Should metadata columns be kept?
#' @param chunks Chunking option. Use `FALSE` for no chunking, `TRUE` for
#' automatic chunking based on size, or a numeric value for a fixed chunk size
#' (in number of cells per side).
#' @param future Logical. Use future for parallelization?
#' @param verbose Verbosity level.
#' @param iter_diagnostics Logical. If `TRUE`, returns iteration diagnostics as an attribute
#' (`iter_diagnostics`) on the output. Only available when chunks are not used.
#'
#' @return An sf object with the supercell polygons and summary statistics.
#' Attributes `step` and `compactness` are attached to the result.
#' @export
sc_slic = function(x, k = NULL, step = NULL, centers = NULL, compactness, dist_fun = "euclidean",
                     avg_fun = "mean", clean = TRUE, iter = 10, transform = NULL,
                     minarea, metadata = FALSE, chunks = FALSE, future = FALSE, verbose = 0,
                     iter_diagnostics = FALSE) {
  prep = .sc_slic_prep(x, k, step, centers, compactness, dist_fun, avg_fun,
                       minarea, chunks, iter_diagnostics)

  if (iter == 0) {
    clean = FALSE
  }

  slic_sf = .sc_create_run(prep$x, prep$chunk_ext, prep$step, compactness,
                           prep$funs$dist_name, prep$funs$dist_fun,
                           prep$funs$avg_fun_name, prep$funs$avg_fun_fun,
                           clean, iter, prep$minarea, transform, prep$input_centers,
                           verbose, future, prep$iter_diagnostics)

  iter_attr = NULL
  if (prep$iter_diagnostics && is.list(slic_sf) && length(slic_sf) > 0) {
    iter_attr = attr(slic_sf[[1]], "iter_diagnostics")
  }
  .sc_create_post(slic_sf, metadata, prep$step, compactness, iter_attr)
}

.sc_slic_prep = function(x, k, step, centers, compactness, dist_fun, avg_fun,
                         minarea, chunks, iter_diagnostics) {
  x = .sc_create_input(x)
  .sc_create_validate(k, step, centers, compactness, chunks)
  step = .sc_create_step(x, k, step)
  input_centers = .sc_create_centers(x, centers)
  funs = .sc_create_funs(dist_fun, avg_fun)
  minarea = .sc_create_minarea(minarea, step)
  chunk_ext = prep_chunks_ext(dim(x), limit = chunks)
  if (iter_diagnostics && nrow(chunk_ext) > 1) {
    warning("Iteration diagnostics are only available when chunks = FALSE (single chunk). Iteration diagnostics were disabled.", call. = FALSE)
    iter_diagnostics = FALSE
  }
  list(x = x, step = step, input_centers = input_centers, funs = funs,
       minarea = minarea, chunk_ext = chunk_ext, iter_diagnostics = iter_diagnostics)
}

.sc_create_input = function(x) {
  if (inherits(x, "SpatRaster")) {
    return(x)
  }
  if (inherits(x, "stars")) {
    return(terra::rast(x))
  }
  stop("The SpatRaster class is expected as an input", call. = FALSE)
}

.sc_create_validate = function(k, step, centers, compactness, chunks) {
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

.sc_create_step = function(x, k, step) {
  if (!is.null(step)) {
    return(step)
  }
  mat = dim(x)[1:2]
  superpixelsize = round((mat[1] * mat[2]) / k + 0.5)
  round(sqrt(superpixelsize) + 0.5)
}

.sc_create_centers = function(x, centers) {
  if (is.null(centers)) {
    return(matrix(c(0L, 0L), ncol = 2))
  }
  centers_to_dims(x, centers)
}

.sc_create_funs = function(dist_fun, avg_fun) {
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

.sc_create_minarea = function(minarea, step) {
  if (missing(minarea)) {
    return(0)
  }
  if (minarea > step^2) {
    warning("The provided minarea value is larger than than the average supercell (step^2). The connectivity cleaning is likely to fail.", call. = FALSE)
  }
  minarea
}

.sc_create_run = function(x, chunk_ext, step, compactness, dist_name, dist_fun,
                          avg_fun_name, avg_fun_fun, clean, iter, minarea,
                          transform, input_centers, verbose, future, iter_diagnostics) {
  if (future) {
    if (in_memory(x)) {
      names_x = names(x)
      x = terra::writeRaster(x, tempfile(fileext = ".tif"))
      names(x) = names_x
    }
    if (!in_memory(x)) {
      x = terra::sources(x)[[1]]
    }
    oopts = options(future.globals.maxSize = +Inf)
    on.exit(options(oopts))
    future.apply::future_apply(chunk_ext, MARGIN = 1, run_slic_chunks, x = x,
                               step = step, compactness = compactness, dist_name = dist_name,
                               dist_fun = dist_fun, avg_fun_fun = avg_fun_fun, avg_fun_name = avg_fun_name,
                               clean = clean, iter = iter, minarea = minarea, transform = transform,
                               input_centers = input_centers, verbose = verbose,
                               iter_diagnostics = iter_diagnostics,
                               future.seed = TRUE)
  } else {
    apply(chunk_ext, MARGIN = 1, run_slic_chunks, x = x,
          step = step, compactness = compactness, dist_name = dist_name,
          dist_fun = dist_fun, avg_fun_fun = avg_fun_fun, avg_fun_name = avg_fun_name,
          clean = clean, iter = iter, minarea = minarea, transform = transform,
          input_centers = input_centers, verbose = verbose,
          iter_diagnostics = iter_diagnostics)
  }
}

.sc_create_post = function(slic_sf, metadata, step, compactness, iter_attr) {
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
  return(slic_sf)
}
