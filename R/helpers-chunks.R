# update supercells ids for chunked sf outputs
.sc_chunk_offset_ids_sf = function(x){
  x = x[lapply(x, length) > 0]
  max_id = 0
  for (i in seq_along(x)) {
    if (max_id > 0) {
      x[[i]][["supercells"]] = x[[i]][["supercells"]] + max_id
    }
    max_id = max(x[[i]][["supercells"]], na.rm = TRUE) + 1
  }
  x
}

# combine chunked sf outputs with unique ids
.sc_chunk_update_ids = function(x){
  x = .sc_chunk_offset_ids_sf(x)
  x = do.call(rbind, x)
  return(x)
}

# update supercells ids for raster chunks
.sc_chunk_offset_ids_raster = function(rasters){
  if (length(rasters) <= 1) {
    return(rasters)
  }
  max_id = 0
  for (i in seq_along(rasters)) {
    r = rasters[[i]]
    if (max_id > 0) {
      r = r + max_id
    }
    curr_max = terra::global(r, "max", na.rm = TRUE)[1, 1]
    if (!is.na(curr_max)) {
      max_id = curr_max + 1
    }
    rasters[[i]] = r
  }
  return(rasters)
}

# predict memory usage in gb (approximate, conservative)
# Accounts for: R values matrix, duplicated values vector in C++,
# and per-pixel distance/cluster buffers (assuming cleaning enabled).
.sc_chunk_mem_gb = function(dim_x){
  ncell = dim_x[1] * dim_x[2]
  bands = dim_x[3]
  bytes_per_cell = (16 * bands) + 16
  mem_bytes = ncell * bytes_per_cell
  mem_bytes / (1024 * 1024 * 1024)
}

# find an approximate optimal chunk size
.sc_chunk_optimize_size = function(dim_x, limit, step = NULL){
  max_dim = max(dim_x[1:2])
  bands = dim_x[3]
  bytes_per_cell = (16 * bands) + 16
  target_cells = floor((limit * 1024 * 1024 * 1024) / bytes_per_cell)
  wsize = floor(sqrt(target_cells))
  if (wsize < 1) {
    wsize = 1
  }
  if (wsize > max_dim) {
    wsize = max_dim
  }
  if (!is.null(step) && is.numeric(step) && length(step) == 1 && step > 0) {
    wsize = floor(wsize / step) * step
    if (wsize < step) {
      wsize = step
    }
    if (wsize > max_dim) {
      wsize = max_dim
    }
  }
  return(wsize)
}

# compute chunk extents for splitting
# if limit = FALSE, the extent of the whole input is returned
# if limit = TRUE, the extent of the input is split into chunks,
#                  where the size of each raster chunk is optimized to be as close to
#                  the (hardcoded) limit of 1GB as possible
# if limit is numeric, the extent of the input is split into chunks,
#                      where the width/height of each chunk is equal to the limit
.sc_chunk_extents = function(dim_x, limit, step = NULL){
  if (is.numeric(limit)){
    wsize = limit
    limit = 0
    dims1 = ceiling(.sc_util_seq_last(0, to = dim_x[1], by = wsize))
    dims2 = ceiling(.sc_util_seq_last(0, to = dim_x[2], by = wsize))
  } else if (!limit){
    limit = Inf
  } else {
    limit = getOption("supercells.chunk_mem_gb", 4)
    if (!is.numeric(limit) || length(limit) != 1 || is.na(limit) || limit <= 0) {
      stop("The 'supercells.chunk_mem_gb' option must be a positive number", call. = FALSE)
    }
    wsize = .sc_chunk_optimize_size(dim_x, limit, step = step)
    dims1 = ceiling(seq.int(0, to = dim_x[1],
                            length.out = as.integer((dim_x[1] - 1) / wsize + 1) + 1))
    dims2 = ceiling(seq.int(0, to = dim_x[2],
                            length.out = as.integer((dim_x[2] - 1) / wsize + 1) + 1))
  }
  if (.sc_chunk_mem_gb(dim_x) > limit){
    row_dims = seq_along(dims1)[-length(dims1)]
    col_dims = seq_along(dims2)[-length(dims2)]
    n_chunks = max(row_dims) * max(col_dims)
    row_cols_chunks = cbind(min_row = integer(length = n_chunks),
                            max_row = integer(length = n_chunks),
                            min_col = integer(length = n_chunks),
                            max_col = integer(length = n_chunks))
    l = 0
    for (i in row_dims){
      for (j in col_dims){
        l = l + 1
        row_cols_chunks[l, 1] = dims1[i] + 1
        row_cols_chunks[l, 2] = dims1[i + 1]
        row_cols_chunks[l, 3] = dims2[j] + 1
        row_cols_chunks[l, 4] = dims2[j + 1]
      }
    }
  } else {
    row_cols_chunks = cbind(min_row = 1,
                            max_row = dim_x[1],
                            min_col = 1,
                            max_col = dim_x[2])
  }
  return(row_cols_chunks)
}
