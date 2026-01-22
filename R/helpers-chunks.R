# updates supercells ids for chunks
offset_supercells_ids_sf = function(x){
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

update_supercells_ids = function(x){
  x = offset_supercells_ids_sf(x)
  x = do.call(rbind, x)
  return(x)
}

# update supercells ids for raster chunks (not needed for a single raster)
offset_supercells_ids_raster = function(rasters){
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

# predicts (rough estimation) memory usage of the algorithm
pred_mem_usage = function(dim_x){
  mem_bytes = dim_x[1] * dim_x[2] * dim_x[3] * 8 #in bytes
  mem_gb = mem_bytes / (1024 * 1024 * 1024)
  mem_gb
}

# looks for the optimal chunk size
optimize_chunk_size = function(dim_x, limit, by = 500){
  min_diff_memory = function(a, dim_x, limit){
    abs((dim_x[3] * a^2 * 8 / (1024 * 1024 * 1024)) - limit)
  }
  opti = stats::optimize(min_diff_memory,
                  interval = c(seq(100, max(dim_x[1:2]), by = by), max(dim_x[1:2])),
                  dim_x, limit)
  return(opti$minimum)
}

# prepares the extents of chunks:
# if limit = FALSE, the extent of the whole input is returned
# if limit = TRUE, the extent of the input is split into chunks,
#                  where the size of each raster chunk is optimized to be as close to
#                  the (hardcoded) limit of 1GB as possible
# if limit is numeric, the extent of the input is split into chunks,
#                      where the width/height of each chunk is equal to the limit
prep_chunks_ext = function(dim_x, limit){
  if (is.numeric(limit)){
    wsize = limit
    limit = 0
    dims1 = ceiling(seq_last(0, to = dim_x[1], by = wsize))
    dims2 = ceiling(seq_last(0, to = dim_x[2], by = wsize))
  } else if (!limit){
    limit = Inf
  } else {
    limit = 1 #hardcoded limit
    wsize = optimize_chunk_size(dim_x, limit, by = 500)
    dims1 = ceiling(seq.int(0, to = dim_x[1],
                            length.out = as.integer((dim_x[1] - 1) / wsize + 1) + 1))
    dims2 = ceiling(seq.int(0, to = dim_x[2],
                            length.out = as.integer((dim_x[2] - 1) / wsize + 1) + 1))
  }
  if (pred_mem_usage(dim_x) > limit){
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
