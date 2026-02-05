# Prototype: merge adjacent supercells based on value-space distance.
# This is a minimal, adjacency-constrained greedy merge using sf input.

.sc_merge_dist_fun = function(dist_fun) {
  if (is.function(dist_fun)) {
    return(function(a, b) sc_dist_vec_cpp(a, b, "", dist_fun))
  }
  if (!is.character(dist_fun) || length(dist_fun) != 1 || is.na(dist_fun)) {
    stop("dist_fun must be a function or a single string", call. = FALSE)
  }
  if (!(dist_fun %in% c("euclidean", "jsd", "dtw", "dtw2d", philentropy::getDistMethods()))) {
    stop("Unsupported dist_fun; provide a function for this distance", call. = FALSE)
  }
  dummy_fun = function() ""
  return(function(a, b) sc_dist_vec_cpp(a, b, dist_fun, dummy_fun))
}

sc_merge_supercells = function(x, dist_fun = "euclidean",
                               method = c("greedy", "fh", "mst"),
                               method_opts = list(),
                               weight = "area", verbose = FALSE) {
  if (!inherits(x, "sf")) {
    stop("x must be an sf object (output of sc_slic)", call. = FALSE)
  }
  if (isTRUE(sf::st_is_longlat(x))) {
    stop("sc_merge_supercells requires projected coordinates; reproject x before merging", call. = FALSE)
  }
  old_s2 = sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  method = match.arg(method)
  if (method != "greedy") {
    stop("Only method = 'greedy' is enabled for now.", call. = FALSE)
  }

  if (nrow(x) < 2) {
    return(x)
  }

  x_df = sf::st_drop_geometry(x)

  skip_cols = c("supercells", "x", "y")
  if (is.character(weight) && weight %in% names(x_df)) {
    skip_cols = c(skip_cols, weight)
  }
  value_cols = setdiff(names(x_df), skip_cols)
  value_cols = value_cols[sapply(x_df[, value_cols, drop = FALSE], is.numeric)]
  if (length(value_cols) == 0) {
    stop("No numeric value columns found for merging", call. = FALSE)
  }
  value_mat = x_df[, value_cols, drop = FALSE]
  na_cols = colSums(is.na(value_mat)) > 0
  if (any(na_cols)) {
    value_cols = value_cols[!na_cols]
    if (length(value_cols) == 0) {
      stop("All numeric value columns contain NA; cannot compute distances", call. = FALSE)
    }
    if (verbose) {
      message("Dropping value columns with NA: ", paste(names(value_mat)[na_cols], collapse = ", "))
    }
  }

  weight_is_col = is.character(weight) && length(weight) == 1 && weight %in% names(x_df)
  if (is.character(weight) && length(weight) == 1) {
    if (weight == "area") {
      w = as.numeric(sf::st_area(x))
    } else if (weight_is_col) {
      w = as.numeric(x_df[[weight]])
    } else {
      stop("weight must be 'area', a column name, or numeric vector", call. = FALSE)
    }
  } else if (is.numeric(weight) && length(weight) == nrow(x)) {
    w = as.numeric(weight)
  } else {
    stop("weight must be 'area', a column name, or numeric vector", call. = FALSE)
  }

  vals = as.matrix(x_df[, value_cols, drop = FALSE])
  has_xy = all(c("x", "y") %in% names(x_df))

  dist_one = .sc_merge_dist_fun(dist_fun)

  get_opt = function(name) {
    if (is.list(method_opts) && name %in% names(method_opts)) {
      return(method_opts[[name]])
    }
    NULL
  }
  target_k = get_opt("target_k")
  tau = get_opt("tau")
  if (is.null(target_k) && is.null(tau)) {
    stop("Provide target_k or tau to control merging", call. = FALSE)
  }
  if (!is.null(target_k) && (!is.numeric(target_k) || length(target_k) != 1 || target_k < 1)) {
    stop("target_k must be a single positive number", call. = FALSE)
  }
  if (!is.null(target_k) && target_k > nrow(x)) {
    stop("target_k cannot exceed the number of supercells", call. = FALSE)
  }

  crs_x = sf::st_crs(x)
  merge_geoms = function(g1, g2) {
    merged = sf::st_combine(sf::st_sfc(g1, g2, crs = crs_x))
    merged[[1]]
  }

  dissolve_geoms = function(geoms) {
    out = vector("list", length(geoms))
    for (k in seq_along(geoms)) {
      g = sf::st_sfc(geoms[[k]], crs = crs_x)
      g_type = sf::st_geometry_type(g)
      if (any(g_type %in% c("MULTIPOLYGON", "GEOMETRYCOLLECTION"))) {
        g = suppressWarnings(sf::st_cast(g, "POLYGON"))
        g = sf::st_union(g)
      }
      if (length(g) == 0) g = sf::st_sfc(geoms[[k]], crs = crs_x)
      out[[k]] = g[[1]]
    }
    sf::st_sfc(out, crs = crs_x)
  }

  n0 = nrow(x)
  alive = rep(TRUE, n0)
  n_alive = n0
  geoms = sf::st_geometry(x)
  neighbors = sf::st_touches(geoms)
  version = integer(n0)
  heap_d = numeric(0)
  heap_i = integer(0)
  heap_j = integer(0)
  heap_vi = integer(0)
  heap_vj = integer(0)
  heap_n = 0L
  heap_push = function(d, i, j, vi, vj) {
    if (!is.finite(d)) return(invisible(NULL))
    heap_n <<- heap_n + 1L
    heap_d[heap_n] <<- d
    heap_i[heap_n] <<- i
    heap_j[heap_n] <<- j
    heap_vi[heap_n] <<- vi
    heap_vj[heap_n] <<- vj
    k = heap_n
    while (k > 1L) {
      p = k %/% 2L
      if (heap_d[p] <= heap_d[k]) break
      tmp = heap_d[p]; heap_d[p] <<- heap_d[k]; heap_d[k] <<- tmp
      tmp = heap_i[p]; heap_i[p] <<- heap_i[k]; heap_i[k] <<- tmp
      tmp = heap_j[p]; heap_j[p] <<- heap_j[k]; heap_j[k] <<- tmp
      tmp = heap_vi[p]; heap_vi[p] <<- heap_vi[k]; heap_vi[k] <<- tmp
      tmp = heap_vj[p]; heap_vj[p] <<- heap_vj[k]; heap_vj[k] <<- tmp
      k = p
    }
    invisible(NULL)
  }
  heap_pop_min = function() {
    if (heap_n == 0L) return(NULL)
    out = list(d = heap_d[1L], i = heap_i[1L], j = heap_j[1L],
               vi = heap_vi[1L], vj = heap_vj[1L])
    if (heap_n == 1L) {
      heap_n <<- 0L
      return(out)
    }
    heap_d[1L] <<- heap_d[heap_n]
    heap_i[1L] <<- heap_i[heap_n]
    heap_j[1L] <<- heap_j[heap_n]
    heap_vi[1L] <<- heap_vi[heap_n]
    heap_vj[1L] <<- heap_vj[heap_n]
    heap_n <<- heap_n - 1L
    k = 1L
    repeat {
      l = k * 2L
      r = l + 1L
      if (l > heap_n) break
      m = l
      if (r <= heap_n && heap_d[r] < heap_d[l]) m = r
      if (heap_d[k] <= heap_d[m]) break
      tmp = heap_d[k]; heap_d[k] <<- heap_d[m]; heap_d[m] <<- tmp
      tmp = heap_i[k]; heap_i[k] <<- heap_i[m]; heap_i[m] <<- tmp
      tmp = heap_j[k]; heap_j[k] <<- heap_j[m]; heap_j[m] <<- tmp
      tmp = heap_vi[k]; heap_vi[k] <<- heap_vi[m]; heap_vi[m] <<- tmp
      tmp = heap_vj[k]; heap_vj[k] <<- heap_vj[m]; heap_vj[m] <<- tmp
      k = m
    }
    out
  }
  for (i in seq_len(n0)) {
    nb = neighbors[[i]]
    nb = nb[nb > i]
    if (length(nb) == 0) next
    d = vapply(nb, function(j) dist_one(vals[i, ], vals[j, ]), numeric(1))
    for (k in seq_along(nb)) {
      heap_push(d[k], i, nb[k], version[i], version[nb[k]])
    }
  }

  repeat {
    if (!is.null(target_k) && n_alive <= target_k) break
    item = NULL
    repeat {
      item = heap_pop_min()
      if (is.null(item)) break
      i = item$i
      j = item$j
      if (!alive[i] || !alive[j]) next
      if (item$vi != version[i] || item$vj != version[j]) next
      if (!(j %in% neighbors[[i]])) next
      break
    }
    if (is.null(item)) break
    min_dist = item$d
    if (!is.null(tau) && min_dist > tau) break

    if (verbose) {
      message(sprintf("Merging %d and %d (dist=%.4f)", i, j, min_dist))
    }

    w_i = w[i]
    w_j = w[j]
    w_new = w_i + w_j
    vals[i, ] = (w_i * vals[i, ] + w_j * vals[j, ]) / w_new
    w[i] = w_new

    geoms[[i]] = merge_geoms(geoms[[i]], geoms[[j]])

    old_i = neighbors[[i]]
    old_j = neighbors[[j]]
    new_neighbors = union(old_i, old_j)
    new_neighbors = setdiff(new_neighbors, c(i, j))
    new_neighbors = new_neighbors[alive[new_neighbors]]

    for (k in old_j) {
      if (!alive[k] || k == i) next
      nk = neighbors[[k]]
      nk = nk[nk != j]
      if (!(i %in% nk)) nk = c(nk, i)
      neighbors[[k]] = nk
    }
    for (k in old_i) {
      if (!alive[k] || k == j) next
      nk = neighbors[[k]]
      nk = nk[nk != j]
      neighbors[[k]] = nk
    }
    neighbors[[i]] = new_neighbors
    neighbors[[j]] = integer(0)

    alive[j] = FALSE
    n_alive = n_alive - 1L
    version[i] = version[i] + 1L

    for (k in new_neighbors) {
      a = min(i, k)
      b = max(i, k)
      d = dist_one(vals[a, ], vals[b, ])
      heap_push(d, a, b, version[a], version[b])
    }
  }

  keep = which(alive)
  out = x[keep, , drop = FALSE]
  out$geometry = dissolve_geoms(geoms[keep])
  out[value_cols] = vals[keep, , drop = FALSE]
  if (has_xy) {
    coords = sf::st_coordinates(sf::st_centroid(out$geometry))
    out$x = coords[, 1]
    out$y = coords[, 2]
  }
  if (weight_is_col) {
    out[[weight]] = w[keep]
  }
  out

# --- Archived FH/MST implementations (commented) ---
# aggregate_components = function(groups) {
#   keep_ids = sort(unique(groups))
#   new_n = length(keep_ids)
#   new_vals = matrix(0, nrow = new_n, ncol = ncol(vals))
#   new_w = numeric(new_n)
#   new_geom = vector("list", new_n)
#   new_centers = if (has_xy) matrix(0, nrow = new_n, ncol = 2) else NULL
#
#   id_map = setNames(seq_len(new_n), keep_ids)
#   for (idx in seq_len(nrow(x))) {
#     gid = id_map[as.character(groups[idx])]
#     wi = w[idx]
#     new_vals[gid, ] = new_vals[gid, ] + wi * vals[idx, ]
#     new_w[gid] = new_w[gid] + wi
#     if (is.null(new_geom[[gid]])) {
#       new_geom[[gid]] = sf::st_geometry(x[idx, ])[[1]]
#     } else {
#       new_geom[[gid]] = merge_geoms(new_geom[[gid]], sf::st_geometry(x[idx, ])[[1]])
#     }
#     if (has_xy) {
#       new_centers[gid, ] = new_centers[gid, ] + wi * centers[idx, ]
#     }
#   }
#
#   new_vals = new_vals / new_w
#   if (has_xy) {
#     new_centers = new_centers / new_w
#   }
#
#   out = x[keep_ids, , drop = FALSE]
#   out$geometry = sf::st_sfc(new_geom, crs = sf::st_crs(x))
#   out[value_cols] = new_vals
#   if (has_xy) {
#     out$x = new_centers[, 1]
#     out$y = new_centers[, 2]
#   }
#   out
# }
#
# new_union_find = function(n) {
#   parent = seq_len(n)
#   find_root = function(i) {
#     while (parent[i] != i) {
#       parent[i] <<- parent[parent[i]]
#       i = parent[i]
#     }
#     i
#   }
#   list(
#     parent = parent,
#     find_root = find_root,
#     union = function(a, b) {
#       ra = find_root(a)
#       rb = find_root(b)
#       if (ra == rb) return(FALSE)
#       parent[rb] <<- ra
#       TRUE
#     }
#   )
# }
#
# if (method == "fh") {
#   kappa = opts$kappa
#   n = nrow(x)
#   pair_data = build_pairs()
#   pairs = pair_data$pairs
#   dists = pair_data$dists
#   if (length(pairs) == 0) return(x)
#
#   ord = order(dists)
#   pairs = pairs[ord]
#   dists = dists[ord]
#
#   parent = seq_len(n)
#   find_root = function(i) {
#     while (parent[i] != i) {
#       parent[i] <<- parent[parent[i]]
#       i = parent[i]
#     }
#     i
#   }
#   comp_size = rep(1L, n)
#   comp_int = rep(0.0, n)
#   union = function(a, b, w) {
#     ra = find_root(a)
#     rb = find_root(b)
#     if (ra == rb) return(invisible(NULL))
#     if (comp_size[ra] < comp_size[rb]) {
#       tmp = ra; ra = rb; rb = tmp
#     }
#     parent[rb] <<- ra
#     comp_size[ra] <<- comp_size[ra] + comp_size[rb]
#     comp_int[ra] <<- max(comp_int[ra], comp_int[rb], w)
#     invisible(NULL)
#   }
#
#   for (k in seq_along(pairs)) {
#     i = pairs[[k]][1]
#     j = pairs[[k]][2]
#     ri = find_root(i)
#     rj = find_root(j)
#     if (ri == rj) next
#     w_ij = dists[k]
#     thresh_i = comp_int[ri] + kappa / comp_size[ri]
#     thresh_j = comp_int[rj] + kappa / comp_size[rj]
#     if (w_ij <= min(thresh_i, thresh_j)) {
#       union(ri, rj, w_ij)
#     }
#   }
#   groups = vapply(seq_len(n), find_root, integer(1))
#   return(aggregate_components(groups))
# }
#
# if (method == "mst") {
#   target_k = opts$target_k
#   tau = opts$tau
#   n = nrow(x)
#   pair_data = build_pairs()
#   pairs = pair_data$pairs
#   dists = pair_data$dists
#   if (length(pairs) == 0) return(x)
#
#   ord = order(dists)
#   pairs = pairs[ord]
#   dists = dists[ord]
#
#   mst_edges = list()
#   mst_weights = numeric(0)
#   uf_mst = new_union_find(n)
#   for (k in seq_along(pairs)) {
#     i = pairs[[k]][1]
#     j = pairs[[k]][2]
#     if (uf_mst$union(i, j)) {
#       mst_edges[[length(mst_edges) + 1]] = c(i, j)
#       mst_weights[length(mst_weights) + 1] = dists[k]
#       if (length(mst_edges) == n - 1) break
#     }
#   }
#   if (length(mst_edges) == 0) return(x)
#
#   keep_edges = rep(TRUE, length(mst_edges))
#   if (!is.null(target_k)) {
#     cuts = target_k - 1
#     if (cuts > 0) {
#       cut_idx = order(mst_weights, decreasing = TRUE)[seq_len(cuts)]
#       keep_edges[cut_idx] = FALSE
#     }
#   } else if (!is.null(tau)) {
#     keep_edges = mst_weights <= tau
#   }
#
#   uf = new_union_find(n)
#   for (k in seq_along(mst_edges)) {
#     if (!keep_edges[k]) next
#     i = mst_edges[[k]][1]
#     j = mst_edges[[k]][2]
#     uf$union(i, j)
#   }
#   groups = vapply(seq_len(n), uf$find_root, integer(1))
#   return(aggregate_components(groups))
# }
}
