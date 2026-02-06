# Prototype: merge adjacent supercells based on value-space distance.
# This is a minimal, adjacency-constrained greedy merge using sf input.

# prepare a distance function wrapper
# input: dist_fun (string name or function)
# output: function(a, b) -> numeric distance
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

# build adjacency pairs and distances
# input: geoms (sfc), vals (matrix), dist_one (function)
# output: list(pairs = int matrix [m x 2], dists = numeric [m])
.sc_merge_pairs_adjacent = function(geoms, vals, dist_one) {
  adj = sf::st_touches(geoms)
  if (length(adj) == 0) {
    return(list(pairs = matrix(integer(0), ncol = 2), dists = numeric(0)))
  }
  pairs_list = lapply(seq_along(adj), function(i) {
    nb = adj[[i]]
    nb = nb[nb > i]
    if (length(nb) == 0) return(NULL)
    cbind(rep.int(i, length(nb)), nb)
  })
  pairs = do.call(rbind, pairs_list)
  if (is.null(pairs)) {
    return(list(pairs = matrix(integer(0), ncol = 2), dists = numeric(0)))
  }
  dists = vapply(seq_len(nrow(pairs)), function(k) {
    dist_one(vals[pairs[k, 1], ], vals[pairs[k, 2], ])
  }, numeric(1))
  list(pairs = pairs, dists = dists)
}

# aggregate values/weights by group label
# input: groups (int vector), vals (matrix), w (numeric)
# output: list(group = reindexed groups, vals, w)
.sc_merge_aggregate_components = function(groups, vals, w) {
  keep_ids = sort(unique(groups))
  id_map = setNames(seq_len(length(keep_ids)), keep_ids)
  group_ids = unname(id_map[as.character(groups)])
  new_vals = matrix(0, nrow = length(keep_ids), ncol = ncol(vals))
  new_w = numeric(length(keep_ids))
  for (idx in seq_along(groups)) {
    gid = group_ids[idx]
    wi = w[idx]
    new_vals[gid, ] = new_vals[gid, ] + wi * vals[idx, ]
    new_w[gid] = new_w[gid] + wi
  }
  new_vals = new_vals / new_w
  list(group = group_ids, vals = new_vals, w = new_w)
}

# finalize sf output from group labels and aggregated stats
# input: x (sf), geoms (sfc), group (int), vals_out (matrix), w_out (numeric)
#        value_cols (char), weight_is_col (bool), has_xy (bool), crs_x (crs),
#        dissolve_geoms (function)
# output: sf with merged geometry and updated attributes
.sc_merge_finalize_groups = function(x, geoms, group, vals_out, w_out,
                                     value_cols, weight_is_col, has_xy,
                                     crs_x, dissolve_geoms) {
  group_levels = seq_len(max(group))
  geom_idx = split(seq_len(nrow(x)), factor(group, levels = group_levels))
  keep = vapply(geom_idx, `[`, integer(1), 1L)
  out = x[keep, , drop = FALSE]
  out_geom = sf::st_sfc(lapply(geom_idx, function(idx) sf::st_union(geoms[idx])[[1]]), crs = crs_x)
  out$geometry = dissolve_geoms(out_geom)
  out[value_cols] = vals_out
  if (weight_is_col) {
    out[[weight]] = w_out
  }
  if (has_xy) {
    coords = sf::st_coordinates(sf::st_centroid(out$geometry))
    out$x = coords[, 1]
    out$y = coords[, 2]
  }
  out
}

# Felzenszwalb-Huttenlocher (FH) merge
# input: geoms (sfc), vals (matrix), w (numeric), dist_one (function), kappa (numeric)
# output: list(group, vals, w)
.sc_merge_fh = function(geoms, vals, w, dist_one, kappa) {
  n = nrow(vals)
  pair_data = .sc_merge_pairs_adjacent(geoms, vals, dist_one)
  pairs = pair_data$pairs
  dists = pair_data$dists
  if (length(pairs) == 0) {
    return(list(group = seq_len(n), vals = vals, w = w))
  }

  ord = order(dists)
  pairs = pairs[ord, , drop = FALSE]
  dists = dists[ord]

  parent = seq_len(n)
  find_root = function(i) {
    while (parent[i] != i) {
      i = parent[i]
    }
    i
  }
  comp_size = rep(1L, n)
  comp_int = rep(0.0, n)

  for (k in seq_len(nrow(pairs))) {
    i = pairs[k, 1]
    j = pairs[k, 2]
    ri = find_root(i)
    rj = find_root(j)
    if (ri == rj) next
    w_ij = dists[k]
    thresh_i = comp_int[ri] + kappa / comp_size[ri]
    thresh_j = comp_int[rj] + kappa / comp_size[rj]
    if (w_ij <= min(thresh_i, thresh_j)) {
      if (comp_size[ri] < comp_size[rj]) {
        tmp = ri; ri = rj; rj = tmp
      }
      parent[rj] = ri
      comp_size[ri] = comp_size[ri] + comp_size[rj]
      comp_int[ri] = max(comp_int[ri], comp_int[rj], w_ij)
    }
  }
  groups = vapply(seq_len(n), find_root, integer(1))
  .sc_merge_aggregate_components(groups, vals, w)
}

# MST-based merge
# input: geoms (sfc), vals (matrix), w (numeric), dist_one (function),
#        target_k (int or NULL), tau (numeric or NULL)
# output: list(group, vals, w)
.sc_merge_mst = function(geoms, vals, w, dist_one, target_k, tau) {
  n = nrow(vals)
  pair_data = .sc_merge_pairs_adjacent(geoms, vals, dist_one)
  pairs = pair_data$pairs
  dists = pair_data$dists
  if (length(pairs) == 0) {
    return(list(group = seq_len(n), vals = vals, w = w))
  }

  ord = order(dists)
  pairs = pairs[ord, , drop = FALSE]
  dists = dists[ord]

  parent = seq_len(n)
  find_root = function(i) {
    while (parent[i] != i) {
      i = parent[i]
    }
    i
  }
  mst_edges = list()
  mst_weights = numeric(0)
  for (k in seq_len(nrow(pairs))) {
    i = pairs[k, 1]
    j = pairs[k, 2]
    ri = find_root(i)
    rj = find_root(j)
    if (ri != rj) {
      parent[rj] = ri
      mst_edges[[length(mst_edges) + 1L]] = c(i, j)
      mst_weights[length(mst_weights) + 1L] = dists[k]
      if (length(mst_edges) == n - 1) break
    }
  }
  if (length(mst_edges) == 0) {
    return(list(group = seq_len(n), vals = vals, w = w))
  }

  keep_edges = rep(TRUE, length(mst_edges))
  if (!is.null(target_k)) {
    cuts = target_k - 1
    if (cuts > 0) {
      cut_idx = order(mst_weights, decreasing = TRUE)[seq_len(cuts)]
      keep_edges[cut_idx] = FALSE
    }
  } else if (!is.null(tau)) {
    keep_edges = mst_weights <= tau
  }

  parent = seq_len(n)
  find_root = function(i) {
    while (parent[i] != i) {
      i = parent[i]
    }
    i
  }
  for (k in seq_along(mst_edges)) {
    if (!keep_edges[k]) next
    i = mst_edges[[k]][1]
    j = mst_edges[[k]][2]
    ri = find_root(i)
    rj = find_root(j)
    if (ri != rj) {
      parent[rj] = ri
    }
  }
  groups = vapply(seq_len(n), find_root, integer(1))
  .sc_merge_aggregate_components(groups, vals, w)
}

# greedy adjacency-constrained merge
# input: geoms (sfc), vals (matrix), w (numeric), target_k/tau,
#        dist_one (function), verbose (bool)
# output: list(group, vals, w)
.sc_merge_greedy = function(geoms, vals, w, target_k, tau, dist_one, verbose) {
  n0 = nrow(vals)
  vals0 = vals
  w0 = w
  alive = rep(TRUE, n0)
  n_alive = n0
  neighbors = sf::st_touches(geoms)
  parent = seq_len(n0)
  find_root = function(i) {
    while (parent[i] != i) {
      i = parent[i]
    }
    i
  }

  repeat {
    if (!is.null(target_k) && n_alive <= target_k) break
    min_dist = Inf
    i = NA_integer_
    j = NA_integer_
    for (idx in which(alive)) {
      nb = neighbors[[idx]]
      if (length(nb) == 0) next
      for (k in nb) {
        if (!alive[k] || k <= idx) next
        d = dist_one(vals[idx, ], vals[k, ])
        if (!is.finite(d)) next
        if (d < min_dist) {
          min_dist = d
          i = idx
          j = k
        }
      }
    }
    if (!is.finite(min_dist)) break
    if (!is.null(tau) && min_dist > tau) break

    if (verbose) {
      message(sprintf("Merging %d and %d (dist=%.4f)", i, j, min_dist))
    }

    w_i = w[i]
    w_j = w[j]
    w_new = w_i + w_j
    vals[i, ] = (w_i * vals[i, ] + w_j * vals[j, ]) / w_new
    w[i] = w_new

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
    parent[j] = i
  }

  groups = vapply(seq_len(n0), find_root, integer(1))
  .sc_merge_aggregate_components(groups, vals0, w0)
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

  get_opt = function(name) {
    if (is.list(method_opts) && name %in% names(method_opts)) {
      return(method_opts[[name]])
    }
    NULL
  }
  target_k = get_opt("target_k")
  tau = get_opt("tau")
  if (is.null(target_k) && is.null(tau) && method != "fh") {
    stop("Provide target_k or tau to control merging", call. = FALSE)
  }
  if (!is.null(target_k) && (!is.numeric(target_k) || length(target_k) != 1 || target_k < 1)) {
    stop("target_k must be a single positive number", call. = FALSE)
  }
  if (!is.null(target_k) && target_k > nrow(x)) {
    stop("target_k cannot exceed the number of supercells", call. = FALSE)
  }
  if (method == "fh") {
    kappa = get_opt("kappa")
    if (is.null(kappa) || !is.numeric(kappa) || length(kappa) != 1) {
      stop("Provide a single numeric kappa for method = 'fh'", call. = FALSE)
    }
  }

  crs_x = sf::st_crs(x)
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

  dist_one = .sc_merge_dist_fun(dist_fun)
  geoms = sf::st_geometry(x)
  finalize_groups = function(res) {
    .sc_merge_finalize_groups(x, geoms, res$group, res$vals, res$w,
                              value_cols, weight_is_col, has_xy,
                              crs_x, dissolve_geoms)
  }

  if (method == "greedy") {
    res = .sc_merge_greedy(geoms, vals, w, target_k, tau, dist_one, verbose)
  } else if (method == "fh") {
    res = .sc_merge_fh(geoms, vals, w, dist_one, kappa)
  } else {
    res = .sc_merge_mst(geoms, vals, w, dist_one, target_k, tau)
  }
  finalize_groups(res)

#
# --- Archived FH/MST implementations (commented) ---
# Removed from inline comments now that FH/MST are implemented as helpers above.
}
