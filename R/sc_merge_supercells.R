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
  method = match.arg(method)
  if (method != "greedy") {
    stop("Only method = 'greedy' is enabled for now.",
         call. = FALSE)
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

  if (is.character(weight) && length(weight) == 1) {
    if (weight == "area") {
      w = as.numeric(sf::st_area(x))
    } else if (weight %in% names(x_df)) {
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
  centers = if (has_xy) as.matrix(x_df[, c("x", "y"), drop = FALSE]) else NULL

  dist_one = .sc_merge_dist_fun(dist_fun)

  get_opt = function(name, default = NULL) {
    if (is.list(method_opts) && name %in% names(method_opts)) {
      return(method_opts[[name]])
    }
    default
  }

  opts = list(
    target_k = get_opt("target_k", NULL),
    tau = get_opt("tau", NULL),
    kappa = get_opt("kappa", 0.5)
  )

  validate_controls = function() {
    target_k = opts$target_k
    tau = opts$tau
    if (is.null(target_k) && is.null(tau)) {
      stop("Provide target_k or tau to control merging", call. = FALSE)
    }
    if (!is.null(target_k) && (!is.numeric(target_k) || length(target_k) != 1 || target_k < 1)) {
      stop("target_k must be a single positive number", call. = FALSE)
    }
    if (!is.null(target_k) && target_k > nrow(x)) {
      stop("target_k cannot exceed the number of supercells", call. = FALSE)
    }
    invisible(TRUE)
  }

  build_pairs = function() {
    n = nrow(x)
    adj = sf::st_touches(x)
    pairs = vector("list", length = n)
    pair_count = 0
    for (i in seq_len(n)) {
      if (length(adj[[i]]) == 0) next
      for (j in adj[[i]]) {
        if (j > i) {
          pair_count = pair_count + 1
          pairs[[pair_count]] = c(i, j)
        }
      }
    }
    if (pair_count == 0) return(list(pairs = list(), dists = numeric(0)))
    pairs = pairs[seq_len(pair_count)]
    dists = vapply(pairs, function(idx) {
      dist_one(vals[idx[1], ], vals[idx[2], ])
    }, numeric(1))
    list(pairs = pairs, dists = dists)
  }

  merge_geoms = function(g1, g2) {
    merged = sf::st_union(sf::st_sfc(g1, crs = sf::st_crs(x)), sf::st_sfc(g2, crs = sf::st_crs(x)))
    if (any(sf::st_geometry_type(merged) %in% c("GEOMETRYCOLLECTION", "MULTIPOLYGON"))) {
      merged = suppressWarnings(sf::st_collection_extract(merged, "POLYGON"))
    }
    if (length(merged) == 0) {
      merged = sf::st_union(sf::st_sfc(g1, crs = sf::st_crs(x)), sf::st_sfc(g2, crs = sf::st_crs(x)))
    }
    merged[[1]]
  }

  validate_controls()

  repeat {
    target_k = opts$target_k
    tau = opts$tau
    n = nrow(x)
    if (!is.null(target_k) && n <= target_k) break

    pair_data = build_pairs()
    pairs = pair_data$pairs
    dists = pair_data$dists
    if (length(pairs) == 0) break

    min_idx = which.min(dists)
    min_dist = dists[min_idx]
    if (!is.null(tau) && min_dist > tau) break

    i = pairs[[min_idx]][1]
    j = pairs[[min_idx]][2]

    if (verbose) {
      message(sprintf("Merging %d and %d (dist=%.4f)", i, j, min_dist))
    }

    w_i = w[i]
    w_j = w[j]
    w_new = w_i + w_j
    vals[i, ] = (w_i * vals[i, ] + w_j * vals[j, ]) / w_new
    w[i] = w_new

    if (has_xy) {
      centers[i, ] = (w_i * centers[i, ] + w_j * centers[j, ]) / w_new
    }

    new_geom = merge_geoms(sf::st_geometry(x[i, ])[[1]], sf::st_geometry(x[j, ])[[1]])
    x$geometry[i] = sf::st_sfc(new_geom, crs = sf::st_crs(x))[[1]]
    x = x[-j, , drop = FALSE]
    w = w[-j]
    vals = vals[-j, , drop = FALSE]
    if (has_xy) {
      centers = centers[-j, , drop = FALSE]
    }

    x[value_cols] = vals
    if (has_xy) {
      x$x = centers[, 1]
      x$y = centers[, 2]
  }
}

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

  x
}
