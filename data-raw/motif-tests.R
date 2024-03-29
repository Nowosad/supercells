devtools::load_all()
library(motif)
library(stars)
library(terra)
library(tmap)

# 1 -----------------------------------------------------------------------
nlcd = read_stars(system.file("raster/nlcd2011.tif", package = "spDataLarge"))
coma_output = lsp_signature(nlcd, type = "cove", window = 10, normalization = "pdf", ordered = FALSE)

slic = supermotifs(coma_output, k = 300, compactness = 0.01, dist_fun = "jensen_shannon")

slic_cent = st_as_sf(st_drop_geometry(slic[2:3]), coords = c("x", "y"), crs = st_crs(slic))

tm_shape(nlcd) +
  tm_raster(legend.show = FALSE, drop.levels = TRUE) +
  tm_shape(slic) +
  tm_borders(col = "red") +
  tm_shape(slic_cent) +
  tm_dots()


# 2 -----------------------------------------------------------------------
landcover = read_stars(system.file("raster/landcover2015.tif", package = "motif"))
coma_output2 = lsp_signature(landcover, type = "cove", window = 10, normalization = "pdf", ordered = FALSE)

slic2 = supermotifs(coma_output2, k = 5000, compactness = 0.1, dist_fun = "jensen_shannon")

slic_cent2 = st_as_sf(st_drop_geometry(slic2[2:3]), coords = c("x", "y"), crs = st_crs(slic2))

tm_shape(landcover) +
  tm_raster(legend.show = FALSE, drop.levels = TRUE) +
  tm_shape(slic2) +
  tm_borders(col = "red") #+
  #tm_shape(slic_cent2) +
  #tm_dots()

write_stars(landcover, "landcover.tf")
write_sf(slic2, "slic2.gpkg")
# old (2021-05-06) --------------------------------------------------------
lsp_restructure = function(x){
  x_attr = attributes(x)
  nc = ncol(x$signature[[1]])

  unnested_signature = matrix(unlist(x$signature, use.names = FALSE),
                              ncol = nc, byrow = TRUE)
  colnames(unnested_signature) = paste0("X", seq_len(nc))
  unnested_signature = tibble::as_tibble(unnested_signature)

  x["signature"] = NULL

  x = tibble::as_tibble(cbind(x, unnested_signature))
  x_attr$names = names(x)
  attributes(x) = x_attr

  x
}
devtools::load_all()
library(motif)
library(stars)
library(terra)
library(tmap)
landcover = read_stars(system.file("raster/landcover2015.tif", package = "motif"))
coma_output = lsp_signature(landcover, type = "cove", window = 15, normalization = "pdf")
coma_output_sig = lsp_restructure(coma_output)
coma_output_stars = lsp_add_stars(coma_output)
mat = coma_output_stars[[1]]
vals = as.matrix(as.data.frame(lsp_add_stars(coma_output_sig))[-c(1:4)])

slic = run_slic(mat, vals, step = 10, nc = 1, con = TRUE, output_type = TRUE, type = "jensen_shannon")

slic = terra::rast(slic)
terra::ext(slic) = terra::ext(as(coma_output_stars, "Raster"))
terra::crs(slic) = terra::crs(as(coma_output_stars, "Raster"))
slic = sf::st_as_sf(terra::as.polygons(slic, dissolve = TRUE))

plot(slic)

tm_shape(landcover) +
  tm_raster(legend.show = FALSE) +
  tm_shape(slic) +
  tm_borders(col = "red")



# nlcd = rast(system.file("raster/nlcd2011.tif", package = "spDataLarge"))
# nlcd_slic = supercells(nlcd, 40, 10, "jenson_shannon")

nlcd = read_stars(system.file("raster/nlcd2011.tif", package = "spDataLarge"))
coma_output = lsp_signature(nlcd, type = "cove", window = 10, normalization = "pdf", ordered = FALSE)
coma_output_sig = lsp_restructure(coma_output)
coma_output_stars = lsp_add_stars(coma_output)
coma_output_sf = lsp_add_sf(coma_output)
mat = dim(coma_output_stars)[2:1]
vals = as.matrix(as.data.frame(lsp_add_stars(coma_output_sig))[-c(1:4)])

slic = run_slic(mat, vals, step = 5, nc = 0.01, con = TRUE, output_type = TRUE, type = "jensen_shannon")

slic = terra::rast(slic)
terra::ext(slic) = terra::ext(as(coma_output_stars, "Raster"))
terra::crs(slic) = terra::crs(as(coma_output_stars, "Raster"))
slic = sf::st_as_sf(terra::as.polygons(slic, dissolve = TRUE))

tm_shape(nlcd) +
  tm_raster(legend.show = FALSE, drop.levels = TRUE) +
  tm_shape(coma_output_sf) +
  tm_borders(col = "black", alpha = 0.15) +
  tm_shape(slic) +
  tm_borders(col = "red")

