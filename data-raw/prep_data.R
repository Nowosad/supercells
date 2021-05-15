library(terra)

# var 1 -------------------------------------------------------------------
# https://gist.github.com/btupper/8e8eb8c0ebf4402a3f87b5638eca954a
v = rev(rast(volcano))
ext(v) = c(2667400, 2668010, 6478705, 6479575)
crs(v) = "EPSG:27200"
names(v) = "elevation"
writeRaster(v, "inst/raster/volcano.tif")

# RGB ---------------------------------------------------------------------
library(rgugik)

geodb_download("wielkopolskie", outdir = "./tdata")
reserves = read_sf("data/PL.PZGiK.201.30/BDOO/PL.PZGIK.201.30__OT_ADMS_P.xml")
p = reserves[1, ]
req_df = ortho_request(p)
tile_download(req_df[3, ], outdir = "./tdata")

r = rast("data/69884_519852_N-33-130-D-d-1-2.tif")
r2 = terra::aggregate(r, fact = 8)
r2 = crop(r2, c(358900, 359776.7, 505400, 505800))
mymask = st_sf(geom = st_as_sfc(st_bbox(raster::extent(c(359000, 359100, 505500, 505600)))))
r3 = mask(r2, vect(mymask), inverse = TRUE)

r3_slic2 = supercells(r3, 1000, 1, clean = TRUE, transform = "to_LAB")
# terra::plotRGB(r3, stretch = "lin")
# lines(vect(r3_slic2), add = TRUE, col = "red")

writeRaster(r2, "inst/raster/ortho.tif", overwrite = TRUE, datatype = "INT1U")
unlink("tdata", recursive = TRUE, force = TRUE)

# landsat = rast(system.file("raster/landsat.tif", package = "spDataLarge"))
#
# .linStretch = function (x) {
#   v_vals = terra::values(x)
#   v = stats::quantile(v_vals, c(0.02, 0.98), na.rm = TRUE)
#   temp = (255 * (x - v[1]))/(v[2] - v[1])
#   temp[temp < 0] = 0
#   temp[temp > 255] = 255
#   return(temp)
# }
# landsat2 = landsat
# landsat2[[1]] = .linStretch(landsat[[1]])
# landsat2[[2]] = .linStretch(landsat[[2]])
# landsat2[[3]] = .linStretch(landsat[[3]])
# landsat2[[4]] = .linStretch(landsat[[4]])
# landsat3 = landsat2[[1:3]]
#
# new_bbox = c(301905, 311905, 4111245, 4121245)
#
# landsat4 = crop(landsat3, new_bbox)
# writeRaster(landsat4, "inst/raster/landsat.tif", overwrite = TRUE, datatype = "INT1U")
