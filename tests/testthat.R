library(testthat)
library(supercells)
library(terra)

# 1 var
# https://gist.github.com/btupper/8e8eb8c0ebf4402a3f87b5638eca954a
v = rev(rast(volcano))
ext(v) = c(2667400, 2668010, 6478705, 6479575)
crs(v) = "EPSG:27200"
names(v) = "elevation"

# ttm()
# qtm(v) + tm_grid()


