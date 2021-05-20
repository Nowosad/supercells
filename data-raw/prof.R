library(jointprof)
devtools::load_all()
library(spDataLarge)
library(terra)
library(sf)
library(tmap)

input = rast(system.file("raster/srtm.tif", package = "spDataLarge"))
input = rast(system.file("raster/ortho.tif", package = "supercells"))

target_file = "Rprof.out"

# Collect profile data
start_profiler(target_file)
## code to be profiled
supercells(input, k = 2000, compactness = 10)
stop_profiler()

# Analyze profile data
# summaryRprof(target_file)
# profvis::profvis(prof_input = target_file)
# proftools::readProfileData(target_file)
# prof.tree::prof.tree(target_file)

# Convert to pprof format and analyze
pprof_target_file = "Rprof.pb.gz"
profile_data = profile::read_rprof(target_file)
profile::write_pprof(profile_data, pprof_target_file)
# system2(
#   find_pprof(),
#   c(
#     "-http",
#     "localhost:8081",
#     shQuote(pprof_target_file)
#   )
# )
# "http://localhost:8080"
