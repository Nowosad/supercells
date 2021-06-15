# remotes::install_github("r-prof/jointprof")
library(jointprof)
devtools::load_all()
library(spDataLarge)
library(terra)
library(sf)
library(tmap)

# dog = raster::raster("dog.png"); crs(dog) = "EPSG:2180"
# dog2 = rast("dog.png"); ext(dog2) = c(0, 320, 0, 240); crs(dog2) = "EPSG:2180"
srtm = rast(system.file("raster/srtm.tif", package = "spDataLarge"))

target_file = "Rprof.out"

# Collect profile data
start_profiler(target_file)
## code to be profiled
supercells(srtm, 5000, 1)
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
system2(
  find_pprof(),
  c(
    "-http",
    "localhost:8080",
    shQuote(pprof_target_file)
  )
)
"http://localhost:8080"



# 2 -----------------------------------------------------------------------
srtm = rast(system.file("raster/landsat.tif", package = "spDataLarge"))
target_file = "Rprof.out"

# Collect profile data
start_profiler(target_file)
supercells(srtm, 1000, 1, dist_fun = "dtw", iter = 1)
stop_profiler()

# Convert to pprof format and analyze
pprof_target_file = "Rprof.pb.gz"
profile_data = profile::read_rprof(target_file)
profile::write_pprof(profile_data, pprof_target_file)

paste0(find_pprof(), " -http", " localhost:8080", " 'Rprof.pb.gz'")

# system2(
#   find_pprof(),
#   c(
#     "-http",
#     "localhost:8080",
#     shQuote(pprof_target_file)
#   )
# )
"http://localhost:8080"
