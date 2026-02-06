# Benchmarks of the supercells package

``` r
library(terra)
library(sf)
library(supercells)
library(bench)
terra::terraOptions(tempdir = "~/tmp/")
```

``` r
library(rsi)
if (file.exists(low_path)) {
  lnd_low = rast(low_path)
} else {
  my_point = st_sfc(st_point(c(22.850000, 54.250000)), crs = "EPSG:4326")
  my_point = st_transform(my_point, crs = "EPSG:2180")
  my_poly = st_sf(geom = st_buffer(my_point, dist = 5000))
  sa_landsat = get_landsat_imagery(my_poly,
                                   start_date = "2023-09-01", end_date = "2023-10-31",
                                   output_filename = tempfile(fileext = ".tif"))
  lnd_low = rast(sa_landsat)
  if (interactive() && dir.exists("data-raw")) {
    writeRaster(lnd_low, data_raw_low_path, overwrite = TRUE, filetype = "COG")
  }
}
lnd_mid = disagg(lnd_low, 4)
lnd_high = disagg(lnd_mid, 8)
#> 
|---------|---------|---------|---------|
=
                                          
```

``` r
lnd_low
#> class       : SpatRaster 
#> size        : 333, 333, 8  (nrow, ncol, nlyr)
#> resolution  : 30, 30  (x, y)
#> extent      : 745716.1, 755706.1, 711383.8, 721373.8  (xmin, xmax, ymin, ymax)
#> coord. ref. : ETRF2000-PL / CS92 (EPSG:2180) 
#> source      : lnd_low.tif 
#> names       :           A,        B,         G,         R,         N,      S1, ... 
#> min values  : -0.00237125, 0.006140, 0.0142250, 0.0096875, 0.0335025, 0.01516, ... 
#> max values  :  0.09818250, 0.112455, 0.1661212, 0.2097775, 0.5576800, 0.41105, ...
```

``` r
results_low = bench::press(
    step = c(10, 30, 100),
    clean = c(TRUE, FALSE),
    iter = c(1, 10),
    {
    bench::mark(
        s = supercells(lnd_low, compactness = 1, step = step, clean = clean, iter = iter),
        iterations = 10
        )
    }
)
#> Running with:
#>     step clean  iter
#>  1    10 TRUE      1
#>  2    30 TRUE      1
#>  3   100 TRUE      1
#>  4    10 FALSE     1
#>  5    30 FALSE     1
#>  6   100 FALSE     1
#>  7    10 TRUE     10
#>  8    30 TRUE     10
#>  9   100 TRUE     10
#> 10    10 FALSE    10
#> 11    30 FALSE    10
#> 12   100 FALSE    10
```

``` r
results_low
#> # A tibble: 12 × 9
#>    expression  step clean  iter      min   median `itr/sec` mem_alloc `gc/sec`
#>    <bch:expr> <dbl> <lgl> <dbl> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#>  1 s             10 TRUE      1  102.3ms    112ms      7.18    21.2MB    19.1 
#>  2 s             30 TRUE      1   74.9ms   76.1ms     13.1     17.7MB    30.6 
#>  3 s            100 TRUE      1   62.2ms   65.4ms     15.0     16.8MB    15.0 
#>  4 s             10 FALSE     1  292.5ms  316.2ms      3.09    22.6MB     4.95
#>  5 s             30 FALSE     1  109.8ms    116ms      8.74    18.9MB    20.4 
#>  6 s            100 FALSE     1   67.5ms   74.1ms     13.8     17.6MB    13.8 
#>  7 s             10 TRUE     10  325.9ms  350.4ms      2.86    20.2MB     5.44
#>  8 s             30 TRUE     10  285.7ms  297.6ms      3.04    17.9MB     5.47
#>  9 s            100 TRUE     10  217.3ms  243.9ms      3.32    16.9MB     4.31
#> 10 s             10 FALSE    10  502.9ms  552.6ms      1.71    22.2MB     5.48
#> 11 s             30 FALSE    10  318.6ms  339.5ms      2.89    19.1MB     5.79
#> 12 s            100 FALSE    10  220.8ms  232.1ms      3.72      18MB     5.21
```

``` r
lnd_mid
#> class       : SpatRaster 
#> size        : 1332, 1332, 8  (nrow, ncol, nlyr)
#> resolution  : 7.5, 7.5  (x, y)
#> extent      : 745716.1, 755706.1, 711383.8, 721373.8  (xmin, xmax, ymin, ymax)
#> coord. ref. : ETRF2000-PL / CS92 (EPSG:2180) 
#> source      : lnd_mid.tif 
#> names       :           A,        B,         G,         R,         N,      S1, ... 
#> min values  : -0.00237125, 0.006140, 0.0142250, 0.0096875, 0.0335025, 0.01516, ... 
#> max values  :  0.09818250, 0.112455, 0.1661212, 0.2097775, 0.5576800, 0.41105, ...
```

``` r
results_mid = bench::press(
    step = c(10, 30, 100),
    clean = c(TRUE, FALSE),
    iter = c(1, 10),
    {
    bench::mark(
        s = supercells(lnd_mid, compactness = 1, step = step, clean = clean, iter = iter),
        iterations = 1
        )
    }
)
#> Running with:
#>     step clean  iter
#>  1    10 TRUE      1
#>  2    30 TRUE      1
#>  3   100 TRUE      1
#>  4    10 FALSE     1
#>  5    30 FALSE     1
#>  6   100 FALSE     1
#>  7    10 TRUE     10
#>  8    30 TRUE     10
#>  9   100 TRUE     10
#> 10    10 FALSE    10
#> 11    30 FALSE    10
#> 12   100 FALSE    10
```

``` r
results_mid
#> # A tibble: 12 × 9
#>    expression  step clean  iter      min   median `itr/sec` mem_alloc `gc/sec`
#>    <bch:expr> <dbl> <lgl> <dbl> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#>  1 s             10 TRUE      1    2.14s    2.14s     0.466     300MB     2.33
#>  2 s             30 TRUE      1    1.32s    1.32s     0.759     272MB     3.04
#>  3 s            100 TRUE      1    1.17s    1.17s     0.858     266MB     4.29
#>  4 s             10 FALSE     1    6.46s    6.46s     0.155     325MB     6.50
#>  5 s             30 FALSE     1    1.57s    1.57s     0.638     275MB     4.47
#>  6 s            100 FALSE     1     1.3s     1.3s     0.768     267MB     3.07
#>  7 s             10 TRUE     10    5.81s    5.81s     0.172     300MB     8.77
#>  8 s             30 TRUE     10    4.79s    4.79s     0.209     273MB     8.55
#>  9 s            100 TRUE     10    5.55s    5.55s     0.180     266MB     6.85
#> 10 s             10 FALSE    10     9.7s     9.7s     0.103     325MB     4.85
#> 11 s             30 FALSE    10    5.45s    5.45s     0.183     276MB     4.40
#> 12 s            100 FALSE    10     4.9s     4.9s     0.204     268MB     4.08
```

``` r
lnd_high
#> class       : SpatRaster 
#> size        : 10656, 10656, 8  (nrow, ncol, nlyr)
#> resolution  : 0.9375, 0.9375  (x, y)
#> extent      : 745716.1, 755706.1, 711383.8, 721373.8  (xmin, xmax, ymin, ymax)
#> coord. ref. : ETRF2000-PL / CS92 (EPSG:2180) 
#> source      : lnd_high.tif 
#> names       :           A,        B,         G,         R,         N,      S1, ... 
#> min values  : -0.00237125, 0.006140, 0.0142250, 0.0096875, 0.0335025, 0.01516, ... 
#> max values  :  0.09818250, 0.112455, 0.1661212, 0.2097775, 0.5576800, 0.41105, ...
```

``` r
results_high = bench::mark(
        s = supercells(lnd_high, compactness = 1, step = 30, clean = FALSE, iter = 1),
        iterations = 10
)
```

``` r
results_high
#> # A tibble: 1 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 s             2.05m    2.12m   0.00776    16.9GB     1.25
```
