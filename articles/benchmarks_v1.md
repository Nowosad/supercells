# Benchmarks of the supercells package

> **Version note:** This vignette documents the
> [`supercells()`](https://jakubnowosad.com/supercells/reference/supercells.md)
> interface as it existed in version 1.0 of the package. Some arguments,
> defaults, and behaviors may differ in newer releases. For up-to-date
> details, see
> [`?supercells`](https://jakubnowosad.com/supercells/reference/supercells.md)
> and the current reference docs.

``` r
library(terra)
library(sf)
library(supercells)
library(bench)
```

``` r
library(rsi)
my_point = st_sfc(st_point(c(22.850000, 54.250000)), crs = "EPSG:4326")
my_point = st_transform(my_point, crs = "EPSG:2180")
my_poly = st_sf(geom = st_buffer(my_point, dist = 5000))
sa_landsat = get_landsat_imagery(my_poly,
                                 start_date = "2023-09-01", end_date = "2023-10-31",
                                 output_filename = tempfile(fileext = ".tif"))
lnd_low = rast(sa_landsat)
lnd_mid = disagg(lnd_low, 4)
lnd_high = disagg(lnd_mid, 8)
```

``` r
lnd_low
#> class       : SpatRaster 
#> dimensions  : 333, 333, 8  (nrow, ncol, nlyr)
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
```

``` r
results_low
#> # A tibble: 12 × 9
#>    expression  step clean  iter      min   median `itr/sec` mem_alloc `gc/sec`
#>    <bch:expr> <dbl> <lgl> <dbl> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#>  1 s             10 TRUE      1    267ms 382.88ms     2.68     62.2MB     8.85
#>  2 s             30 TRUE      1    194ms 328.96ms     3.24     57.6MB     9.08
#>  3 s            100 TRUE      1    179ms 196.08ms     4.28     56.7MB     9.85
#>  4 s             10 FALSE     1    524ms 569.85ms     1.67     82.1MB     6.53
#>  5 s             30 FALSE     1    238ms 319.95ms     3.16     60.8MB     7.59
#>  6 s            100 FALSE     1    175ms 177.26ms     4.92     57.5MB     8.86
#>  7 s             10 TRUE     10    637ms 664.22ms     1.46     60.1MB     2.78
#>  8 s             30 TRUE     10    572ms  585.1ms     1.56     57.7MB     2.96
#>  9 s            100 TRUE     10    490ms 505.76ms     1.80     56.8MB     3.78
#> 10 s             10 FALSE    10    916ms    1.03s     0.990    82.4MB     3.86
#> 11 s             30 FALSE    10    624ms 644.35ms     1.47     61.1MB     2.79
#> 12 s            100 FALSE    10    476ms 486.17ms     2.00     57.9MB     3.40
```

``` r
lnd_mid
#> class       : SpatRaster 
#> dimensions  : 1332, 1332, 8  (nrow, ncol, nlyr)
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
        iterations = 10
        )
    }
)
```

``` r
results_mid
#> # A tibble: 12 × 9
#>    expression  step clean  iter      min   median `itr/sec` mem_alloc `gc/sec`
#>    <bch:expr> <dbl> <lgl> <dbl> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#>  1 s             10 TRUE      1    3.64s     3.8s    0.262   936.35MB    1.54 
#>  2 s             30 TRUE      1    2.53s     2.6s    0.380   908.11MB    1.75 
#>  3 s            100 TRUE      1    2.49s    2.72s    0.370   902.22MB    2.41 
#>  4 s             10 FALSE     1    8.86s    9.15s    0.108     1.36GB    1.94 
#>  5 s             30 FALSE     1    2.48s    2.67s    0.371   949.17MB    1.41 
#>  6 s            100 FALSE     1    1.79s    2.13s    0.471   906.74MB    1.84 
#>  7 s             10 TRUE     10   10.99s   11.05s    0.0904  940.44MB    0.443
#>  8 s             30 TRUE     10    9.62s    9.75s    0.103   909.09MB    0.329
#>  9 s            100 TRUE     10    9.37s    9.48s    0.106   902.34MB    0.327
#> 10 s             10 FALSE    10    15.9s   16.58s    0.0604    1.37GB    1.33 
#> 11 s             30 FALSE    10    9.82s    9.87s    0.101   951.88MB    0.456
#> 12 s            100 FALSE    10    8.57s    9.02s    0.111    907.1MB    0.411
```

``` r
lnd_high
#> class       : SpatRaster 
#> dimensions  : 10656, 10656, 8  (nrow, ncol, nlyr)
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
#> 1 s             3.54m    3.61m   0.00432    59.7GB    0.214
```
