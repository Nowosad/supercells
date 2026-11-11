# Use adaptive compactness mode

Creates a compactness mode object for adaptive compactness. The
`"local_max"` method scales value distances separately for each center
using the largest feature-space distance within its local search window.
The scale is updated before each assignment iteration, while spatial
distances remain scaled by `step`.

## Usage

``` r
use_adaptive(method = "local_max")
```

## Arguments

- method:

  Adaptive compactness method. Currently only `"local_max"` is
  supported.

## Value

An adaptive compactness mode object for `compactness` arguments.

## Examples

``` r
use_adaptive()
#> $method
#> [1] "local_max"
#> 
#> attr(,"class")
#> [1] "sc_adaptive"
```
