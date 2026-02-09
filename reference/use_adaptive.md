# Use adaptive compactness mode

Creates a compactness mode object for adaptive compactness. The
`"local_max"` method corresponds to SLIC0-style local scaling, where
compactness is adapted using local maximum value distances.

## Usage

``` r
use_adaptive(method = "local_max")
```

## Arguments

- method:

  Adaptive compactness method. Currently only `"local_max"` is supported
  (SLIC0-style).

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
