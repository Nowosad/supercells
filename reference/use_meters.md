# Mark step values as meters

Creates a units value in meters for use in `step` arguments. Use plain
numerics for cell units, and `use_meters()` for map-distance steps.

## Usage

``` r
use_meters(x)
```

## Arguments

- x:

  A single positive numeric value.

## Value

A
[units::units](https://r-quantities.github.io/units/reference/units.html)
object in meters (`m`).

## Examples

``` r
use_meters(100)
#> 100 [m]
```
