# Get stored `sc_slic()` parameters

Returns key
[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md)
parameters stored as attributes on a supercells object.

## Usage

``` r
sc_slic_get_params(sc)
```

## Arguments

- sc:

  An sf object returned by
  [`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md).

## Value

A one-row data.frame with columns: `step`, `compactness`,
`compactness_method`, and `dist_fun`. The `dist_fun` column is
character; custom distance functions are stored as `NA`.

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md),
[`sc_slic_set_params()`](https://jakubnowosad.com/supercells/reference/sc_slic_set_params.md)
