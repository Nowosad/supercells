# Set stored `sc_slic()` parameters

Writes key
[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md)
parameters to attributes on a supercells object.

## Usage

``` r
sc_slic_set_params(sc, params)
```

## Arguments

- sc:

  An sf object.

- params:

  A data.frame, typically from
  [`sc_slic_get_params()`](https://jakubnowosad.com/supercells/reference/sc_slic_get_params.md).
  Only the first row is used.

## Value

The input object with updated attributes.

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md),
[`sc_slic_get_params()`](https://jakubnowosad.com/supercells/reference/sc_slic_get_params.md)
