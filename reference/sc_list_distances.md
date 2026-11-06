# List supported distance measures

Returns information about the distance measures that can be provided by
name to the `dist_fun` argument. The table describes their input
requirements, what they compare, the effect of rescaling individual
raster layers, and practical guidance.

## Usage

``` r
sc_list_distances()
```

## Value

A data frame with one row per supported measure and columns:

- `name`:

  Value accepted by `dist_fun`.

- `label`:

  Measure name.

- `input`:

  Type of the input vector.

- `compares`:

  What aspect of the vectors the measure compares.

- `effect_of_layer_rescaling`:

  Effect of multiplying individual raster layers by positive constants.

- `guidance`:

  Practical usage guidance.

## Details

Other measures can be supplied explicitly as a function to `dist_fun`.
In that case, the function must take two numeric vectors and return one
finite, nonnegative value, with smaller values meaning greater
similarity.

For `"cosine_distance"`, two zero profiles have distance zero. A zero
profile and a nonzero profile have distance one.

## Examples

``` r
sc_list_distances()
#>                 name                                label
#> 1          euclidean                            Euclidean
#> 2  squared_euclidean                    Squared Euclidean
#> 3          manhattan                            Manhattan
#> 4          chebyshev                            Chebyshev
#> 5           canberra                             Canberra
#> 6         lorentzian                           Lorentzian
#> 7            soergel                              Soergel
#> 8         wavehedges                          Wave Hedges
#> 9        czekanowski                          Czekanowski
#> 10     squared_chord                        Squared chord
#> 11       squared_chi                   Squared chi-square
#> 12             clark                                Clark
#> 13               jsd            Jensen-Shannon divergence
#> 14  non-intersection                     Non-intersection
#> 15      wasserstein1                        Wasserstein-1
#> 16   cosine_distance                      Cosine distance
#> 17               dtw                 Dynamic time warping
#> 18             dtw2d Two-dimensional dynamic time warping
#>                                             input
#> 1                                  numeric values
#> 2                                  numeric values
#> 3                                  numeric values
#> 4                                  numeric values
#> 5                                  numeric values
#> 6                                  numeric values
#> 7                              nonnegative values
#> 8                              nonnegative values
#> 9                              nonnegative values
#> 10                             nonnegative values
#> 11                             nonnegative values
#> 12                             nonnegative values
#> 13                     probabilities summing to 1
#> 14                     probabilities summing to 1
#> 15 probabilities in equal-width bins summing to 1
#> 16                                 numeric values
#> 17                       ordered numeric sequence
#> 18          two equal-length coordinate sequences
#>                                           compares effect_of_layer_rescaling
#> 1                root-sum-square layer differences       can change distance
#> 2                 summed squared layer differences       can change distance
#> 3                summed absolute layer differences       can change distance
#> 4                largest absolute layer difference       can change distance
#> 5                       relative layer differences                 no effect
#> 6              summed log-scaled layer differences       can change distance
#> 7             difference relative to larger values       can change distance
#> 8                summed relative layer differences                 no effect
#> 9           difference relative to total magnitude       can change distance
#> 10    differences after square-root transformation       can change distance
#> 11 squared differences relative to combined values       can change distance
#> 12             combined relative layer differences                 no effect
#> 13         difference between probability profiles  requires renormalization
#> 14  probability assigned differently across layers  requires renormalization
#> 15    probability mass shifted across ordered bins  requires renormalization
#> 16   pattern across layers, ignoring overall scale       can change distance
#> 17             sequence difference after alignment       can change distance
#> 18        2D trajectory difference after alignment       can change distance
#>                                                                                guidance
#> 1                                           Scale layers to comparable units or ranges.
#> 2                                           Scale layers to comparable units or ranges.
#> 3                                           Scale layers to comparable units or ranges.
#> 4                                           Scale layers to comparable units or ranges.
#> 5                                                          Layer scaling has no effect.
#> 6                                           Scale layers to comparable units or ranges.
#> 7                                                  Scale layers to set their influence.
#> 8                                                          Layer scaling has no effect.
#> 9                                                  Scale layers to set their influence.
#> 10                                                 Scale layers to set their influence.
#> 11                                                 Scale layers to set their influence.
#> 12                                                         Layer scaling has no effect.
#> 13                               Normalize every cell to sum to 1 and use mean centres.
#> 14                               Normalize every cell to sum to 1 and use mean centres.
#> 15 Normalize every cell to sum to 1, use mean centres, and arrange layers in bin order.
#> 16                                         Scale layers only to change their influence.
#> 17                                                    Arrange layers in sequence order.
#> 18              Arrange layers by dimension, with sequence order within each dimension.
```
