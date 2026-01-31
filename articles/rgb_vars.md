# Supercells of an RGB raster

Superpixels are a collection of segmentation concepts of grouping pixels
with similar characteristics. In this package, we refer to them as
supercells. It is often used in computer vision to delineate parts of
RGB images that are more meaningful and easier to analyze. When applied
to RGB images, each superpixel contains similar colors that also could
represent real-world objects. A large number of methods for creating
superpixels were developed in the last decades, with the SLIC algorithm
(Achanta et al. (2012), <doi:10.1109/TPAMI.2012.120>) being the most
prominent.

The **supercells** package aims to utilize the concept of supercells for
a variety of spatial data. This package works on spatial data with one
variable (e.g., continuous raster), many variables (e.g., RGB rasters),
and spatial patterns (e.g., areas in categorical rasters). Therefore, it
enables not only to find areas that look similar on an RGB (satellite)
image, but also to regionalize areas with comparable values of one or
more variables.

This vignette shows how to use the **supercells** package on an RGB
raster dataset. To reproduce the following results on your own computer,
install and attach the packages:

``` r
library(supercells)    # supercells for spatial data
library(terra)         # spatial raster data reading and handling
library(sf)            # spatial vector data reading and handling
```

The first step is to read the input data. This time our input data
`ortho.tif`, included in the **supercells** package, contains three
layers representing red, green, and blue satellite bands[¹](#fn1).

``` r
ortho = rast(system.file("raster/ortho.tif", package = "supercells"))
plot(ortho)
```

![](rgb_vars_files/figure-html/setup-1.png)

The
[`supercells()`](https://jakubnowosad.com/supercells/reference/supercells.md)
function will be used to delineate areas with similar colors - they
could potentially represent the same objects (e.g., buildings) or land
covers (e.g., grasses). For an RGB image we can use one additional
argument called `transform = "to_LAB"`. This means that internal
calculation will be done on the LAB color space instead on the RGB input
one.

``` r
ortho_slic1 = supercells(ortho, k = 2000, compactness = 10,
                         transform = "to_LAB")
ortho_slic1
```

The `ortho_slic1` output is an `sf` object, where each row stores
supercell id (`supercells`), coordinates of the supercell centers (`x`
and `y`), and an average of all of the input variables (red, green, and
blue colors).

This allows us to create two types of visualizations. The first one is
just an overlay of the supercells borders on top of the original RGB
image.

``` r
plot(ortho)
plot(st_geometry(ortho_slic1), add = TRUE)
```

![](rgb_vars_files/figure-html/unnamed-chunk-4-1.png)

The second visualization type requires converting the average RGB values
into their *hex*adecimal representation first.

``` r
rgb_to_hex = function(x){
  apply(t(x), 2, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
}
avg_colors = rgb_to_hex(st_drop_geometry(ortho_slic1[4:6]))
```

Next, we can plot each superpixel with its average color, but without
the border lines:

``` r
plot(st_geometry(ortho_slic1), border = avg_colors, col = avg_colors)
```

![](rgb_vars_files/figure-html/unnamed-chunk-6-1.png)

Note that the above visualization is not an image - it is a set of
colored supercells. Therefore, instead of representing this area by
87,600 cells, we are using (just) 1,675 supercells.

The `ortho_slic1` object can be next use for clustering similar objects
or labeling them.

## References

Achanta, Radhakrishna, Appu Shaji, Kevin Smith, Aurelien Lucchi, Pascal
Fua, and Sabine Süsstrunk. 2012. “SLIC Superpixels Compared to
State-of-the-Art Superpixel Methods.” *IEEE Transactions on Pattern
Analysis and Machine Intelligence* 34 (11): 2274–82.

------------------------------------------------------------------------

1.  It also has an empty square on the left part of the image for
    testing purposes.
