# The supercells() function

The main function in this package (Nowosad and Stepinski 2022) is called
[`supercells()`](https://jakubnowosad.com/supercells/reference/supercells.md).
An overview of its arguments is shown in the table below.

This function expects raster data with one or more layers (representing,
for example, different bands, variables, or dates) in the form of a
**terra**’s `SpatRaster` object or a **stars**’s `stars` object. The
resulting superpixels are stored as `sf` polygons with four or more
columns containing identification numbers of superpixels, y and x
coordinates of their centroids, and one or more columns with average
values of variables for each superpixel.

The number of resulting superpixels can be specified with either `k` or
`step` arguments. `k` relates to the desired number of superpixels. When
the `k` value is set, the algorithm automatically calculates the value
of `step`. As, by default, cluster centroids are located regularly, the
resulting number of superpixels may slightly differ from the number
provided as `k`, e.g., a square raster cannot be divided into five equal
square areas. `step` is the expected distance, in the number of cells,
between initial superpixels’ centroids. This parameter also defines a
zone of influence of each cluster center ($2S \times 2S$ region).

In our software, it is also possible to provide a set of points (an `sf`
object) as `k` together with the `step` value. This way, custom cluster
centers are used, with each of them attracting cells in the surrounding
$2S \times 2S$ region.

While `k` or `step` determines the number of superpixels, the
`compactness` argument controls their spatial shape, with its large
value giving more importance to spatial distances between cells and
superpixels’ centers and its smaller value putting more weight
(importance) to the value distance. The impact of the compactness value
depends on the range of input cell values and the selected distance
measure.

| argument          | description                                                                                                                                                                   |
|:------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **`x`**           | An object of class SpatRaster (terra) or class stars (stars)                                                                                                                  |
| **`k`**           | A number of superpixels desired by the user. It is also possible to provide a set of points (an sf object) as k together with the step value to create custom cluster centers |
| **`compactness`** | A compactness value. Larger values cause clusters to be more compact/even (squarish)                                                                                          |
| `dist_fun`        | A distance function used                                                                                                                                                      |
| `avg_fun`         | An averaging function specifying how the values of the superpixels’ centers are recalculated                                                                                  |
| `clean`           | A boolean specifying if the additional process of connectivity enforcement should be performed                                                                                |
| `iter`            | A number of iterations performed to create the final output                                                                                                                   |
| `step`            | A distance, in the number of cells, between initial superpixels’ centers                                                                                                      |
| `minarea`         | A minimal size of the output superpixels in cells                                                                                                                             |
| `chunks`          | A boolean or numeric value specifying if the input (x) should be split into chunks before deriving superpixels                                                                |
| `future`          | A boolean specifying if the future package should be used for parallelization of the calculations                                                                             |
| `verbose`         | An integer specifying the verboseness of text messages printed during calculations                                                                                            |

Table: Overview of the arguments accepted by the
[`supercells()`](https://jakubnowosad.com/supercells/reference/supercells.md)
function

By default, the `supercells` function behaves accordingly to the
original algorithm described by Achanta et al. (2012) with Euclidean
distance used to calculate the distance between values, and the
arithmetic mean used to calculate an average value of each superpixel.
However, in this package, both of the above parameters can be customized
with the `dist_fun` and `avg_fun` arguments.

The role of the `dist_fun` argument is to specify the distance function
used to obtain the distance between values. It can be done with one of
three mechanisms. The first is to use one of internal C++ functions,
such as `"euclidean"`, `"jsd"` (the Jensen-Shannon distance, Lin
(1991)), and `"dtw"` (dynamic time warping). The second mechanism allows
selecting one of 46 distance and similarity measures implemented in the
R package `philentropy` (Dorst 2018). Thirdly, this argument also
accepts any user-defined R function that returns one value based on
provided two vectors.

The `avg_fun` function, on the other hand, specifies how the values of
the superpixels’ centers are calculated. It has two internal functions
implemented in C++ - `"mean"` and `"median"`, but also accepts any
fitting R function such as,
[`base::mean()`](https://rdrr.io/r/base/mean.html) or
`psych::geometric.mean()`. This also allows providing any other
user-defined R function that returns one value based on an R vector. For
example, the `"median"` function can be used when our input data is
categorical.

Due to its simplicity, the SLIC algorithm does not consider spatial
connectivity directly. For example, it makes it possible to create a
superpixel consisting of two or more distinct patches, where one patch
is large while additional ones are distinctly smaller. The `supercells`
function has two arguments, `clean` and `minarea`, allowing to enforce
connectivity of the superpixels. The first argument, `clean`, is a
boolean (`TRUE` by default) specifying if the connectivity should be
enforced by removing small disconnected patches (by merging with larger
neighborhood ones) or promoting small patches to new superpixels. The
role of the second one, `minarea`, is to control how large disconnected
patches must be not be removed. By default, when `clean = TRUE`, an
average area ($A$) is calculated automatically based on the total number
of cells divided by a number of superpixels, and next, the minimal size
of a superpixel equals $A/\left( 2^{2} \right)$. Alternatively, users
can also specify the minimum size of a superpixel by themselves by
providing an area in the unit of a number of cells in the `minarea`
argument.

The next argument is `iter` - specifying the number of iterations to
create the output. Its default value of 10 follows the advice by Achanta
et al. (2012) This argument defines how many times superpixels (and
their centers) are recalculated before obtaining the final results.

By default, the `supercells` function, as R language, is single-threaded
– runs only on a single thread on the CPU and reads the input raster
values into the computer memory. These features may limit the function’s
usability for raster datasets with millions or more cells and many
variables, for which calculations can either be too slow or require more
memory that is available. To overcome the aforementioned issues, the
`supercells` function has two related optional arguments - `chunks` and
`future`.

The `chunks` argument is set by default to `FALSE`. However, when it is
either `TRUE` or some numerical value, then the split, apply, combine
procedure is used: input raster is divided into several chunks, each
chunk is read into RAM independently and has a set of superpixels
derived; this process is repeated for every chunk, and all the results
are combined into one final object. When the user sets this argument to
`TRUE`, the chunks’ sizes are calculated automatically, while when the
user provides a numerical value to this argument, then the input raster
data is split into chunks with user-defined side length (in the number
of pixels). The `future` argument is also set by default to `FALSE`. If
it is `TRUE`, the user also needs to specify how parallel processing
should be performed and on how many CPU threads with
[`future::plan()`](https://future.futureverse.org/reference/plan.html).

The final argument is called `verbose`, which takes an integer value of
0 or larger, where 0 means no additional messages during the
calculations, and 1 provides basic messages about the current
calculation stage.

## References

Achanta, Radhakrishna, Appu Shaji, Kevin Smith, Aurelien Lucchi, Pascal
Fua, and Sabine Süsstrunk. 2012. “SLIC Superpixels Compared to
State-of-the-Art Superpixel Methods.” *IEEE Transactions on Pattern
Analysis and Machine Intelligence* 34 (11): 2274–82.

Dorst, HG. 2018. “Philentropy: Information Theory and Distance
Quantification with R.” *Journal of Open Source Software* 3 (26): 765.

Lin, Jianhua. 1991. “Divergence Measures Based on the Shannon Entropy.”
*IEEE Transactions on Information Theory* 37 (1): 145–51.
<https://doi.org/djxkkh>.

Nowosad, Jakub, and Tomasz F. Stepinski. 2022. “Extended SLIC
Superpixels Algorithm for Applications to Non-Imagery Geospatial
Rasters.” *International Journal of Applied Earth Observation and
Geoinformation* 112 (August): 102935.
<https://doi.org/10.1016/j.jag.2022.102935>.
