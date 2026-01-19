# supercells: Superpixels of Spatial Data

Creates superpixels based on input spatial data. This package works on
spatial data with one variable (e.g., continuous raster), many variables
(e.g., RGB rasters), and spatial patterns (e.g., areas in categorical
rasters). It is based on the SLIC algorithm (Achanta et al. (2012)
[doi:10.1109/TPAMI.2012.120](https://doi.org/10.1109/TPAMI.2012.120) ),
and readapts it to work with arbitrary dissimilarity measures.

## See also

Useful links:

- <https://jakubnowosad.com/supercells/>

- Report bugs at <https://github.com/Nowosad/supercells/issues>

## Author

**Maintainer**: Jakub Nowosad <nowosad.jakub@gmail.com>
([ORCID](https://orcid.org/0000-0002-1057-3721))

Other contributors:

- Pascal Mettes (Author of the initial C++ implementation of the SLIC
  Superpixel algorithm for image data) \[contributor\]

- Charles Jekel (Author of underlying C++ code for dtw) \[contributor\]
