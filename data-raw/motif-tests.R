lsp_restructure = function(x){
  x_attr = attributes(x)
  nc = ncol(x$signature[[1]])

  unnested_signature = matrix(unlist(x$signature, use.names = FALSE),
                              ncol = nc, byrow = TRUE)
  colnames(unnested_signature) = paste0("X", seq_len(nc))
  unnested_signature = tibble::as_tibble(unnested_signature)

  x["signature"] = NULL

  x = tibble::as_tibble(cbind(x, unnested_signature))
  x_attr$names = names(x)
  attributes(x) = x_attr

  x
}
library(motif)
library(stars)
library(tmap)
landcover = read_stars(system.file("raster/landcover2015.tif", package = "motif"))
coma_output = lsp_signature(landcover, type = "cove", window = 200, normalization = "pdf")
coma_output_sig = lsp_restructure(coma_output)
coma_output_stars = lsp_add_stars(coma_output_sig)
coma_output_stars2 = as.matrix(as.data.frame(coma_output_stars)[-c(1:4)])

