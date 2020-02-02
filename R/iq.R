#' IQ Quality index
#'
#' @param original Raster brick object with the original multispectral bands
#' @param modified Raster brick object with the bands obtained by the fusion process
#' @param method Resampling method, should be ‘"bilinear"’ for bilinear interpolation, or ‘"ngb"’
#'
#' @return Quality index
#'
#' @importFrom raster nlayers
#'
#' @export
#'
iq <- function(original, modified, method = "bilinear") {
   iq <- rep(NA, nlayers(original) * 4)
   dim(iq) <- c(4, nlayers(original))
   colnames(iq) <- names(original)
   rownames(iq) <- c("Correlation", "Mean", "Standard Deviation", "IQ")
   for (b in 1:nlayers(original)) {
     iq[, b] <- f_iq(original[[b]], modified[[b]], method = "bilinear")
   }
   return(iq)
}
