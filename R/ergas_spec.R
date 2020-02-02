#' ERGAS Spectral Quality index
#'
#' @param original Raster brick object with the original multispectral bands
#' @param modified Raster brick object with the bands obtained by the fusion process
#' @param pan Raster object with the panchromatic band.
#' @param method Resampling method, should be ‘"bilinear"’ for bilinear interpolation, or ‘"ngb"’
#' @param rp Resolution of the original panchromatic band
#' @param rms Resolution of the original multispectral bands
#'
#' @return Quality index
#'
#' @import stats
#' @importFrom raster values
#'
#' @export
#'
#' @examples
#' pca <- pca_fusion(mis=mis, pan=pan, method="bilinear", bits=16)
#' ergas_spec(original=mis, modified=pca, pan=pan, method="bilinear")
#'
ergas_spec <- function (original, modified, pan, method,
                        rp=mean(res(pan)), rms=mean(res(original))) {

   resamp <- resample(original, modified, method = method)
   me     <- values(resamp)
   fus    <- values(modified)
   rmse2  <- multi <- c()

   for (b in 1:dim(me)[2]) {
      rmse2[b] <- mean((me[, b] - fus[, b])^2, na.rm = TRUE)
      multi[b] <- mean(me[, b], na.rm = TRUE)^2
   }

   return(100 * rp * sqrt(mean(rmse2 / multi, na.rm = TRUE)) / rms)
}
