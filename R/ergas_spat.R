#' ERGAS Spatial Quality index
#'
#' @param rp Resolution of the original panchromatic band
#' @param rms Resolution of the original multispectral bands
#' @param original Raster brick object with the original multispectral bands
#' @param modified Raster brick object with the bands obtained by the fusion process
#' @param pan Raster layer object with the panchromatic band
#' @param method Resampling method, should be ‘"bilinear"’ for bilinear interpolation, or ‘"ngb"’
#'
#' @return The value of the quality index
#'
#' @import stats
#' @importFrom raster values nlayers raster
#'
#' @export
#'
#' @examples
#' pca <- pca_fusion(mis=mis, pan=pan, method="bilinear", bits=16)
#' ergas_spat(original=mis, modified=pca, pan=pan, method="bilinear")
#'
ergas_spat <- function (original, modified, pan, method,
                        rp = mean(res(pan)), rms = mean(res(original))) {
    resamp <- resample(original, modified, method = method)
    me     <- values(resamp)
    rmse2  <- multi <- c()

    for (b in 1:nlayers(original)) {
       mean_values <- mean(values(raster(resamp, b)), na.rm = TRUE)
       mean_pan    <- mean(values(pan), na.rm = TRUE)
       sd_values   <- sd(values(raster(resamp, b)), na.rm = TRUE)
       sd_pan      <- sd(values(pan), na.rm = TRUE)

       gain1   <- sd_values / sd_pan
       bias    <- mean_values - gain1 * mean_pan
       pan_adj <- pan * gain1 + bias

       rmse2[b] <- mean(values(pan_adj - modified[[b]])^2, na.rm = TRUE)
       multi[b] <- mean(me[, b], na.rm = TRUE)^2
    }
    return(100 * rp * sqrt(mean(rmse2 / multi, na.rm = TRUE)) / rms)
}
