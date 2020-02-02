#' High Pass Filter Transformation
#'
#' @param mis Raster brick object with the original multispectral bands
#' @param pan Raster layer object with the panchromatic band
#' @param method Resampling method, should be ‘"bilinear"’ for bilinear interpolation, or ‘"ngb"’
#' @param bits Radiometric resolution of the original multispectral bands
#'
#' @return Raster brick object with the bands obtained by the fusion process
#'
#' @import stats
#' @importFrom raster res focal trim cellStats resample values values<- nlayers raster brick setValues
#'
#' @export
#'
#' @examples
#' hpf <- hpf_fusion(mis=mis, pan=pan, method="bilinear", bits=16)
#'
hpf_fusion <- function(mis, pan, method="bilinear", bits) {

   r_res <- bits

   # Get the filter matrix
   res_ratio <- res(mis)[1] / res(pan)[1]  #(eq. 1)
   matriz <- get_matrix(res_ratio)
   w <- matriz$hpm

   # Filtering and trimming
   pan_fil <- focal(pan, w)
   pan_trim <- trim(pan_fil)

   # Get standard deviation of filtered-trimmed panchromatic image
   sd_pan <- cellStats(pan_trim, "sd")

   # For each layer in panchromatic:
   #   1: resample, 2: calculate sd, 3: calculate w (eq. 2),
   #   4: calculate Pout layer (eq. 3), 5: calculate gain and bias (eq. 5),
   #   6: calculate Pfus (eq. 4), 7: Level correction
   pfus <- list()
   rr <- 2^r_res - 1

   for (b in 1:nlayers(mis)) {
      capa <- resample(raster(mis, layer = b), pan_trim, method = method)
      sd_mis <- cellStats(capa, "sd")
      w <- matriz$m * sd_mis / sd_pan
      pout <- capa + w * pan_trim
      gain <- sd_mis / cellStats(pout, "sd")
      bias <- cellStats(capa, "mean") - gain * cellStats(pout, "mean")

      pfus[[b]] <- round(pout * gain + bias)
      values(pfus[[b]])[which(values(pfus[[b]]) < 0)]  <- 0
      values(pfus[[b]])[which(values(pfus[[b]]) > rr)]  <- rr
   }
   return(brick(pfus))
}
