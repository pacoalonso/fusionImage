#' Gram-Schmidt Transformation
#'
#' @param mis Raster brick object with the original multispectral bands
#' @param pan Raster layer object with the panchromatic band
#' @param method Resampling method, should be ‘"bilinear"’ for bilinear interpolation, or ‘"ngb"’
#' @param bits Radiometric resolution of the original multispectral bands
#'
#' @return Raster brick object with the bands obtained by the fusion process
#'
#' @import stats
#' @importFrom raster res focal trim cellStats resample values nlayers raster layerStats brick setValues mean
#'
#' @export
#'
#' @examples
#' gs <- gs_fusion(mis=mis, pan=pan, method="bilinear", bits=16)
#'
gs_fusion <- function(mis, pan, method="bilinear", bits) {

   r_res <- bits
   gs <- gs_inv <- bt <- nb <- list()

   # Step 1 Simulated panchromatic
   gs[[1]] <- mean(mis)

   # Step 2 Gram-Schmidt transformation
   for (b in 1:nlayers(mis)) {
     k <- b + 1
     rs <- raster(mis, layer = b)
     gs[[k]] <- rs - cellStats(rs, "mean")
     for (j in 1:b) {
       cv <- layerStats(brick(raster(mis, layer = b), gs[[j]]), "cov", na.rm = TRUE)
       cat(dim(cv), "\n")
       phi <- cv$covariance[1, 2] / cv$covariance[2, 2]
       gs[[k]] <- gs[[k]] - phi * gs[[j]]
     }
   }

   # Step 3 Adjust the pan layer to have the same mean and sd than GS[[1]]
   gain <- sd(values(gs[[1]]), na.rm = TRUE) / sd(values(pan), na.rm = TRUE)
   bias <- mean(values(gs[[1]]), na.rm = TRUE) - gain * mean(values(pan), na.rm = TRUE)
   pan <- pan * gain + bias

   # Paso 4a: Resampling of layers GS2 a GS(N+1) and originals
   gs_inv[[1]] <- pan
   for (b in 1:nlayers(mis)) {
      gs_inv[[b+1]] <- resample(gs[[b+1]], pan, method = method)
   }
   bt <- resample(mis, pan, method = method)

   # Step 5: Calculation of GS2 a GS(N+1) adding the average of the original band
   for (b in 1:nlayers(mis)) {
      k <- b + 1
      rs <- raster(bt, layer = 1)
      nb[[k]] <- gs_inv[[k]] + cellStats(rs, "mean")
   }

   # Step 6: Inverse of GS
   rr <- 2^r_res - 1
   for (b in 1:nlayers(mis)) {
      rs <- raster(mis, layer = b)
      for (j in 1:b) {
         cv <- layerStats(brick(rs, gs[[j]]), "cov", na.rm = TRUE)
         phi <- cv$covariance[1, 2] / cv$covariance[2, 2]
         if (j == 1) {
            suma <- gs_inv[[j]] * phi
         } else {
            suma <- suma + phi * gs_inv[[j]]
         }
      }
      v <- round(gs_inv[[b+1]] + cellStats(rs, "mean", na.rm = TRUE) + suma)
      v[which(values(v) < 0)]  <- 0
      v[which(values(v) > rr)]  <- rr
      nb[[b]] <- v
   }

   kk <- brick(nb[-length(nb)])
   return(kk)
}
