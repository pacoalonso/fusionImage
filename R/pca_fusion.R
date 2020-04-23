#' Principal Components Analysis Transformation
#'
#' @param mis Raster brick object with the original multispectral bands
#' @param pan Raster layer object with the panchromatic band
#' @param method Resampling method, should be ‘"bilinear"’ for bilinear interpolation, or ‘"ngb"’
#' @param bits Radiometric resolution of the original multispectral bands
#' @param matrix PCA matrix. It is the pca matrix of variable loadings. The usual working mode is to allow the function calculates it, but it can be introduced by the user. 
#' @param mode With mode=-1 the first component is multiplied by -1
#'
#' @return Raster brick object with the bands obtained by the fusion process
#'
#' @import stats
#' @importFrom raster res focal trim cellStats resample values values<- nlayers brick setValues
#'
#' @export
#'
#' @examples
#' pca <- pca_fusion(mis=mis, pan=pan, method="bilinear", bits=16)
#'
pca_fusion <- function (mis, pan, method, bits, matrix = NULL, mode = 1) {
    r_res <- bits
    capa <- list()
    for (b in 1:nlayers(mis)) {
        capa[[b]] <- resample(mis[[b]], pan, method = method)
    }
    pixels <- nrow(capa[[1]]) * ncol(capa[[1]])
    mx <- matrix(rep(NA, nlayers(mis) * pixels), ncol = nlayers(mis))
    for (b in 1:nlayers(mis)) mx[, b] <- values(capa[[b]])
    if (!is.null(matrix)) {
        rotacion <- matrix
    } else {
        pca <- prcomp(na.omit(mx), center = TRUE)
        rotacion <- t(pca$rotation)
        if (mode == -1 & sum(sign(rotacion[1, ])) < 0) {
           rotacion[1, ] <- -1 * rotacion[1, ]
        }
    }
    cat("PCA Rotation:\n\n")
    for (b in 1:nlayers(mis)) cat(paste0("PC", b), rotacion[b, ], "\n")
    if (sum(sign(rotacion[1, ])) < 0) {
       cat("\nWARNING: First component negatively correlated with bands.\n")
       cat("Use mode=-1\n")
    }
    prediccion <- matrix(rep(0, nlayers(mis) * pixels), ncol = nlayers(mis))
    for (j in 1:nlayers(mis)) {
       for (i in 1:nlayers(mis)) {
          prediccion[, i] <- prediccion[, i] + rotacion[i, j] * values(capa[[j]])
       }
    }
    pc1 <- prediccion[, 1]
    gain <- sd((pc1), na.rm = T) / sd(values(pan), na.rm = T)
    bias <- mean((pc1), na.rm = T) - gain * mean(values(pan), na.rm = T)
    pan_ajust <- pan * gain + bias
    prediccion[, 1] <- values(pan_ajust)

    xx <- matrix(rep(0, nlayers(mis) * pixels), ncol = nlayers(mis))
    for (i in 1:nlayers(mis)) {
       for (b in 1:nlayers(mis)) {
          xx[, i] <- xx[, i] + rotacion[b, i] * prediccion[, b]
       }
    }
    capa2 <- capa
    rr <- 2^r_res - 1
    for (b in 1:nlayers(mis)) {
        values(capa2[[b]]) <- round(xx[, b])
        values(capa2[[b]]) [which(values(capa2[[b]])  < 0)]  <- 0
        values(capa2[[b]]) [which(values(capa2[[b]])  > rr)]  <- rr
    }
    return(brick(capa2))
}
