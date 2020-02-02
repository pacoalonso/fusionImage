#' IQ helpler function
#'
#' @param original Raster layer object with the original image
#' @param modified Raster layer object with the modified image
#' @param method Resampling method, should be ‘"bilinear"’ for bilinear interpolation, or ‘"ngb"’
#'
#' @return Raster brick object with the bands obtained by the fusion process
#'
#' @import stats
#' @importFrom raster values
#'
f_iq <- function(original, modified, method = "bilinear") {
    resamp <- resample(original, modified, method = method)
    b11 <- values(resamp)
    b12 <- values(modified)
    w <- which(!is.na(b11) & !is.na(b12))
    r <- cor(b11[w], b12[w])
    m <- 2 * mean(b11[w]) * mean(b12[w]) / (mean(b11[w])^2 + mean(b12[w])^2)
    s <- 2 * sd(b11[w]) * sd(b12[w]) / (var(b11[w]) + var(b12[w]))
    q <- r * m * s
    iq <- c(r, m, s, q)
    names(iq) <- c("Correlation", "Mean", "Standard Deviation", "IQ")
    return(iq)
}
