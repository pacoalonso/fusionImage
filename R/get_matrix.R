#' Calculate filtering matrix for High Pass Filter. Internally used by HPF_fusion function
#'
#' @param r A number with the size of the matriz
#'
#' @return A matrix object with the filtering matrix whose size depends on r
#'
#' @export
#'
#' @examples
#' get_matrix(3)
#'
get_matrix <- function(r) {
    if (r >= 1 & r < 2.5) {
        hpk <- 5
        vc <- 24
        m <- 0.25
    }
    if (r >= 2.5 & r < 3.5) {
        hpk <- 7
        vc <- 48
        m <- 0.5
    }
    if (r >= 3.5 & r < 5.5) {
        hpk <- 9
        vc <- 80
        m <- 0.5
    }
    if (r >= 5.5 & r < 7.5) {
        hpk <- 11
        vc <- 120
        m <- 0.65
    }
    if (r >= 7.5 & r < 9.5) {
        hpk <- 13
        vc <- 168
        m <- 1
    }
    if (r >= 9.5) {
        hpk <- 15
        vc <- 336
        m <- 1.35
    }
    hpm <- matrix(rep(-1, hpk^2), ncol = hpk)
    hpm[ceiling(hpk / 2), ceiling(hpk / 2)] <- vc
    return(list(hpk = hpk, vc = vc, m = m, hpm = hpm))
}
