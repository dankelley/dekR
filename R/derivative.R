#' Numerical derivative, using centred differences for interior points
#'
#' If the lengths of x and y exceed 2, then `derivative()` computes centred
#' derivatives for interior points, and first differences at the endpoints.
#' If the lengths equal 2, then `derivative()` computes the first difference
#' and repeats the value. Otherwise, an error is reported.
#'
#' @param x a numerical vector holding the independent coordinate. An error is
#' reported if `x` holds under 3 points.
#'
#' @param y a numerical vector of the same length as `x`, holding the dependent
#' coordinate, y=y(x).
#'
#' @examples
#'
#' # Demonstrate with a known derivative
#' x <- seq(0, 2*pi, pi/16)
#' y <- sin(x)
#' dydx <- derivative(x, y)
#' layout(rbind(1,2), heights=c(0.6, 0.4))
#' par(mar=c(3,3,1,1), mgp=c(2,0.7,0))
#' plot(x, cos(x), ylab="True, Estimate")
#' points(x, dydx, col=2, pch=3)
#' legend("top", pch=c(1, 3), col=c(1, 2), legend=c("True", "Estimate"))
#' plot(x, cos(x) - dydx, ylab="True - Estimate")
#'
#' @return derivative dy/dx
#' @author Dan Kelley
#' @export
derivative <- function(x, y)
{
    n <- length(x)
    if (length(y) != n)
        stop("x and y must have equal length, but the lengths are ", n, " and ", length(y))
    if (n > 2) {
        c((y[2L]-y[1L]) / (x[2L]-x[1L]),
            diff(y, lag=2L) / diff(x, lag=2L),
            (y[n]-y[n-1L]) / (x[n]-x[n-1L]))
    } else {
        rep((y[2]-y[1]) / (x[2]-x[1]), 2L)
    }
}

