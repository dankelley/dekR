#' Numerical derivative, using centred differences for interior points
#'
#' This computes the numerical derivative of y=y(x), using 
#' centred derivatives for interior points, and first differences for
#' the endpoints.  (If there are only 2 points, then first differences
#' are used, then repeated so the result has the same length as `x`.)
#' If `y` is not provided, then the derivative of `x` with respect to
#' index is computed.
#'
#' @param x a numerical vector.
#'
#' @param y optional numerical vector of the same length as `x`. See
#' \sQuote{Details} for behaviour if `y` is not given.
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
    if (n < 1L)
        stop("Must provide x")
    if (missing(y)) {
        y <- x
        x <- seq_len(n)
    } else {
        if (length(y) != n)
            stop("x and y must have equal length, but the lengths are ", n, " and ", length(y))
    }
    if (n > 2) {
        c((y[2L]-y[1L]) / (x[2L]-x[1L]),
            diff(y, lag=2L) / diff(x, lag=2L),
            (y[n]-y[n-1L]) / (x[n]-x[n-1L]))
    } else {
        rep((y[2]-y[1]) / (x[2]-x[1]), 2L)
    }
}

