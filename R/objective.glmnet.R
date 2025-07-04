#' Penalized log-likelihood for a fitted glmnet model
#'
#' Computes the penalized log-likelihood sequence for a fitted `glmnet`
#' object. The value for each lambda is obtained from the stored
#' deviance and the elastic net penalty used in the fit.
#'
#' @param object fitted \code{glmnet} object.
#' @return Numeric vector of penalized log-likelihood values, one for each
#'   value of \code{lambda} in the model.
#' @examples
#' x <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' fit <- glmnet(x, y)
#' objective.glmnet(fit)
#' @export objective.glmnet
objective.glmnet <- function(object) {
    dev_ratio <- object$dev
    if (is.null(dev_ratio)) dev_ratio <- object$dev.ratio
    nulldev <- object$nulldev
    lambda <- object$lambda
    if (is.null(dev_ratio) || is.null(nulldev) || is.null(lambda))
        stop("object does not contain deviance information")
    dev_val <- (1 - dev_ratio) * nulldev / 2
    alpha <- 1.0
    if (!is.null(object$call$alpha)) {
        alpha <- eval(object$call$alpha)
    }
    vp <- NULL
    if (!is.null(object$call$penalty.factor)) {
        vp <- eval(object$call$penalty.factor)
    } else if (!is.null(object$warm_fit$vp)) {
        vp <- object$warm_fit$vp
    }
    if (is.null(vp)) {
        if (is.list(object$beta)) vp <- rep(1.0, nrow(object$beta[[1]]))
        else vp <- rep(1.0, nrow(object$beta))
    }
    if (is.list(object$beta)) {
        nlambda <- length(lambda)
        penalty <- numeric(nlambda)
        for (k in seq_len(nlambda)) {
            bet <- unlist(lapply(object$beta, function(b) b[, k]))
            penalty[k] <- pen_function(bet, alpha, rep(vp, length.out = length(bet)))
        }
    } else {
        penalty <- apply(object$beta, 2, pen_function, alpha = alpha, vp = vp)
    }
    dev_val + lambda * penalty
}
