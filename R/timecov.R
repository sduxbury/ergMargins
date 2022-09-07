
#'this is an internal helper function imported from btergm
#'Original authors are Philip Leifeld and Bruce Desmairis
#'Not intended for external use
#'
timecov2<-function (covariate, minimum = 1, maximum = length(covariate),
          transform = function(t) 1 + (0 * t) + (0 * t^2), onlytime = FALSE)
{
  if (!class(covariate)%in%"list") {
    stop("'covariate' must be a list of matrices or network objects.")
  }
  for (i in 1:length(covariate)) {
    if (network::is.network(covariate[[i]])) {
      covariate[[i]] <- as.matrix(covariate[[i]])
    }
    else if (!is.matrix(covariate[[i]])) {
      stop("'covariate' must be a list of matrices or network objects.")
    }
  }
  if (is.null(minimum) || is.null(maximum) || !is.numeric(minimum) ||
      !is.numeric(maximum) || length(minimum) > 1 || length(maximum) >
      1) {
    stop("'minimum' and 'maximum' must be single numeric values.")
  }
  if (is.null(transform)) {
    transform <- function(t) 1 + (0 * t) + (0 * t^2)
  }
  else if (!is.function(transform)) {
    stop("'transform' must be a function.")
  }
  l <- 1:length(covariate)
  values <- transform(l)
  if (is.null(values) || any(is.null(values)) || any(!is.finite(values)) ||
      any(is.na(values)) || any(!is.numeric(values))) {
    stop("The 'transform' function produces non-numeric values.")
  }
  values <- values * (l >= minimum) * (l <= maximum)
  timecov <- list()
  for (i in 1:length(l)) {
    if (onlytime == FALSE) {
      timecov[[i]] <- covariate[[i]] * matrix(values[i],
                                              nrow = nrow(covariate[[i]]), ncol = ncol(covariate[[i]]))
    }
    else {
      timecov[[i]] <- matrix(values[i], nrow = nrow(covariate[[i]]),
                             ncol = ncol(covariate[[i]]))
    }
  }
  for (i in 1:length(timecov)) {
    rownames(timecov[[i]]) <- rownames(covariate[[i]])
    colnames(timecov[[i]]) <- colnames(covariate[[i]])
  }
  return(timecov)
}
