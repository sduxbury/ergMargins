#' Function to compute edge probabilities from a converged ERGM.
#' this is largely a wrapper function to btergm::edge.prob.The only differnece is built in error
#' handling for curved ergm--the native edge.prob returns an error when used on curved ergms
#' Original authors of edgeprob are Philip Leifeld, imported from btergm package in R
#' @param model An ergm object.
#' @param verbose Determines whether to print updates as the function progresses. Appealing in large networks.
#' edge.prob2()

getformula<-function(x){x@formula}
setMethod("getformula", signature = className("btergm", "btergm"),
          definition = function(x) x@formula)

setMethod("getformula", signature = className("mtergm", "btergm"),
          definition = function(x) x@formula)

setMethod("getformula", signature = className("ergm", "ergm"),
          definition = function(x) x$formula)


edge.prob2<-function (model, verbose = FALSE)
{
  object<-model
  if (class(object) == "ergm") {
    tergm <- FALSE
  }
  else if (class(object) %in% c("btergm", "mtergm")) {
    tergm <- TRUE
  }
  else {
    stop(paste("The edgeprob function is only applicable to ergm, btergm, and",
               "mtergm objects."))
  }

  ##check for duplicated names
  if(class(object)=="ergm"){
    if(any(duplicated(model$network%v%"vertex.names"))){
      stop(paste("Duplicated names detected in network; edgeprob function will fail.",
                 "Please assign unique vertex names before using edgeprob.",
                 "Recommended code:",
                 "model$network%v%'vertex.names'<-1:network.size(model$network)"))
    }
  }


  l <- tergmprepare2(formula = getformula(object), offset = FALSE,
                    blockdiag = FALSE, verbose = FALSE)
  for (cv in 1:length(l$covnames)) {
    assign(l$covnames[cv], l[[l$covnames[cv]]])
  }
  assign("offsmat", l$offsmat)
  form <- stats::as.formula(l$form)
  covnames <- l$covnames[-1]
  coefs <- stats::coef(object)
  if (verbose == TRUE) {
    message("Creating data frame with predictors...")
  }
  Y <- NULL
  dyads <- NULL
  for (i in 1:length(l$networks)) {
    mat <- as.matrix(l$networks[[i]])
    imat <- matrix(rep(1:nrow(mat), ncol(mat)), nrow = nrow(mat))
    if ((class(l$networks[[i]]) == "network" && network::is.bipartite(l$networks[[i]])) ||
        (class(l$networks[[i]]) == "matrix" && xergm.common::is.mat.onemode(l$networks[[i]]) ==
         FALSE)) {
      mn <- nrow(mat) + 1
      mx <- nrow(mat) + ncol(mat)
      jmat <- matrix(rep(mn:mx, nrow(mat)), nrow = nrow(mat),
                     byrow = TRUE)
    }
    else {
      jmat <- matrix(rep(1:ncol(mat), nrow(mat)), nrow = nrow(mat),
                     byrow = TRUE)
    }
    f <- stats::as.formula(paste(l$form, " + edgecov(imat) + edgecov(jmat)"))
    mpli <- ergm::ergmMPLE(f)
    Y <- c(Y, mpli$response)
    dyads <- rbind(dyads, cbind(mpli$predictor, i))
  }
  term.names <- colnames(dyads)[-(length(colnames(dyads)):(length(colnames(dyads)) -
                                                             2))]
  term.names <- c(term.names, "i", "j", "t")
  dyads <- data.frame(dyads)
  colnames(dyads) <- term.names
  dyads <- cbind(Y, dyads)
  colnames(dyads)[1] <- "tie"
  class(dyads[, length(colnames(dyads))]) <- "integer"
  class(dyads[, length(colnames(dyads)) - 1]) <- "integer"
  class(dyads[, length(colnames(dyads)) - 2]) <- "integer"
  cf <- stats::coef(object)
  cf.length <- length(cf)
  cf <- cf[!cf %in% c(Inf, -Inf)]
  if (length(cf) != cf.length) {
    warning(paste("There are structural zeros or ones. For these dyads, the",
                  "predicted probabilities are not valid and must be manually replaced",
                  "by 0 or 1, respectively."))
  }
  cbcoef <- cbind(cf)
  chgstat <- dyads[, 2:(ncol(dyads) - 3)]
  ##handle decay term in curved ergms
  if(ergm::is.curved(object)){
    curved.term<-vector(length=length(object$etamap$curved))
    for(i in 1:length(object$etamap$curved)){
    curved.term[i]<-object$etamap$curved[[i]]$from[2]
    }
    cbcoef<-cbcoef[-c(curved.term)]
  }

  lp <- apply(chgstat, 1, function(x) t(x) %*% cbcoef)
  result <- c(1/(1 + exp(-lp)))
  i.name <- numeric(nrow(dyads))
  j.name <- numeric(nrow(dyads))
  for (t in 1:length(l$networks)) {
    vnames.t <- colnames(l$networks[[t]][, ])
    dyads.t <- dyads[which(dyads$t == t), ]
    i.name.t <- vnames.t[dyads.t$i]
    j.name.t <- vnames.t[dyads.t$j]
    i.name[which(dyads$t == t)] <- i.name.t
    j.name[which(dyads$t == t)] <- j.name.t
  }
  dyads$i.name <- i.name
  dyads$j.name <- j.name
  dyads <- cbind(dyads, result)
  colnames(dyads)[ncol(dyads)] <- "probability"
  dyads <- dyads[order(dyads$t, dyads$i, dyads$j), ]
  rownames(dyads) <- NULL
  return(dyads)
}
