#'this is an internal workhorse function imported from btergm
#'Original authors are Philip Leifeld and Bruce Desmairis
#'Not intended for external use
#'
tergmprepare2<-function (formula, offset = TRUE, blockdiag = FALSE, verbose = TRUE)
{
  l <- list()
  l$lhs.original <- deparse(formula[[2]])
  l$networks <- eval(parse(text = deparse(formula[[2]])), envir = environment(formula))
  if (class(l$networks) %in% "list" || class(l$networks) %in%
      "network.list") {
  } else {
    l$networks <- list(l$networks)
  }
  if("mlnet"%in%class(l$networks[[1]])){
    class(l$networks[[1]])<-"network"
  }

  for (i in 1:length(l$networks)) {
    if (!class(l$networks[[i]]) %in% c("network", "matrix",
                                       "list")) {
      tryCatch({
        l$networks[[i]] <- as.matrix(l$networks[[i]])
      }, error = function(cond) {
        stop(paste("Object", i, "could not be converted to a matrix."))
      })
    }
  }
  l$num.vertices <- max(sapply(l$networks, function(x) network::get.network.attribute(network::network(x),
                                                                                      "n")))
  if (network::is.network(l$networks[[1]])) {
    l$directed <- network::is.directed(l$networks[[1]])
    l$bipartite <- network::is.bipartite(l$networks[[1]])
  } else {
    if (is.mat.directed(as.matrix(l$networks[[1]]))) {
      l$directed <- TRUE
    }else {
      l$directed <- FALSE
    }
    if (is.mat.onemode(as.matrix(l$networks[[1]]))) {
      l$bipartite <- FALSE
    } else {
      l$bipartite <- TRUE
    }
  }
  l$form <- stats::update.formula(formula, networks[[i]] ~ .)
  l$time.steps <- length(l$networks)
  tilde <- deparse(l$form[[1]])
  lhs <- deparse(l$form[[2]])
  rhs <- paste(deparse(l$form[[3]]), collapse = "")
  rhs <- gsub("\\s+", " ", rhs)
  rhsterms <- strsplit(rhs, "\\s*\\+\\s*")[[1]]
  if (length(rhsterms) > 1) {
    for (i in length(rhsterms):2) {
      leftbracketmatches <- gregexpr("\\(", rhsterms[i])[[1]]
      leftbracketmatches <- leftbracketmatches[leftbracketmatches !=
                                                 -1]
      leftbracketmatches <- length(leftbracketmatches)
      rightbracketmatches <- gregexpr("\\)", rhsterms[i])[[1]]
      rightbracketmatches <- rightbracketmatches[rightbracketmatches !=
                                                   -1]
      rightbracketmatches <- length(rightbracketmatches)
      if (leftbracketmatches != rightbracketmatches) {
        rhsterms[i - 1] <- paste(rhsterms[i - 1], rhsterms[i],
                                 sep = " + ")
        rhsterms <- rhsterms[-i]
      }
    }
  }
  l$rhs.terms <- rhsterms
  rhs.operators <- rep("+", length(l$rhs.terms) - 1)
  covnames <- character()
  for (k in 1:length(l$rhs.terms)) {
    if (grepl("((edge)|(dyad))cov", l$rhs.terms[k])) {
      if (grepl(",\\s*?((attr)|\\\")", l$rhs.terms[k])) {
        s <- "((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,\\s*a*.*?)\\)(?:\\))?)"
      }else {
        s <- "((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,*\\s*a*.*?)\\)(?:\\))?)"
      }
      x1 <- sub(s, "\\1", l$rhs.terms[k], perl = TRUE)
      x2 <- sub(s, "\\5", l$rhs.terms[k], perl = TRUE)
      if (grepl("\\[.*\\]", x2)) {
        stop(paste0("Covariate names are not allowed to have indices: ",
                    x2, ". Please prepare a list object before estimation."))
      }
      if (grepl("^\"", x2))
        next
      x3 <- sub(s, "\\6", l$rhs.terms[k], perl = TRUE)
      x.current <- eval(parse(text = x2), envir = environment(formula))
      type <- class(x.current)
      l$covnames <- c(l$covnames, x2)
      l[[x2]] <- x.current
      if (grepl("\\[i\\]+$", x2)) {
        stop(paste0("Error in the following model term: ",
                    l$rhs.terms[k], ". The index 'i' is used internally by btergm. Please use a ",
                    "different index, for example 'j'."))
      }
      if (grepl("[^\\]]\\]$", x2)) {
        l$rhs.terms[k] <- paste0(x1, x2, x3)
        if (type[1] %in% c("matrix", "network",
                        "dgCMatrix", "dgTMatrix", "dsCMatrix",
                        "dsTMatrix", "dgeMatrix")) {
          x.current <- list(x.current)
          l[[x2]] <- x.current
        }
        if (length(x.current) != l$time.steps) {
          stop(paste(x2, "has", length(x.current),
                     "elements, but there are", l$time.steps,
                     "networks to be modeled."))
        }
        if (blockdiag == TRUE) {
        }else {
          x2 <- paste0(x2, "[[i]]")
        }
      }else if (type[1] %in% c("matrix", "network",
                           "dgCMatrix", "dgTMatrix", "dsCMatrix",
                           "dsTMatrix", "dgeMatrix")) {
        if (!type[1] %in% c("matrix", "network")) {
          x.current <- as.matrix(x.current)
        }
        l[[x2]] <- list()
        for (i in 1:l$time.steps) {
          l[[x2]][[i]] <- x.current
        }
        if (blockdiag == TRUE) {
        }else {
          x2 <- paste0(x2, "[[i]]")
        }
        l$rhs.terms[k] <- paste(x1, x2, x3, sep = "")
      }else if (type == "list" || type == "network.list") {
        if (length(x.current) != l$time.steps) {
          stop(paste(x2, "has", length(get(x2)),
                     "elements, but there are", l$time.steps,
                     "networks to be modeled."))
        }
        if (blockdiag == TRUE) {
        }else {
          x2 <- paste0(x2, "[[i]]")
        }
        l$rhs.terms[k] <- paste0(x1, x2, x3)
      }else {
        tryCatch({
          l[[x2]] <- list(rep(as.matrix(x.current)),
                          l$time.steps)
        }, error = function(cond) {
          stop(paste0("Object '", x2, "' could not be converted to a matrix."))
        })
      }
    }else if (grepl("memory", l$rhs.terms[k])) {
      s <- "(?:memory\\((?:.*type\\s*=\\s*)?(?:\"|'))(\\w+)(?:(\"|').*\\))"
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        type <- "stability"
      }else {
        type <- sub(s, "\\1", l$rhs.terms[k], perl = TRUE)
      }
      s <- "(?:memory\\(.*lag\\s*=\\s*)(\\d+)(?:.*\\))"
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        lag <- 1
      }else {
        lag <- as.integer(sub(s, "\\1", l$rhs.terms[k],
                              perl = TRUE))
      }
      if (lag > length(l$networks) - 1) {
        stop("The 'lag' argument in the 'memory' term is too large.")
      }
      mem <- l$networks[-(length(l$networks):(length(l$networks) -
                                                lag + 1))]
      mem <- lapply(mem, as.matrix)
      memory <- list()
      for (i in 1:length(mem)) {
        if (type == "autoregression") {
          memory[[i]] <- mem[[i]]
        }else if (type == "stability") {
          mem[[i]][mem[[i]] == 0] <- -1
          memory[[i]] <- mem[[i]]
        }else if (type == "innovation") {
          memory[[i]] <- mem[[i]]
          memory[[i]][mem[[i]] == 0] <- 1
          memory[[i]][mem[[i]] == 1] <- 0
        }else if (type == "loss") {
          memory[[i]] <- mem[[i]]
          memory[[i]][mem[[i]] == 0] <- 0
          memory[[i]][mem[[i]] == 1] <- -1
        }else {
          stop("'type' argument in the 'memory' term not recognized.")
        }
      }
      rm(mem)
      l[["memory"]] <- memory
      if (blockdiag == TRUE) {
        l$rhs.terms[k] <- "edgecov(memory)"
      }else {
        l$rhs.terms[k] <- "edgecov(memory[[i]])"
      }
      l$covnames <- c(l$covnames, "memory")
    }else if (grepl("delrecip", l$rhs.terms[k])) {
      s <- "(?:delrecip\\((?:.*mutuality\\s*=\\s*)?)((TRUE)|(FALSE)|T|F)(?:.*\\))"
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        mutuality <- FALSE
      }else {
        mutuality <- as.logical(sub(s, "\\1", l$rhs.terms[k],
                                    perl = TRUE))
      }
      s <- "(?:delrecip\\(.*lag\\s*=\\s*)(\\d+)(?:.*\\))"
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        lag <- 1
      }else {
        lag <- as.integer(sub(s, "\\1", l$rhs.terms[k],
                              perl = TRUE))
      }
      if (lag > length(l$networks) - 1) {
        stop("The 'lag' argument in the 'delrecip' term is too large.")
      }
      dlr <- l$networks[-(length(l$networks):(length(l$networks) -
                                                lag + 1))]
      dlr <- lapply(dlr, function(x) t(as.matrix(x)))
      delrecip <- list()
      for (i in 1:length(dlr)) {
        delrecip[[i]] <- dlr[[i]]
        if (mutuality == TRUE) {
          delrecip[[i]][dlr[[i]] == 0] <- -1
        }
      }
      rm(dlr)
      l[["delrecip"]] <- delrecip
      if (blockdiag == TRUE) {
        l$rhs.terms[k] <- "edgecov(delrecip)"
      }else {
        l$rhs.terms[k] <- "edgecov(delrecip[[i]])"
      }
      l$covnames <- c(l$covnames, "delrecip")
    }else if (grepl("timecov", l$rhs.terms[k])) {
      s <- "(?:timecov\\((?:.*x\\s*=\\s*)?)(\\w+)(?:.*\\))"
      if (sub(s, "\\1", l$rhs.terms[k], perl = TRUE) %in%
          c("minimum", "maximum", "transform",
            "min", "max", "trans")) {
        s <- "(?:timecov\\(?:.*x\\s*=\\s*)(\\w+)(?:.*\\))"
      }
      countprevtc <- 1
      if (k > 1) {
        for (i in (k - 1):1) {
          if (grepl("timecov", l$rhs.terms[i])) {
            countprevtc <- countprevtc + 1
          }
        }
      }
      if (countprevtc > 0) {
        countprevtc <- as.character(countprevtc)
      }else {
        countprevtc <- ""
      }
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        x <- NULL
        label <- paste0("timecov", countprevtc)
      }else {
        x <- sub(s, "\\1", l$rhs.terms[k], perl = TRUE)
        label <- paste0("timecov", countprevtc,
                        ".", x)
      }
      s <- "(?:timecov\\(.*minimum\\s*=\\s*)(\\d+)(?:.*\\))"
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        minimum <- 1
      }else {
        minimum <- as.integer(sub(s, "\\1", l$rhs.terms[k],
                                  perl = TRUE))
      }
      s <- "(?:timecov\\(.*maximum\\s*=\\s*)(\\d+)(?:.*\\))"
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        maximum <- l$time.steps
      }else {
        maximum <- as.integer(sub(s, "\\1", l$rhs.terms[k],
                                  perl = TRUE))
      }
      s <- "(?:timecov\\(.*transform\\s*=\\s*)(.+?)(?:(?:,|\\)$)]*.*)"
      if (grepl(s, l$rhs.terms[k]) == FALSE) {
        transform <- function(t) t
      }else {
        transform <- eval(parse(text = sub(s, "\\1",
                                           l$rhs.terms[k], perl = TRUE)))
      }
      if (is.null(x)) {
        covariate <- l[["networks"]]
        onlytime <- TRUE
      }else {
        onlytime <- FALSE
        covariate <- get(x)
      }
      tc <- timecov2(covariate = covariate, minimum = minimum,
                    maximum = maximum, transform = transform, onlytime = onlytime)
      l[[label]] <- tc
      labelsuffix <- sub(s, "\\1", l$rhs.terms[k],
                         perl = TRUE)
      labelsuffix <- if (blockdiag == TRUE) {
        l$rhs.terms[k] <- paste0("edgecov(", label,
                                 ")")
      }else {
        l$rhs.terms[k] <- paste0("edgecov(", label,
                                 "[[i]])")
      }
      l$covnames <- c(l$covnames, label)
      if (verbose == TRUE) {
        timecovreporting <- matrix(sapply(tc, function(x) mean(x[1,
                                                                 2])), nrow = 1)
        colnames(timecovreporting) <- paste0("t=",
                                             1:length(l$networks))
        rownames(timecovreporting) <- ""
        message("Mean transformed timecov values:")
        print(timecovreporting)
      }
    }
  }
  l$covnames <- c("networks", l$covnames)
  lengths <- sapply(l$covnames, function(cn) length(l[[cn]]))
  mn <- max(lengths)
  if (length(table(lengths)) > 1) {
    mn <- min(lengths)
    l$time.steps <- mn
    for (i in 1:length(l$covnames)) {
      cn <- l$covnames[[i]]
      ll <- l[[cn]]
      difference <- length(ll) - mn
      if (difference > 0) {
        l[[cn]] <- ll[(difference + 1):length(ll)]
      }
    }
  }

  t.end <- max(lengths)
  t.start <- t.end - mn + 1
  if (verbose == TRUE) {
    if (length(l$covnames) > 1) {
      dimensions <- lapply(lapply(l$covnames, function(x) l[[x]]),
                           function(y) sapply(y, function(z) dim(as.matrix(z))))
      rownames(dimensions[[1]]) <- paste(l$lhs.original,
                                         c("(row)", "(col)"))
      for (i in 2:length(dimensions)) {
        rownames(dimensions[[i]]) <- c(paste(l$covnames[i],
                                             "(row)"), paste(l$covnames[i], "(col)"))
      }
      dimensions <- do.call(rbind, dimensions)
      colnames(dimensions) <- paste0("t=", t.start:t.end)
      message("\nInitial dimensions of the network and covariates:")
      print(dimensions)
    }else {
      message("\nNo covariates provided.")
    }
  }

  l$auto.adjust <- FALSE
  if (length(l$covnames) > 1) {
    nr <- lapply(lapply(l$covnames, function(x) l[[x]]),
                 function(y) sapply(y, function(z) nrow(as.matrix(z))))
    nr <- do.call(rbind, nr)
    nc <- lapply(lapply(l$covnames, function(x) l[[x]]),
                 function(y) sapply(y, function(z) ncol(as.matrix(z))))
    nc <- do.call(rbind, nc)
    for (i in 1:ncol(nr)) {
      if (length(unique(nr[, i])) > 1) {
        l$auto.adjust <- TRUE
      }
    }
    for (i in 1:ncol(nc)) {
      if (length(unique(nc[, i])) > 1) {
        l$auto.adjust <- TRUE
      }
    }
    if (verbose == TRUE && l$auto.adjust == TRUE) {
      message(paste("\nDimensions differ across networks within time steps."))
    }
    if (l$auto.adjust == TRUE) {
      for (i in 1:length(l$covnames)) {
        for (t in 1:l$time.steps) {
          if (is.null(rownames(as.matrix(l[[l$covnames[i]]][[t]]))) ||
              is.null(colnames(as.matrix(l[[l$covnames[i]]][[t]])))) {
            stop(paste0("The dimensions of the covariates differ, but ",
                        "covariate '", l$covnames[i], " does not have node labels at t = ",
                        t, ". Automatic adjustment of dimensions is therefore not ",
                        "possible."))
          }
        }
      }
    }
    if (l$auto.adjust == FALSE) {
      for (t in 1:l$time.steps) {
        rlabels.i <- list()
        clabels.i <- list()
        for (i in 1:length(l$covnames)) {
          rlabels.i[[i]] <- rownames(as.matrix(l[[l$covnames[i]]][[t]]))
          clabels.i[[i]] <- colnames(as.matrix(l[[l$covnames[i]]][[t]]))
        }
        rlabels.i <- do.call(rbind, rlabels.i)
        clabels.i <- do.call(rbind, clabels.i)
        flag <- FALSE
        if (!is.null(rlabels.i)) {
          for (j in 1:ncol(rlabels.i)) {
            if (length(unique(rlabels.i[, j])) > 1) {
              l$auto.adjust <- TRUE
              flag <- TRUE
              break
            }
          }
        }
        if (!is.null(clabels.i)) {
          for (j in 1:ncol(clabels.i)) {
            if (length(unique(clabels.i[, j])) > 1) {
              l$auto.adjust <- TRUE
              flag <- TRUE
              break
            }
          }
        }
      }
      if (verbose == TRUE && flag == TRUE) {
        message(paste("\nSame dimensions but different labels across",
                      "networks within time steps."))
      }
    }
  }
  if (verbose == TRUE && l$auto.adjust == TRUE) {
    message("Trying to auto-adjust the dimensions of the networks. ",
            "If this fails, provide conformable matrices or network objects.")
  }else if (verbose == TRUE) {
    message("\nAll networks are conformable.")
  }
  structzero.df <- data.frame(label = character(), time = integer(),
                              object = character(), where = character())




  ###error below

  if (length(l$covnames) > 0 && l$auto.adjust == TRUE) {
    for (i in 1:l$time.steps) {
      for (j in 1:length(l$covnames)) {
        for (k in 1:length(l$covnames)) {
          if (j != k) {
            nw.j <- l[[l$covnames[j]]][[i]]
            rn.j <- rownames(as.matrix(nw.j))
            cn.j <- colnames(as.matrix(nw.j))
            nr.j <- nrow(as.matrix(nw.j))
            nc.j <- ncol(as.matrix(nw.j))
            nw.k <- l[[l$covnames[k]]][[i]]
            nr.k <- nrow(as.matrix(nw.k))
            nc.k <- ncol(as.matrix(nw.k))
            cn.k <- colnames(as.matrix(nw.k))
            if(nr.k==nc.k){
              rn.k<-cn.k
            }else{
            rn.k <- rownames(as.matrix(nw.k))
            }

            if (is.null(rn.j) || is.null(cn.j)) {
              stop(paste0("Missing row or column labels in object '",
                          l$covnames[j], "'. Provide row and column ",
                          "labels for all networks and covariates."))
            }
            else if (is.null(rn.k) || is.null(cn.k)) {
              stop(paste0("Missing row or column labels in object '",
                          l$covnames[k], "'. Provide row and column ",
                          "labels for all networks and covariates."))
            }
            else {
              if (is.null(rn.j) && !is.null(rn.k) &&
                  nr.j == nr.k) {
                if (class(nw.j) %in% "network") {
                  network::set.vertex.attribute(nw.j,
                                                "vertex.names", rn.k)
                }
                else {
                  rownames(nw.j) <- rn.k
                }
              }
              else if (is.null(rn.k) && !is.null(rn.j) &&
                       nr.j == nr.k) {
                if (class(nw.k) %in% "network") {
                  network::set.vertex.attribute(nw.k,
                                                "vertex.names", rn.j)
                }
                else {
                  rownames(nw.k) <- rn.j
                }
              }
              else if ((is.null(rn.k) || is.null(rn.j)) &&
                       nr.j != nr.k) {
                stop(paste0("Object '", l$covnames[j],
                            "' is incompatible with object '",
                            l$covnames[k], "' at t = ", i,
                            "."))
              }
              nw.j.labels <- btergm::adjust(nw.j, nw.k, remove = FALSE,
                                    value = 1, returnlabels = TRUE)
              nw.j <- btergm::adjust(nw.j, nw.k, remove = FALSE,
                             value = 1)
              l[[l$covnames[j]]][[i]] <- nw.j
              ro <- nw.j.labels$added.row
              co <- nw.j.labels$added.col
              if (length(ro) > 0) {
                ro <- data.frame(label = ro, time = rep(i,
                                                        length(ro)), object = rep(l$covnames[j],
                                                                                  length(ro)), where = rep("row",
                                                                                                           length(ro)))
                structzero.df <- rbind(structzero.df,
                                       ro)
              }
              if (length(co) > 0) {
                co <- data.frame(label = co, time = rep(i,
                                                        length(co)), object = rep(l$covnames[j],
                                                                                  length(co)), where = rep("col",
                                                                                                           length(co)))
                structzero.df <- rbind(structzero.df,
                                       co)
              }
              nw.k.labels <- btergm::adjust(nw.k, nw.j, remove = FALSE,
                                    value = 1, returnlabels = TRUE)
              nw.k <- btergm::adjust(nw.k, nw.j, remove = FALSE,
                             value = 1)
              l[[l$covnames[k]]][[i]] <- nw.k
              ro <- nw.k.labels$added.row
              co <- nw.k.labels$added.col
              if (length(ro) > 0) {
                ro <- data.frame(label = ro, time = rep(i,
                                                        length(ro)), object = rep(l$covnames[j],
                                                                                  length(ro)), where = rep("row",
                                                                                                           length(ro)))
                structzero.df <- rbind(structzero.df,
                                       ro)
              }
              if (length(co) > 0) {
                co <- data.frame(label = co, time = rep(i,
                                                        length(co)), object = rep(l$covnames[j],
                                                                                  length(co)), where = rep("col",
                                                                                                           length(co)))
                structzero.df <- rbind(structzero.df,
                                       co)
              }
            }
          }
        }
      }
    }
  }


  nr.net <- sapply(l$networks, function(x) nrow(as.matrix(x)))
  for (i in 1:length(l$covnames)) {
    nr <- sapply(l[[l$covnames[i]]], function(x) {
      nrow(as.matrix(x))
    })
    for (j in 1:l$time.steps) {
      if (nr[j] != nr.net[j]) {
        stop(paste0("Covariate object '", l$covnames[i],
                    "' does not have the same number of rows as the dependent ",
                    "network at time step ", j, "."))
      }
    }
  }
  nc.net <- sapply(l$networks, function(x) ncol(as.matrix(x)))
  for (i in 1:length(l$covnames)) {
    nc <- sapply(l[[l$covnames[i]]], function(x) {
      ncol(as.matrix(x))
    })
    for (j in 1:l$time.steps) {
      if (nc[j] != nc.net[j]) {
        stop(paste0("Covariate object '", l$covnames[i],
                    "' does not have the same number of columns as the dependent ",
                    "network at time step ", j, "."))
      }
    }
  }
  if (verbose == TRUE) {
    if (l$auto.adjust == TRUE) {
      sz.row <- unique(structzero.df[structzero.df$where ==
                                       "row", -3])
      szrownum <- numeric(length(l$networks))
      for (i in 1:length(l$networks)) {
        szrownum[i] <- nrow(sz.row[sz.row$time == i,
                                   ])
      }
      sz.col <- unique(structzero.df[structzero.df$where ==
                                       "col", -3])
      szcolnum <- numeric(length(l$networks))
      for (i in 1:length(l$networks)) {
        szcolnum[i] <- nrow(sz.col[sz.col$time == i,
                                   ])
      }
      totrow <- sapply(l$networks, function(x) nrow(as.matrix(x)))
      totcol <- sapply(l$networks, function(x) ncol(as.matrix(x)))
      if (offset == TRUE) {
        dimensions <- rbind(totrow, totcol, szrownum,
                            szcolnum, totrow - szrownum, totcol - szcolnum)
        rownames(dimensions) <- c("total number of rows",
                                  "total number of columns", "row-wise structural zeros",
                                  "column-wise structural zeros", "remaining rows",
                                  "remaining columns")
      }
      else {
        dimensions <- rbind(szrownum, szcolnum, totrow -
                              szrownum, totcol - szcolnum)
        rownames(dimensions) <- c("maximum deleted nodes (row)",
                                  "maximum deleted nodes (col)", "remaining rows",
                                  "remaining columns")
      }
      colnames(dimensions) <- paste0("t=", t.start:t.end)
      if (nrow(structzero.df) > 0) {
        if (offset == TRUE) {
          message("\nNodes affected completely by structural zeros:")
        }
        else {
          message("\nAbsent nodes:")
        }
        szcopy <- structzero.df
        szcopy$time <- szcopy$time - 1 + t.start
        print(unique(szcopy))
      }
      else {
        message("\nAll nodes are retained.")
      }
      message("\nNumber of nodes per time step after adjustment:")
      print(dimensions)
    }
  }
  l$nvertices <- sapply(l$networks, function(x) c(nrow(as.matrix(x)),
                                                  ncol(as.matrix(x))))
  rownames(l$nvertices) <- c("row", "col")
  colnames(l$nvertices) <- paste0("t=", t.start:t.end)
  l$offsmat <- list()
  for (i in 1:l$time.steps) {
    mat <- matrix(0, nrow = nrow(as.matrix(l$networks[[i]])),
                  ncol = ncol(as.matrix(l$networks[[i]])))
    rownames(mat) <- rownames(as.matrix(l$networks[[i]]))
    colnames(mat) <- colnames(as.matrix(l$networks[[i]]))
    l$offsmat[[i]] <- mat
  }
  if (nrow(structzero.df) > 0) {
    for (i in 1:nrow(structzero.df)) {
      if (structzero.df$where[i] == "row") {
        index <- which(rownames(l$offsmat[[structzero.df$time[i]]]) ==
                         structzero.df$label[i])
        l$offsmat[[structzero.df$time[i]]][index, ] <- 1
      }
      else {
        index <- which(colnames(l$offsmat[[structzero.df$time[i]]]) ==
                         structzero.df$label[i])
        l$offsmat[[structzero.df$time[i]]][, index] <- 1
      }
    }
  }
  if (offset == TRUE) {
    l$rhs.terms[length(l$rhs.terms) + 1] <- "offset(edgecov(offsmat[[i]]))"
    rhs.operators[length(rhs.operators) + 1] <- "+"
  }  else {
    if (l$auto.adjust == TRUE) {
      l$offsmat <- suppressMessages(btergm::handleMissings(l$offsmat,
                                                   na = 1, method = "remove"))
      for (j in 1:length(l$covnames)) {
        l[[l$covnames[j]]] <- btergm::adjust(l[[l$covnames[j]]],
                                     l$offsmat)
      }
    }
  }
  if (verbose == TRUE && length(l$covnames) > 1) {
    dimensions <- lapply(lapply(l$covnames, function(x) l[[x]]),
                         function(y) sapply(y, function(z) dim(as.matrix(z))))
    rownames(dimensions[[1]]) <- paste(l$lhs.original, c("(row)",
                                                         "(col)"))
    for (i in 2:length(dimensions)) {
      rownames(dimensions[[i]]) <- c(paste(l$covnames[i],
                                           "(row)"), paste(l$covnames[i], "(col)"))
    }
    dimensions <- do.call(rbind, dimensions)
    colnames(dimensions) <- paste0("t=", t.start:t.end)
    message("\nDimensions of the network and covariates after adjustment:")
    print(dimensions)
  }
  rhs <- l$rhs.terms[1]
  if (length(rhs.operators) > 0) {
    for (i in 1:length(rhs.operators)) {
      rhs <- paste(rhs, rhs.operators[i], l$rhs.terms[i +
                                                        1])
    }
  }
  f <- paste(lhs, tilde, rhs)
  l$form <- stats::as.formula(f, env = environment())
  if (blockdiag == TRUE) {
    if (l$bipartite == TRUE) {
      stop(paste("MCMC estimation is currently only supported for one-mode",
                 "networks. Use the btergm function instead."))
    }
    l$form <- stats::update.formula(l$form, networks ~ .)
    l$form <- paste(deparse(l$form), collapse = "")
    l$form <- paste(l$form, "+ offset(edgecov(offsmat))")
    l$form <- stats::as.formula(l$form, env = environment())
    if (length(l$covnames) > 1) {
      for (j in 2:length(l$covnames)) {
        l[[l$covnames[j]]] <- as.matrix(Matrix::bdiag(lapply(l[[l$covnames[j]]],
                                                             as.matrix)))
      }
    }
    l$offsmat <- as.matrix(Matrix::bdiag(l$offsmat))
    bdoffset <- lapply(l$networks, as.matrix)
    for (i in 1:length(bdoffset)) {
      bdoffset[[i]][, ] <- 1
    }
    bdoffset <- as.matrix((Matrix::bdiag(bdoffset) - 1) *
                            -1)
    l$offsmat <- l$offsmat + bdoffset
    rm(bdoffset)
    l$offsmat[l$offsmat > 0] <- 1
    if (class(l$networks[[1]]) %in% "network") {
      attrnames <- network::list.vertex.attributes(l$networks[[1]])
      attributes <- list()
      for (i in 1:length(l$networks)) {
        attrib <- list()
        for (j in 1:length(attrnames)) {
          attrib[[j]] <- network::get.vertex.attribute(l$networks[[i]],
                                                       attrnames[j])
        }
        attributes[[i]] <- attrib
        l$networks[[i]] <- as.matrix(l$networks[[i]])
      }
      l$networks <- network::network(as.matrix(Matrix::bdiag(l$networks)),
                                     directed = l$directed, bipartite = l$bipartite)
      for (i in 1:length(attrnames)) {
        attrib <- unlist(lapply(attributes, function(x) x[[i]]))
        network::set.vertex.attribute(l$networks, attrnames[i],
                                      attrib)
      }
    }
    else {
      l$networks <- network::network(as.matrix(Matrix::bdiag(l$networks)),
                                     directed = l$directed, bipartite = l$bipartite)
    }
    if (verbose == TRUE) {
      cat("\n")
    }
  }
  form3 <- paste(deparse(l$form[[3]]), collapse = "")
  form3 <- gsub("\\s+", " ", form3)
  l$form <- paste(deparse(l$form[[2]]), deparse(l$form[[1]]),
                  form3)
  return(l)
}



#' Check if a matrix is a one-mode matrix
#'
#' Check if a matrix is a one-mode matrix.
#'
#' @param mat A matrix object containing zeros and ones.
#' @return \code{TRUE} if the input matrix \code{mat} represents a one-mode
#'   network and \code{FALSE} otherwise.
#'
#' @noRd
is.mat.onemode <- function(mat) {
  if (nrow(mat) != ncol(mat)) {
    return(FALSE)
  } else if (!is.null(rownames(mat)) && !is.null(colnames(mat))
             && any(rownames(mat) != colnames(mat))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' Check if a matrix represents a directed network
#'
#' Check if a matrix represents a directed network.
#'
#' @param mat A matrix object containing zeros and ones.
#' @return \code{TRUE} if the input matrix \code{mat} represents a directed
#'   network and \code{FALSE} otherwise.
#'
#' @noRd
is.mat.directed <- function(mat) {
  if (nrow(mat) != ncol(mat)) {
    return(FALSE)
  } else if (!is.null(rownames(mat)) && !is.null(colnames(mat))
             && any(rownames(mat) != colnames(mat), na.rm = TRUE)) {
    return(FALSE)
  } else {
    if (any(as.matrix(mat) != t(as.matrix(mat)), na.rm = TRUE)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}

