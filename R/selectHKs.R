###############################################################################
## selection of housekeepers (HKs)
###############################################################################
setMethod("selectHKs", signature = "matrix", definition =
            function(qPCRBatch, group, method = "geNorm", minNrHKs = 2, 
                     log = TRUE, Symbols, trace = TRUE, na.rm = TRUE){
              x <- qPCRBatch
              n <- ncol(x)
              if(n < 3)
                stop("you need data from at least 3 variables/columns")
              if(minNrHKs >= n)
                stop("'minNrHKs' must be smaller than 'ncol(x)'")
              if(minNrHKs < 2){
                warning("'minNrHKs' < 2 => 'minNrHKs' is set to 2")
                minNrHKs <- 2
              }
              if(missing(Symbols))
                stop("'Symbols' has to be specified")
              if(length(Symbols) != n)
                stop("'Symbols' has wrong length")
              
              if(method == "geNorm"){
                V <- numeric(n-minNrHKs)
                names(V) <- paste(((n-1):minNrHKs), "/", (n:(minNrHKs+1)), sep = "")
                meanM <- numeric(n-minNrHKs+1)
                names(meanM) <- as.character(n:minNrHKs)
                R <- character(n)
                names(R) <- as.character(c(rep(1, minNrHKs),(minNrHKs+1):length(R)))
                for(i in n:minNrHKs){
                  M <- stabMeasureM(x, log = log, na.rm = na.rm)
                  names(M) <- Symbols
                  ind <- which.max(M)
                  meanM[n-i+1] <- mean(M)
                  if(i == minNrHKs)
                    R[1:minNrHKs] <- Symbols
                  else
                    R[i] <- Symbols[ind]
                  
                  if(i > 2){
                    if(log){
                      NF.old <- rowMeans(x)
                      NF.new <- rowMeans(x[,-ind])
                      V[n-i+1] <- sd(NF.new - NF.old, na.rm = na.rm)
                    }else{
                      NF.old <- apply(x, 1, geomMean, na.rm = na.rm)
                      NF.new <- apply(x[,-ind], 1, geomMean, na.rm = na.rm)
                      V[n-i+1] <- sd(log2(NF.new/NF.old), na.rm = na.rm)
                    }
                  }
                  
                  if(trace){
                    message("###############################################################\n")
                    message("Step ", n-i+1, ":\n")
                    message("stability values M:\n")
                    print(sort(M))
                    message("average stability M:\t", meanM[n-i+1], "\n")
                    if(i > 2){
                      message("variable with lowest stability (largest M value):\t", Symbols[ind], "\n")
                      message("Pairwise variation, (", i-1, "/", i, "):\t", V[n-i+1], "\n")
                    }
                  }
                  x <- x[,-ind]
                  Symbols <- Symbols[-ind]
                }
                return(list(ranking = R, variation = V, meanM = meanM))
              }
              if(method == "NormFinder"){
                
                NF <- stabMeasureRho(x, group = group, log = log, na.rm = na.rm, returnAll = TRUE)
                k <- length(NF$rho)
                R <- character(minNrHKs)
                rho <- NF$rho
                R[1] <- Symbols[which.min(rho)]
                b <- integer(minNrHKs)
                b[1] <- which.min(rho)
                rho.min <- numeric(minNrHKs)
                rho.min[1] <- rho[b[1]]
                
                if(trace){
                  message("###############################################################\n")
                  message("Step ", 1, ":\n")
                  message("stability values rho:\n")
                  print(sort(rho))
                  message("variable with highest stability (smallest rho value):\t", Symbols[b[1]], "\n")
                }
                for(i in 2:minNrHKs){
                  rho[b[i-1]] <- NA
                  
                  for (j in (c(1:k)[-b])){
                    a <- c(b,j)
                    a1 <- NF$d[a,]
                    a2 <- colMeans(a1)*sqrt(k/(k-i))
                    b1 <- NF$v[a,]
                    b2 <- colMeans(b1)/i
                    rho[j] <- mean(abs(a2)+sqrt(b2))
                  }
                  b[i] <- which.min(rho)
                  R[i] <- Symbols[b[i]]
                  rho.min[i] <- rho[b[i]]
                  if(trace){
                    message("###############################################################\n")
                    message("Step ", i, ":\n")
                    message("stability values rho:\n")
                    print(sort(rho[!is.na(rho)]))
                    message("variable with highest stability (smallest rho value):\t", Symbols[b[i]], "\n")
                  }
                }
                names(R) <- 1:minNrHKs
                names(rho.min) <- 1:minNrHKs
                return(list(ranking = R, rho = rho.min))
              }else{
                stop("specified method not yet implemented")
              }
            }
)

setMethod("selectHKs", signature = "qPCRBatch", definition =
  function(qPCRBatch, group, method = "geNorm", minNrHKs = 2, 
 log = TRUE, Symbols, trace = TRUE, na.rm = TRUE){
    x <- t(exprs(qPCRBatch))
    x <- data.matrix(x)
    selectHKs(x, group = group, method = method, minNrHKs = minNrHKs,
              log = log, Symbols = Symbols, trace = trace, na.rm = na.rm)
  }
)
