setMethod("combineTechRepsWithSD", signature = "qPCRBatch", definition =
  function(qPCRBatch, calc = "arith") {
    expM <- exprs(qPCRBatch)
    origDetectors <- row.names(expM)
    if (FALSE %in% grepl("_TechReps", origDetectors)) stop("These are not tech reps")
    newDetectors <- unique(gsub("_TechReps.\\d","", origDetectors))
    NewExpM <- matrix(nrow = length(newDetectors), ncol = dim(expM)[2], dimnames = list(newDetectors,colnames(expM)))
    NewSeExpM <- matrix(nrow = length(newDetectors), ncol = dim(expM)[2], dimnames = list(newDetectors,colnames(expM)))
    if (calc == "arith") {
      for (detector in newDetectors) {
        tmp <- expM[gsub("_TechReps.\\d", "", origDetectors) %in% detector, ]
        dValues <- colMeans(tmp, na.rm = TRUE)
        NewExpM[detector, ] <- dValues
        sdValues <- apply(tmp, 2, sd, na.rm = TRUE)
        NewSeExpM[detector, ] <- sdValues
      }
    }
    else if (calc == "geom") {
      for (detector in newDetectors) {
        tmp <- expM[gsub("_TechReps.\\d", "", origDetectors) %in% detector, ]
        dValues <- apply(tmp, 2, geomMean, na.rm = TRUE)
        NewExpM[detector, ] <- dValues
        sdValues <- exp(apply(log(tmp), 2, sd, na.rm = TRUE))
        NewSeExpM[detector, ] <- sdValues        
      }
    }
    else if (calc == "median") {
      for (detector in newDetectors) {
        tmp <- expM[gsub("_TechReps.\\d", "", origDetectors) %in% detector, ]
        dValues <- apply(tmp, 2, median, na.rm = TRUE)
        NewExpM[detector, ] <- dValues
        sdValues <- apply(tmp, 2, mad, na.rm = TRUE)
        NewSeExpM[detector, ] <- sdValues
      }
    }
    else {
      stop("ensure you have specified the correct centrality measure, 'arith', 'geom' or 'median'")
    }
    NewExpM[is.na(NewExpM)] <- NA # make NAs real NAs
    qPCRBatch <- new("qPCRBatch", exprs = NewExpM, se.exprs = NewSeExpM, phenoData = phenoData(qPCRBatch))
    return(qPCRBatch)
  }
)
