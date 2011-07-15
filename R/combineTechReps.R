setMethod("combineTechReps", signature = "qPCRBatch", definition =
  function(qPCRBatch, calc = "arith") {
    expM <- exprs(qPCRBatch)
    origDetectors <- row.names(expM)
    if (FALSE %in% grepl("_TechReps", origDetectors)) stop("These are not tech reps")
    newDetectors <- unique(gsub("_TechReps.\\d","", origDetectors))
    NewExpM <- matrix(nrow = length(newDetectors), ncol = dim(expM)[2], dimnames = list(newDetectors,colnames(expM)))
    if (calc == "arith") {
      for (detector in newDetectors) {
        dValues <- colMeans(expM[gsub("_TechReps.\\d", "", origDetectors) %in% detector, ], na.rm = TRUE)
        NewExpM[detector, ] <- dValues
      }
    }
    else if (calc == "geom") {
      for (detector in newDetectors) {
        dValues <- apply(expM[gsub("_TechReps.\\d", "", origDetectors) %in% detector, ], 2, geomMean, na.rm = TRUE)
        NewExpM[detector, ] <- dValues
      }
    }
    else if (calc == "median") {
      for (detector in newDetectors) {
        dValues <- apply(expM[gsub("_TechReps.\\d", "", origDetectors) %in% detector, ], 2, median, na.rm = TRUE)
        NewExpM[detector, ] <- dValues
      }
    }
    else {
      stop("ensure you have specified the correct centrality measure, 'arith', 'geom' or 'median'")
    }
    NewExpM[is.na(NewExpM)] <- NA # make NAs real NAs
    qPCRBatch <- new("qPCRBatch", exprs = NewExpM, phenoData = phenoData(qPCRBatch))
    return(qPCRBatch)
  }
)
