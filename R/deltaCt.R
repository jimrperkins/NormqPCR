setMethod("deltaCt", signature = "qPCRBatch", definition =
  function(qPCRBatch, hkgs, combineHkgs=FALSE, calc="arith") {
    hkgs <- make.names(hkgs)
    if(FALSE %in% (hkgs %in% featureNames(qPCRBatch))) stop ("given housekeeping gene, ", hkgs," not found in file. Ensure entered housekeeping genes appear in the file")
    expM <- exprs(qPCRBatch)
    hkgM <- expM[hkgs, ]

    if(length(hkgs) > 1) {
      if (TRUE %in% apply(hkgM, 1, is.na))  {
        warning("NAs present in housekeeping genes readings")
        if (0 %in% apply(! apply(hkgM, 1, is.na),2,sum)) stop("Need at least 1 non NA for each housekeeper")
      }
      hkgV <- vector(length = dim(hkgM)[2])
      if(calc == "arith") {
        for(i in 1:dim(hkgM)[2]) {
          if(! FALSE %in% is.na(hkgM[,i])) next # go to next sequence if we have only NAs
	  hkgV[i] <- mean(hkgM[,i], na.rm=TRUE)
        }
      }
      else {
        for(i in 1:dim(hkgM)[2]) {
          if(! FALSE %in% is.na(hkgM[,i])) next # go to next sequence if we have only NAs
          hkgV[i] <- geomMean(hkgM[,i], na.rm=TRUE)
        }
      }
    } 
    else {
      if(TRUE %in% is.na(hkgM)) {
       warning("NAs present in housekeeping gene readings")
       if(! FALSE %in% is.na(hkgM)) stop("Need at least 1 non NA for the housekeeper")
      }
      hkgV <- hkgM # Because it's a vector really anyway: we only have one NA
    }
    exprs(qPCRBatch) <- t(t(exprs(qPCRBatch)) - hkgV)
    return(qPCRBatch)
  }
)
