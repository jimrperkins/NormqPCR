setMethod("deltaDeltaCt", signature = "qPCRBatch", definition =
  function(qPCRBatch, maxNACase=0, maxNAControl=0, hkgs, contrastM, case, control, paired=TRUE, hkgCalc="arith", statCalc="arith") {
    hkgs <- make.names(hkgs)
    if(length(hkgs) < 1) stop("Not enough hkgs given")

    if(FALSE %in% (hkgs %in% featureNames(qPCRBatch))) stop("invalid housekeeping gene given")
    for(hkg in hkgs){
        if(sum(is.na(hkg)) > 0) warning(hkg, " May be a bad housekeeping gene to normalise with since it did not produce a reading ", sum(is.na(hkg)), "times out of", length(hkg))
    }
    cases <- row.names(contrastM)[contrastM[,case] == 1]
    controls <- row.names(contrastM)[contrastM[,control] == 1]
    expM <- exprs(qPCRBatch)
    caseM <- expM[,cases]
    controlM <- expM[,controls]

    if(length(hkgs) > 1) {
	hkgMCase <- caseM[hkgs, ]
        hkgMControl <- controlM[hkgs, ]
	if(hkgCalc == "arith") {
		hkgVCase <- apply(hkgMCase, 2, mean, na.rm=TRUE)
		hkgVControl <- apply(hkgMControl, 2, mean, na.rm=TRUE)
	} else {
		hkgVCase <- apply(hkgMCase, 2, geomMean, na.rm=TRUE)
		hkgVControl <- apply(hkgMControl, 2, geomMean, na.rm=TRUE)
	}
    } else { # Just use the first HKG
        hkgVCase <- caseM[hkgs[1], ]
        hkgVControl <- controlM[hkgs[1], ]
    }

    sdHkgCase <- sd(hkgVCase, na.rm=TRUE)
    sdHkgControl <- sd(hkgVControl, na.rm=TRUE)

    if(! FALSE %in% is.na(hkgVCase) || ! FALSE %in% is.na(hkgVControl)) stop("Need at least 1 non NA for the housekeeper")

    ddCts <- vector(length=length(featureNames(qPCRBatch)))
    dCtCases <- vector(length=length(featureNames(qPCRBatch)))
    dCtControls <- vector(length=length(featureNames(qPCRBatch)))
    sdCtCases <- vector(length=length(featureNames(qPCRBatch)))
    sdCtControls <- vector(length=length(featureNames(qPCRBatch)))
    minddCts <- vector(length=length(ddCts))
    maxddCts <- vector(length=length(ddCts))

    i <- 1
    for (detector in featureNames(qPCRBatch)) {
        VCase <- caseM[detector,]
        VControl <- controlM[detector,]
        if(length(VCase) == 1) {
          warning("Only one Detector for Control")
          dCtCase <- VCase
          sdCase <- NA
        } else if(! FALSE %in% is.na(VCase)) {
          warning("No Detector for Case")
          dCtCase <- rep(NA, length = VCase)
          sdCase <- NA
        } else {
          if(statCalc == "geom") {
            dCtCase <- mean(2^-(VCase - hkgVCase), na.rm=TRUE)
            sdCase <- sd(2^-(VCase - hkgVCase), na.rm=TRUE)
          }
          if(statCalc == "arith") {
            dCtCase <- mean(VCase, na.rm=TRUE) - mean(hkgVCase, na.rm=TRUE)
	    if (paired == TRUE) {
	      sdCase <- sd(VCase - hkgVCase, na.rm=TRUE)
	    } else  {
	      sdCase <- sqrt(sd(VCase, na.rm=TRUE)^2 + sdHkgCase^2)
	    }
          }
        }

        if(length(VControl) == 1) {
          warning("Only one Detector for Control")
          dCtControl <- VControl
          sdControl <- NA
        } else if(! FALSE %in% is.na(VControl)) {
          warning("No Detector for Control")
          dCtControl <- rep(NA, length = VControl)
          sdControl <- NA
        } else {
          if(statCalc == "geom") {
            dCtControl <- mean(2^-(VControl - hkgVControl), na.rm=TRUE)
            sdControl <- sd(2^-(VControl - hkgVControl), na.rm=TRUE)
          }
          if(statCalc == "arith") {
            dCtControl <- mean(VControl, na.rm=TRUE) - mean(hkgVControl, na.rm=TRUE)
            if (paired == TRUE) {
              sdControl <- sd(VControl - hkgVControl, na.rm=TRUE)
            } else  {
              sdControl <- sqrt(sd(VControl, na.rm=TRUE)^2 + sdHkgControl^2)
            }
          }
        }
        if(sum(is.na(VCase)) > maxNACase) {
          dCtCase <- NA
        }
        if(sum(is.na(VControl)) > maxNAControl) {
          dCtControl <- NA
        }
        if(statCalc == "geom") {
          ddCt <- dCtCase / dCtControl
        }
        if(statCalc == "arith") {
          ddCt <- dCtCase - dCtControl
        }
        if(is.na(ddCt)) {
          if(is.na(dCtCase) && ! is.na(dCtControl)) ddCt <- "-"
          else if(is.na(dCtControl) && ! is.na(dCtCase)) ddCt <- "+"
          else if(is.na(dCtControl) && is.na(dCtCase)) ddCt <- NA
          minddCts[i] <- NA
          maxddCts[i] <- NA
          ddCts[i] <- ddCt
        }
        else {
          if(is.na(sdCase)) {
            minddCts[i] <- NA
            maxddCts[i] <- NA
          }
          else {
            if(statCalc == "geom") {
              minddCts[i] <- NA
              maxddCts[i] <- NA
              ddCts[i] <- ddCt
            }
            if(statCalc == "arith") {
              minddCts[i] <- 2 ^ -(ddCt + sdCase)
              maxddCts[i] <- 2 ^ -(ddCt - sdCase)
              ddCts[i] <- 2^-ddCt
            }

          }
        }
        if(statCalc == "geom") {
	  dCtCases[i] <- dCtCase
	  sdCtCases[i] <- sdCase
	  dCtControls[i] <- dCtControl
	  sdCtControls[i] <- sdControl
        }
        if(statCalc == "arith") {
	  dCtCases[i] <- 2^-dCtCase
	  sdCtCases[i] <- sdCase
	  dCtControls[i] <- 2^-dCtControl
	  sdCtControls[i] <- sdControl

        }
        i <- i+1
    }

    ddCtTable <- as.data.frame(cbind(featureNames(qPCRBatch),format(dCtCases, digits=4),format(sdCtCases, digits=4),format(dCtControls, digits=4),format(sdCtControls, digits=4),format(ddCts,digits=4),format(minddCts, digits=4),format(maxddCts, digits=4)))
    names(ddCtTable) <- c("ID", paste("2^-dCt",case,sep="."), paste(case,"sd",sep="."), paste("2^-dCt",control,sep="."), paste(control,"sd",sep="."),"2^-ddCt","2^-ddCt.min", "2^-ddCt.max")
    return(ddCtTable)
  }
)
