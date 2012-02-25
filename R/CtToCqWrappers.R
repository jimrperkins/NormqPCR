deltaDeltaCq <- function(qPCRBatch, maxNACase=0, maxNAControl=0, hkgs, contrastM, case, control, 
                         paired=TRUE, hkgCalc="arith", statCalc="arith"){
  deltaDeltaCt(qPCRBatch=qPCRBatch, maxNACase=maxNACase, maxNAControl=maxNAControl, hkgs=hkgs, 
               contrastM=contrastM, case=case, control=control, paired=paired, hkgCalc=hkgCalc, 
               statCalc=statCalc)
}

deltaCq <- function(qPCRBatch, hkgs, combineHkgs=FALSE, calc="arith"){
  deltaCt(qPCRBatch=qPCRBatch, hkgs=hkgs, combineHkgs=combineHkgs, calc=calc)
}
