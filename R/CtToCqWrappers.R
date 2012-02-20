deltaDeltaCq <- function(qPCRBatch, maxNACase, maxNAControl, hkgs, contrastM, case, control, paired, hkgCalc, statCalc){
  deltaDeltaCt(qPCRBatch, maxNACase=0, maxNAControl=0, hkgs, contrastM, case, control, paired=TRUE, hkgCalc="arith", statCalc="arith")
}

deltaCq <- function(qPCRBatch, hkgs, combineHkgs, calc){
  deltaCt(qPCRBatch, hkgs, combineHkgs=FALSE, calc="arith")
}
