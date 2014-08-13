#######################################################################
## Compute efficiencies
#########################################################################

## calculate efficiencies and Cq values for CyclesSet
setMethod("CqValues", signature = "CyclesSet", definition = 
            function(object, Effmethod = "expfit", group = NULL, model = l5, 
                     check = "uni2", checkPAR = parKOD(), remove = "none", 
                     exclude = NULL, type = "cpD2", labels = NULL, 
                     norm = FALSE, baseline = FALSE, basefac = 1, 
                     smooth = NULL, smoothPAR = list(span = 0.1), factor = 1, 
                     opt = FALSE, optPAR = list(sig.level = 0.05, crit = "ftest"), 
                     plot = FALSE, verbose = FALSE, ...){
                if(missing(Effmethod))
                    Effmethod <- "sigfit"
                else{
                    if(length(Effmethod) > 1){
                        Effmethod <- Effmethod[1]
                        warning("Only first element of 'Effmethod' is used!")
                    }
                }
                if(Effmethod != "sigfit")
                    methods <- c("sigfit", Effmethod)
                else
                    methods <- Effmethod
                
                fluoData <- exprs(object)
                x <- data.frame(Cycles = fData(object)[,"Cycle number"], fluoData)
                names(x) <- c("Cycles", pData(object)[,"Sample position"])
                
                if(missing(group))
                    
                res <- pcrbatch(x = x, cyc = 1, fluo = NULL, methods = methods,
                                model = model, check = check, checkPAR = checkPAR, 
                                remove = remove, exclude = exclude, type = type, 
                                labels = labels, 
                                norm = norm, baseline = baseline, basefac = basefac, 
                                smooth = smooth, smoothPAR = smoothPAR, 
                                factor = factor, opt = opt, optPAR = optPAR, 
                                group = group, 
                                names = "first", plot = plot, verbose = verbose, ...)
                ## extract Cq values and efficiencies
                CqVals <- as.matrix(as.numeric(res[res[,"Vars"] == paste("sig", type, sep = "."),-1]))
                if(Effmethod == "sigfit"){
                    Effs1 <- res[res[,"Vars"] == "sig.eff", -1]
                    se.Effs1 <- as.matrix(sqrt(as.numeric(res[res[,"Vars"] == "sig.resVar", -1])))
                }
                if(Effmethod == "sliwin"){
                    Effs1 <- res[res[,"Vars"] == "sli.eff", -1]
                    se.Effs1 <- NULL
                }
                if(Effmethod == "expfit"){
                    Effs1 <- res[res[,"Vars"] == "exp.eff", -1]
                    se.Effs1 <- as.matrix(sqrt(as.numeric(res[res[,"Vars"] == "exp.resVar", -1])))
                }
                if(Effmethod == "LRE"){
                    Effs1 <- res[res[,"Vars"] == "LRE.eff", -1]
                    se.Effs1 <- NULL
                }
                Effs1 <- as.matrix(as.numeric(Effs1))
                metData <- data.frame(labelDescription = res[,"Vars"],
                                      row.names = res[,"Vars"], check.names = FALSE,
                                      stringsAsFactors = FALSE)
                pD <- t(res[,-1])
                colnames(pD) <- res[,"Vars"]
                metData1 <- varMetadata(phenoData(object))
                metData <- rbind(metData1, metData)
                pData(object) <- cbind(pData(object), pD)
                varMetadata(phenoData(object)) <- metData
#                effs <- assayDataNew("environment", effs = Effs1)
#                se.effs <- assayDataNew("environment", effs = se.Effs1)
                res <- new("qPCRBatch", exprs = CqVals, featureData = phenoData(object))
                effs(res) <- Effs1
                se.effs(res) <- se.Effs1
                res
            }
)
