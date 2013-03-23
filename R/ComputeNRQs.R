## Compute NRQs
setMethod("ComputeNRQs", signature = "qPCRBatch", definition = 
            function(qPCRBatch, hkgs){
                warning("This function is still experimental!")
                
                RefCq <- colMeans(exprs(qPCRBatch))
                if(is.null(se.exprs(qPCRBatch)))
                  warning("'se.exprs(qPCRBatch)' is NULL! You could use function 'deltaCt' or 'deltaDeltaCt'.")
                SD.ofCq <- se.exprs(qPCRBatch)
                
                dCq <- t(RefCq-t(exprs(qPCRBatch)))
                
                if(is.null(effs(qPCRBatch)))
                  stop("Efficiencies are missing! You could use function 'deltaCt' or 'deltaDeltaCt'.")
                Ejl <- as.vector(effs(qPCRBatch))
                
                if(is.null(se.effs(qPCRBatch)))
                  warning("'se.effs(qPCRBatch)' is NULL! You could use function 'deltaCt' or 'deltaDeltaCt'.")
                SEjl <- as.vector(se.effs(qPCRBatch))
                
                RQ <- exp(dCq*log(Ejl))
                SD.ofRQ <- sqrt(RQ^2*((dCq*SEjl/Ejl)^2 + (log(Ejl)*SD.ofCq)^2))
                
                ind.HKs <- rownames(RQ) %in% hkgs
                NF <- apply(RQ[ind.HKs,], 2, geomMean)
                SD.ofNF <- NF*sqrt(mean((SD.ofRQ[ind.HKs,]/RQ[ind.HKs,])^2))
                
                NRQ <- t(t(RQ)/NF)
                SD.ofNRQ <- NRQ*sqrt((SD.ofNF/NF)^2 + (SD.ofRQ/RQ)^2)
                
                exprs(qPCRBatch) <- NRQ
                se.exprs(qPCRBatch) <- SD.ofNRQ
                qPCRBatch
            })
