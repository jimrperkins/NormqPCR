## x: matrix or data.frame containing the data
## group: factor
## log: data on log-scale
## na.rm: remove NA values

setMethod("stabMeasureRho", signature(x = "matrix"), definition =
  function(x, group, log = TRUE, na.rm = TRUE, returnAll = FALSE){
    if(is(x, "qPCRBatch")) {
        x <- t(exprs(x))
    }

    if(!is.data.frame(x) & !is.matrix(x))
        stop("'x' has to be of class qPCRBatch, matrix or data.frame")

    if(is.data.frame(x)) x <- data.matrix(x)
    k <- ncol(x) # number of variables
    n <- nrow(x) # number of samples
    if(k == 1) 
        stop("you need at least two variables (i.e., columns) for this computation")

    if(!log) x <- log2(x)
    if(!is.factor(group)){
        group <- factor(group)
        warning("Argument 'group' is transformed to a factor vector.")
    }
    m <- nlevels(group) # number of groups
    Levels <- levels(group)
    ngr <- as.integer(table(group))
    
    mej <- rowMeans(x, na.rm = na.rm)
    
    if(m == 1){
        warning("There is only one group. Only the variance estimates are returned!")
        mei <- colMeans(x, na.rm = na.rm)
        me <- mean(mej, na.rm = na.rm)
        N <- colSums(!is.na(x))
        N[N < 1] <- NA
        a <- rowSums((t(x-mej)-mei+me)^2, na.rm = na.rm)/(N-1)
        b <- sum(a, na.rm = na.rm)
        var.no.group <- (a - b/(k*(k-1)))/(1-2/k)
        return(var.no.group)
    }else{
        meigr <- as.matrix(aggregate(x, by = list(group), mean, na.rm = na.rm)[,-1])
        megr <- rowMeans(meigr, na.rm = na.rm)
        x.split <- split(x, f = group)
        var.group.all <- matrix(0, ncol = k, nrow = m)
        for(i in 1:m){
            x.temp <- matrix(x.split[[i]], nrow = ngr[i])
            N <- colSums(!is.na(x.temp))
            N[N < 1] <- NA
            a <- colSums((t(t(x.temp)-meigr[i,])-mej[group == Levels[i]]+megr[i])^2, na.rm = na.rm)/(N-1)
            b <- sum(a, na.rm = na.rm)
            var.group.all[i,] <- (a - b/(k*(k-1)))/(1-2/k)
        }
        if(any(var.group.all < 0)) var.group.all <- pmax(var.group.all, 0)
        
        m1i <- colMeans(meigr, na.rm = na.rm)
        m1j <- rowMeans(meigr, na.rm = na.rm)
        m1 <- mean(m1j, na.rm = na.rm)
        dif <- t(meigr - m1j) - m1i + m1
        va <- var.group.all/ngr
        tau <- max(sum(dif*dif)/((m-1)*(k-1))-mean(va), 0)
        dnew <- dif*tau/(tau+t(va))
        vanew <- t(va+tau*va/(tau+va))
        rownames(vanew) <- rownames(dnew)
        qm <- abs(dnew)+sqrt(vanew)
        qmaal <- rowMeans(qm)
        if(returnAll)
            return(list(rho = qmaal, d = dnew, v = vanew))
        else
            return(qmaal)
    }
  }
)

setMethod("stabMeasureRho", signature(x = "qPCRBatch"), definition =
  function(x, group, log = TRUE, na.rm = TRUE, returnAll = FALSE){
    x <- t(exprs(x))
    return(stabMeasureRho(x, group, log, na.rm, returnAll))
  }
)
