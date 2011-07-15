################################################################################
## Vandesompele et al. (2002): gene-stability measure M
################################################################################
## x: matrix or data.frame containing data
## log: logical, data on log-scale
## na.rm: remove NA values
stabMeasureM <- function(x, log = TRUE, na.rm = TRUE){
    if(!is.data.frame(x) & !is.matrix(x))
        stop("'x' has to of class matrix or data.frame")

    if(is.data.frame(x)) x <- data.matrix(x)
    
    n <- ncol(x)
    if(n == 1) 
        stop("you need at least two variables (i.e., columns) for this computation")

    M <- numeric(n)
    for(j in 1:n){
        if(log)
            A <- x[,j] - x[,-j]
        else
            A <- log2(x[,j]/x[,-j])
    if(n > 2){
        N <- colSums(!is.na(A))
        N[N < 1] <- NA
        Mean <- colMeans(A, na.rm = na.rm)
        M[j] <- mean(sqrt(rowSums((t(A) - Mean)^2, na.rm = na.rm)/(N - 1)))
    }else
        M[j] <- sd(A, na.rm = na.rm)
    }
    names(M) <- colnames(x)
    M
}
