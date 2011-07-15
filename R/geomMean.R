## Compute geometric mean for numeric vector
## x: numeric vector of non-negative Reals
## na.rm: remove NA values
geomMean <- function(x, na.rm = TRUE){
  if(!is.numeric(x) && !is.complex(x) && !is.logical(x)){
    warning("argument is not numeric or logical: returning NA")
    return(as.numeric(NA))
  }
  if(all(is.na(x))){
    warning("no non-NAs in 'x'"); return (NA)
  }
  if(na.rm) x <- x[!is.na(x)]
  if(length(x) == 0)
    stop("x is a vector of length 0")
  if(any(x < 0))
    stop("'x' contains negative value(s)")

  return(prod(x)^(1/length(x)))
}
