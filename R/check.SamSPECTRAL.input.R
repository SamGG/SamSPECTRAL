check.SamSPECTRAL.input <- function(data.points,dimensions=1:ncol(data.points),replace.inf.with.extremum=FALSE){
    ## Input checking. data.points should be a matrix, without NA and Inf.
    result <- list()
    infinite <- FALSE
    if( class(data.points)!="matrix" ){  ##Check if there are enough number of points, stored in the data.matrix.
        stop("BAD input for SamSPECTRAL!, maybe it is not a matrix.")
    }## else, the input class is good.    

    ## Read data files and transform them using log transform
    data.matrix	<- data.points[,dimensions,drop=FALSE]
    
    ## Contains NA ?
    if (!prod(!is.na(data.matrix))){	## The input contains NA or NaN, so bad!
        stop("BAD input for SamSPECTRAL!, it expects only numbers not any NA or NaN.")
    }

    ## Contains +/- Inf ?
    if ( !prod(!is.infinite(data.matrix)) ){	## The input contains +/-Inf,
        if(!replace.inf.with.extremum ){
            stop("BAD input for SamSPECTRAL!, it expects only numbers not Inf or -Inf.")
        }
	infinite <- TRUE
        ##replacing with extremum,
        positive.inf.indices <- which(is.infinite(data.matrix) & data.matrix >0)
        negative.inf.indices <- which(is.infinite(data.matrix) & data.matrix <0)
        data.matrix[is.infinite(data.matrix)] <- NA
        max.value <- max(data.matrix, na.rm=TRUE)
        min.value <- min(data.matrix, na.rm=TRUE)
        data.matrix[positive.inf.indices] <- max.value
        data.matrix[negative.inf.indices] <- min.value
    }## else, the input is good.
    result[["data.matrix"]] <- data.matrix
    result[["dimensions"]] <- dimensions
    result[["infinite"]] <- infinite
    return(result)
}
