
##  To run:
##  L <- SamSPECTRAL(FC.FCSRobj ,dimensions)

# Call libraries
# library(flowCore)

# Call functions:
.First.lib <-
function(lib, pkg)
{
 library.dynam("SamSPECTRAL", package = pkg, lib.loc = lib)
 provide(SamSPECTRAL) 
 return(invisible(0))
}

   	#dyn.load("SamSPECTRAL.so")                     # C function, packaging
    #dyn.load("maximum_of_rows.so")				    # C function
    #dyn.load("conductance_computation.so")			# C function
    #source("functions1.R")
    #packaging
    #source("Reading_Data.R")
    #source("Building_Communities.R")
    #source("Conductance_Calculation.R")
    #source("Civilized_Spectral_Clustering.R")
    #source("Connecting.R")

# Setting Parameters ----> Parameters are being set in the following R code.
#source("R/Parameters.R")


### MAIN function ###
	################################################################ S T A R T ########################################################################
SamSPECTRAL <- function (data.points, dimensions=1:dim(data.points)[2], normal.sigma, separation.factor, number.of.clusters=NA, 
						scale=rep(1,dim(data.points)[2]),talk=TRUE, precision=6, eigenvalues.num=NA, return_only.labels=TRUE, 
						do.sampling=TRUE, beta=4,stabilizer=500){ 

	if( class(data.points)!="matrix" ){  #Check if there are enough number of points, stored in the data.matrix.
	    if(talk) message("BAD input for SamSPECTRAL!, maybe it is not a matrix.")
		return(NULL)
	}# else, the input is good.    

	# Setting internal parameters: (These are rarely needed to be changed.)
		community.weakness.threshold <- 1     			# All communities with population less than this, will be omitted.
		maximum.number.of.clusters <- 30			
		# We use this parameter to find the "knee spot" and then estimate the bumber of clusters for the spectral clustering algorithm.
		m <- 3000				# upper bound on the number of clusters (almost)
		# precision <- 6		# In R, .Machine$double.eps is about 10^(-16) so a precision of 10^(-6) will be safe in our calculations. => precision=6.
		# iterations <- 200		# Maximum number of iterations for k-means clustering		


    # Read data files and transform them using log transform
	#Data <- FSC.file@exprs    
	#data.mtrix	<- log10(Data[,dimensions])
	data.matrix	<- as.matrix(data.points[,dimensions])	
		# If dim of the above matrix is n*1, then R will do stupid trasformation from matrix to vector.

	# Start reporting:
    t.SamSPECTRAL<-Sys.time()
	if(talk) message("***********************************************")
	if(talk) message(paste("SamSPECTRAL started at: ", t.SamSPECTRAL))
	if(talk) message("***********************************************")

	if(talk) message(dim(data.matrix))
    if(talk) message("^^^^^^^^^^^^^^ Number of original points and clustering dimentionality.")



	#Scaling all the dimensions such that each coordinate of full is in range: [0,1].
		full <- data.matrix
		for (i.column in 1:dim(data.matrix)[2]){#For all columns
			ith.column <- data.matrix[,i.column]
			ith.column.length <- max(ith.column) - min(ith.column)
			if( ith.column.length!=0) 
				full[,i.column] <- ((ith.column - min(ith.column) )/ith.column.length) * scale[i.column]   # This is the scaled column.
		}#End for (i.column.
	# Therefor, 	
		space.length <- 1

    # Sample the data and build the communities
	if(talk) message("Faithful sampling...")
    society <- Building_Communities(full,m, space.length, community.weakness.threshold, talk=talk, do.sampling=do.sampling)
    
    # Compute conductance between communities
	if(talk) message("Conductance computation...")
    conductance <- Conductance_Calculation(full, normal.sigma, space.length, society, precision, talk=talk, beta=beta)
    
    # Use spectral clustering to cluster the data
    clust_result <- Civilized_Spectral_Clustering(full, maximum.number.of.clusters, society, conductance, 
					number.of.clusters=number.of.clusters, talk=talk, eigenvalues.num = eigenvalues.num, stabilizer=stabilizer)
    
    number.of.clusters <- clust_result@number.of.clusters
    labels.for_num.of.clusters <- clust_result@labels.for_num.of.clusters
    
    # Connect components
    component.of <- Connecting(full, society, conductance, number.of.clusters, labels.for_num.of.clusters, separation.factor, talk=talk)$label


	################################################################# E N D ########################################################################
	if(talk) message("Total time:"); if(talk) message(Sys.time()-t.SamSPECTRAL)
	if (return_only.labels){	  
	    return(component.of)
	} else {
	    return(list(labels=component.of, society=society, conductance=conductance, clustering_result=clust_result) )
	}
}


















