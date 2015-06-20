
## Call libraries
## library(flowCore)

## Call functions:
.First.lib <-
    function(lib, pkg)
{
    library.dynam("SamSPECTRAL", package = pkg, lib.loc = lib)
    provide(SamSPECTRAL) 
    return(invisible(0))
}


## C functions
###dyn.load("SamSPECTRAL.so")                
###dyn.load("maximum_of_rows.so")			 
###dyn.load("conductance_computation.so")

### MAIN function ###
####__________________________________ S T A R T ___________________________________________
SamSPECTRAL <- function (data.points, dimensions=1:dim(data.points)[2], normal.sigma, separation.factor,
                         number.of.clusters="NA", scale=rep(1,dim(data.points)[2]),talk=TRUE, precision=6,
                         eigenvalues.num=NA,
                         return_only.labels=TRUE, do.sampling=TRUE, beta=4,stabilizer=100,
                         k.for_kmeans="NA",maximum.number.of.clusters=30,m=3000,minimum.eigenvalue="NA",
                         previous.result=NULL,
                         replace.inf.with.extremum=TRUE, minimum.degree=0, one.line=FALSE, doOrderLabels=TRUE){ 
    ## INPUT
    ## maximum.number.of.clusters=30			
    ## We use this parameter to find the "knee spot" and then estimate 
    ##the number of clusters for the spectral clustering algorithm.
    ## m <- 3000		## upper bound on the number of clusters (almost)
    ## precision <- 6		## In R, .Machine$double.eps is about 10^(-16) so a precision of 10^(-6) 
    ##^ will be safe in our calculations. => precision=6.
    ## iterations <- 200		## Maximum number of iterations for k-means clustering	
    ## minimum.degree <- 0	## If a node in the graph has total edge sum less than this threshold, 
    ## it will be considered as an isolated community.

    ## Input checking:
    checked <- check.SamSPECTRAL.input(data.points=data.points,dimensions=dimensions,
                                       replace.inf.with.extremum=replace.inf.with.extremum)
    data.matrix <- checked$data.matrix

    ## keeping input.
    input= list()	## We would linke to return back all the input at the end.
    input[["data.points"]] <- data.points
    input[["dimensions"]] <- dimensions
    input[["normal.sigma"]] <- normal.sigma
    input[["separation.factor"]] <- separation.factor
    input[["number.of.clusters"]] <- number.of.clusters
    input[["scale"]] <- scale
    input[["talk"]] <- talk
    input[["precision"]] <- precision
    input[["eigenvalues.num"]] <- eigenvalues.num
    input[["return_only.labels"]] <- return_only.labels
    input[["do.sampling"]] <- do.sampling
    input[["beta"]] <- beta
    input[["stabilizer"]] <- stabilizer
    input[["k.for_kmeans"]] <- k.for_kmeans
    input[["maximum.number.of.clusters"]] <- maximum.number.of.clusters
    input[["m"]] <- m
    input[["minimum.eigenvalue"]] <- minimum.eigenvalue
    input[["minimum.degree"]] <- minimum.degree
    ## input[["previous.result"]] <-  previous.result	
    ## This will cause a recursive data production which is memory intensive
    
    
    ## Setting internal parameters: (These are rarely needed to be changed.)
    community.weakness.threshold <- 1
    ##^ All communities with population less than this, will be omitted.	


    ## Start reporting:
    t.SamSPECTRAL<-Sys.time()
    if(talk) message("***********************************************")
    if(talk) message(paste("SamSPECTRAL started at: ", t.SamSPECTRAL))
    if(talk) message("***********************************************")

    if(talk) message(dim(data.matrix))
    if(talk) message("^^^^^^^^^^^^^^ Number of original points and clustering dimentionality.")

    ##Scaling all the dimensions such that each coordinate of full is in range: [0,1].
    full <- data.matrix
    for (i.column in 1:dim(data.matrix)[2]){##For all columns
        ith.column <- data.matrix[,i.column]
        ith.column.length <- max(ith.column) - min(ith.column)
        if( ith.column.length!=0) 
            full[,i.column] <- ((ith.column - min(ith.column) )/ith.column.length) * scale[i.column]
        ##^ This is the scaled column.
    }##End for (i.column.
    ## Therefor, 	
    space.length <- 1

    ## Detrmining which computations to repeat:RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
    rebuild.communities <- TRUE
    recalculate.conductance <- TRUE
    re_civilized.spectral.clustering <- TRUE
    if(!is.null(previous.result)){ 
        if(previous.result$input$m == m){	
            ## We do not to repeat sampling
            society <- previous.result$society
            rebuild.communities <- FALSE
        }##End if.	
        if(!rebuild.communities &
           previous.result$input$normal.sigma == normal.sigma &
           previous.result$input$beta==beta &
           previous.result$input$precision==precision){	
            ## We do not to repeat connductance calculation
            conductance <- previous.result$conductance	
            recalculate.conductance <- FALSE
        }##End if.
        if(!recalculate.conductance & 
           previous.result$input$maximum.number.of.clusters == maximum.number.of.clusters &
           previous.result$input$k.for_kmeans ==k.for_kmeans&
           previous.result$input$number.of.clusters ==number.of.clusters &
           previous.result$input$stabilizer==stabilizer &
           previous.result$input$minimum.eigenvalue==minimum.eigenvalue){
            ## We do not repeat the civilized spectral clustering.	
            clust_result <- previous.result$clustering_result
            re_civilized.spectral.clustering <- FALSE			
        }##Enf if.
    }##End if(!is.null(previous.result)).RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR

    ## Sample the data and build the communities
    if(talk) message("Faithful sampling...")
    if(rebuild.communities) 
        society <- Building_Communities(full,m, space.length, community.weakness.threshold, talk=talk,
                                        do.sampling=do.sampling)
    
    ## Compute conductance between communities
    if(talk) message("Conductance computation...")
    if(recalculate.conductance)
    	conductance <- Conductance_Calculation(full, normal.sigma, space.length, society, precision, talk=talk, beta=beta)
    
    ## Use spectral clustering to cluster the data
    if(re_civilized.spectral.clustering)
        clust_result <- Civilized_Spectral_Clustering(full=full,maximum.number.of.clusters=maximum.number.of.clusters,
                                                      society=society,conductance=conductance, k.for_kmeans=k.for_kmeans,
                                                      number.of.clusters=number.of.clusters,
                                                      talk=talk, eigenvalues.num = eigenvalues.num, stabilizer=stabilizer,
                                                      minimum.eigenvalue=minimum.eigenvalue, minimum.degree=minimum.degree,
                                                      one.line=one.line)
    number.of.clusters <- clust_result@number.of.clusters
    labels.for_num.of.clusters <- clust_result@labels.for_num.of.clusters
    
    ## Connect components
    component.of <- Connecting(full, society, conductance, number.of.clusters, labels.for_num.of.clusters,
                               separation.factor, talk=talk)$label

    ## Ordering labels:
    labels <- component.of
    if(doOrderLabels){
        if(talk) message("Relabeling...")
        orderedComponents <- as.numeric(names(sort(table(component.of),decreasing=TRUE)))
        for(ind in 1:length(orderedComponents)){
            c1 <- orderedComponents[ind] ## the component
            labels[component.of==c1] <- ind
        }
    }
    sizes <- table(labels)
    timeTaken <- Sys.time()-t.SamSPECTRAL
    ##_______________________________ O U T P U T _______________________________________________
    if(talk) message("Total time:"); if(talk) message(timeTaken)
    if (return_only.labels){	  
        return(labels)
    } else {
        return(list(labels=labels, society=society, conductance=conductance, clustering_result=clust_result,
                    input=input,component.of=component.of, timeTaken=timeTaken, sizes=sizes))
    }
}
