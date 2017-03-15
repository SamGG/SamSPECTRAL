## INTENTION:
## Here, we get the selected points (~2000 points) with the conductance matrix
## to build the kernel matrix. Then, we'll apply the spectral algorithm.
##
## IMPORTANT: For heavy nodes we have to consider that "they represent a community" so the similarity should be computed
## by a very shorter distance. In other words, we are looking for the conductance between two communities or a estimation for that.
## Then, for the heavy points, we MAY multiply the weights on their edges by their density (or take other strategies).
##------------------------------------------------------------------------------------------------------------------------


## The class 
##defined for output of the function in this file, namely: Civilized_Spectral_Clustering().=
##---------------------------------------------------------------------------------------------------------
## ClusteringResult
## The eigen function used to output an object of class "list" in R 3.3.0 but this may  change to "eigne" class according to the current devel verion of R.
## --Habil Zare, 2017-03-15
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
setClass("ClusteringResult",
         representation(labels.for_num.of.clusters="list", number.of.clusters="numeric", eigen.space=class(eigen(cbind(1, c(1,-1)))), minimum.degree.sum="numeric"))
##===========================================================================================================================


## The function:
##--------------------------------------------------------------------------------------------------------------##
Civilized_Spectral_Clustering <- function(full, maximum.number.of.clusters=30, society, conductance, iterations=200,
                                          number.of.clusters="NA",
                                          k.for_kmeans="NA",minimum.eigenvalue="NA", minimum.degree=0,
                                          eigenvalues.num = NA, talk=TRUE, stabilizer=1000, one.line=FALSE,
                                          replace.inf.with.extremum=TRUE)
{

    t1<-Sys.time()
    if(talk) message(t1)
    ##--------------------------------------------## S T A R T ------------------------------------------------##


    ## society
    nbhood <- society$nbhood
    repres.indices <-society$representatives
    community <- society$communities
    num.of.cummnities <- length(community)


    ## Initialization
    n <- dim(full)[1] ## Total number of points.

    ## Input checking:
    checked <- check.SamSPECTRAL.input(data.points=full,replace.inf.with.extremum=replace.inf.with.extremum)
    full <- checked$data.matrix

    ## Weightening based on community conductance
    conductance.matrix <- conductance$conductance.matrix
    sigma <- conductance$sigma
    if(talk) message(paste("sigma = ", sigma))
    weight.matrix <- conductance.matrix
    ## We rename it just to make more sence for clustering algorithm, we could also call it similarity matrix.


    l <-c()	## This will be the matrix for which the eigenspace will be computed. We may need it to debuge.


    ##----------------------------------------## Runing the spectrul clustering algorith------------------------

    ## Runing the spectrul clustering algorithm
    ##if(talk) message("Weighted spectral clustering is runnig!!")

    ## Renaming varibles for the sake of easy use.
    x <- weight.matrix
    rown <- rownames(x)
    km <-x
    diag(km) <- 0

    ## Computing eigenspace
    if(talk) message("Computing eigenspace...")
    if(talk) message(Sys.time())
    
    d <- 1/sqrt(rowSums(km))
    isolated.community.indices <- which(d >= 1/minimum.degree)
    ##isolated.community.indices <- which(d==Inf)
    if(length(isolated.community.indices) >0){
	## Removing isolated outliers:
        ## If a point is so far from the rest, then it has little weight on all edges.
        ## Because we have limitted percision in R, then the corresponding rows and columns 
        ## will beceom zero. => computational problem.
        ## Solution: we consider these as "isolated outliers" and omit them for now.
        km <- km[-isolated.community.indices, -isolated.community.indices]
        d <- d[-isolated.community.indices]	
        if(talk) message(paste( length(isolated.community.indices), " communities are considered as isolated outliers." ));
    }		
    ##
    minimum.degree.sum <- min(1/d)
    if(talk) message("Minimum degree (sum of edges): ", minimum.degree.sum)		
    l <- d * km %*% diag(d)
    l <- (l + t(l))/2	## symetrizing might help to prevent complex numbers probably caused by round-off errors.
    ##
    ##
    ##*****************************************************
    eigen.space <-eigen(l)		## EIGEN SPACE COMPUTATION.
    ##*****************************************************
    ##
    if(talk) message(Sys.time())
    ##	
    ##----------------------------------------## K_MEANS --------------------------------------------------------##
    ##try({sample(1:100,10); rm(.Random.seed,inherits=TRUE)},silent=TRUE)	
    ## To remove the random seed which might be loaded from previouse workspace.	
    labels.for_num.of.clusters <- list()
    ## Each element contains the labels for that number of clusters.		
    
    ##------------------------ Estimating the number of clusters based on the "knee spot" ------------------------
    if(one.line){ ## It needs to be determined "Automatically".
        ## linear regression with 1 line.
        x.range <- maximum.number.of.clusters:(2*maximum.number.of.clusters)
        slop <- eigen.space$values[x.range]
        ## Built in R regresion	---> The line in known: ax+b
        f <- lm(slop ~ x.range)
        coeff <- coef(f); 
        a <- as.numeric(as.character(coeff[2]))
        b <- as.numeric(as.character(coeff[1]))	
    
        ## way(1): Finiding the last point in left of the regresion line in the graph of eigenvalues.
        ##  If a point is bellow the line, y-(ax+b) will be negative for that point.
        last.value <-1
        ## way(2); We have: a * x_1 + b = 1 =>  x_1 = (1-b)/a, so:
        ## number.of.clusters <- round((1-b)/a); if(talk) message("way(2), after regresion")
    } else {
        ## linear segmented regression with 2 lines.
        x.range <- 1:(2*maximum.number.of.clusters)
        knee.points <- eigen.space$values[x.range]
        kneepointDetection.result <- kneepointDetection(knee.points) 
        line1 <- kneepointDetection.result$l1
        line2 <- kneepointDetection.result$l2
        coeff <- coefficients(line2); 
        a <- as.numeric(as.character(coeff[2]))
        b <- as.numeric(as.character(coeff[1]))	
        
        ## Finiding the last point in left of the regresion line in the graph of eigenvalues.
        ##  If a point is bellow the line, y-(ax+b) will be negative for that point.
        last.value <- kneepointDetection.result$MinIndex
    }
    
    ## Finding the knee point based on the line(s):
    while (eigen.space$values[last.value] - (a*last.value +b ) < 0 )
        last.value <- last.value +1 ## We have not reached the line yet.


    if(number.of.clusters=="NA"){ ## It needs to be determined "Automatically".
        if(minimum.eigenvalue!="NA"){	## We pick all the values above this threshold.
            large.eigenvalues.num <- length(which(eigen.space$values>=minimum.eigenvalue))
            number.of.clusters <- large.eigenvalues.num 
        } else {	## We use the segmented regression.	
            number.of.clusters <- last.value+1										
            ## This is the first value in right of the regesion line.
            if(number.of.clusters <15){ 
                if(talk) message("number.of.clusters is artificially forced to be >= 15")
                number.of.clusters <- 15
            }##End if(number.of.clusters <15). 	
            if(talk) message("number.of.clusters is automatically calculated to be: ")
            if(talk) message(paste("----------------------------------------------------> ", number.of.clusters))
        }##End elseif(!is.na(minimum.eigenvalue)).	
    }##End if(is.na(number.of.clusters)).

    ##Showing eigen values:
    if (!is.na(eigenvalues.num)){
        title <- paste("normal.sigma= ",sigma, ", knee point: ", last.value+1, ", number of large values:",
                       number.of.clusters, sep='' ) 
        plot(eigen.space$values[1:eigenvalues.num], main=title)
        if(!one.line)
            abline(a=coef(line1)[1],b=coef(line1)[2],col='blue')
        abline(a=b,b=a,col='red')
        abline(v=last.value+1, col='green')
        abline(h=eigen.space$values[number.of.clusters], col='orange')
    }##End if (!is.na(eigenvalues.num)).
    ##----------------------------------------------------------------------------

    ## RUNNING KMEANS: ##
    ##
    for (i.number.of.clusters in 1: length(number.of.clusters)){	## For all different number of clusters,
        if (k.for_kmeans=="NA"){	## It is not given from outside.
            centers.num <- number.of.clusters[i.number.of.clusters]	
        } else {	
            centers.num <- k.for_kmeans
        }##End elseif (!is.na(k.for_kmeans)).  
        
        while (centers.num<n ){	## We will continue this until the output of kmeans is not "try_error".
            try_kmeans <- try(silent=TRUE, {
                if(talk) message("Runing kmeans...")
                xi <- eigen.space$vectors[,1:centers.num]
                yi <- xi/sqrt(rowSums(xi^2))
                res <- kmeans(yi, centers.num, iterations, nstart=stabilizer)
                ## Returning the output for this number of clusters, EXTRA.
                cent <- matrix(unlist(lapply(1:centers.num,ll<- function(l){colMeans(x[which(res$cluster==l),])})), 
                               ncol=dim(x)[2], byrow=TRUE)
                withss <- unlist(lapply(1:centers.num,ll<- function(l){sum((x[which(res$cluster==l),] - cent[l,])^2)}))
                names(res$cluster) <- rown
            })##End try.
            if (is(try_kmeans,"try-error")){
                centers.num <- centers.num + 1
                if(talk) message("centers.num was increased to: ",centers.num)				
            } else break ##End if (try_kmeans=="try-error").
        }##End while ( (try_kmeans=="start" | try_kmeans=="try-error" ) & centers.num<n ).
        number.of.clusters[i.number.of.clusters] <- centers.num	## This may be changed in above.
        ##----------------------------## End of spectral ------------------------------------


        ##------------------------##  start of saving results ------------------------##

        sc <- res$cluster	## This vector contains the cluster labels we get from spectral clustering.

        ## Some communities are deleted while removing isolated outliers, we label them by 0 here.
        community.labels <- c()
        num.of.added.isolated.communities <- 0
        for (i in 1:n){
            if (num.of.added.isolated.communities < length(isolated.community.indices)){
                ## there are more isolated communies to be added.
                if (isolated.community.indices[num.of.added.isolated.communities+1] != i)
                    ## i.e. the community i, has not been isolated,
                    community.labels[i] <- sc[i-num.of.added.isolated.communities]
                ## To see that this is true, 
                ## check it for when num.of.added.isolated.communities=0
                else{	## the i'th community has been considered as isolated,
                    ## so we label it by 0.
                    community.labels[i] <- 0
                    num.of.added.isolated.communities <- num.of.added.isolated.communities +1
                }
            }
            else	## We are done with isolated ones,
                community.labels[i] <- sc[i-num.of.added.isolated.communities]	
        }
	
        ## Labeling all members according to their community lable.
        num.of.comu <- length(society$representatives)
        
        label <- 0	## label of individual points
        for (i in 1:num.of.comu){
            label[society$communities[[i]]$members] <- community.labels[i]
        }
	
        ##if(talk) message(paste(centers, "clusters were distinguished!"))
        labels.for_num.of.clusters[[centers.num]] <- label

    }## End of for loop for clustering.


    result <- new('ClusteringResult', labels.for_num.of.clusters = labels.for_num.of.clusters, 
                  number.of.clusters = number.of.clusters, eigen.space = eigen.space,
                  minimum.degree.sum=minimum.degree.sum)
    ##----------------------------------------## E N D ----------------------------------------
    if(talk) message(Sys.time()-t1)

    return(result) 
}












