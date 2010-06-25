

### ### INTENTION:	######################################################
#
# 	Here, we look at the clusters 
# as outputs of civilized spectral clustering 
# and try to figure out which population is 
# divided to two or more by MISTAKE. 
#
# 	The CRITERIA we use is the connectivity between
# clusters compared to within connectivity.
# Technicaly, we look at some part of conductance matrix
# as an estimate for "connectivity".
#
# data structure:
# cluster[[]]	:	the list of indices for each cluster.
#					e.g. cluster[[3]] is th evector of indices
#					in cluster 3.
#
# components	:	a list, each element is a vector that includes 
#					the indices of clusters in each component 
#
###########################################################################





Connecting <- function(full, society, conductance, number.of.clusters, labels.for_num.of.clusters, separation.factor, talk=TRUE)
{

	# Compute indices of points in cluster j from the output of civilized spcetral clustering
	number.of.clusters.from.spectral <- number.of.clusters			

	### The clusters
	labels <-labels.for_num.of.clusters[[number.of.clusters.from.spectral]]
	repres.cluster	<-list()
	cluster<- list()
	
	for(i in 1:number.of.clusters.from.spectral){
		repres.cluster[[i]] <- which(labels[society$representatives]==i) 
		cluster[[i]] <- which(labels==i) 
	}

	### Loading the original data
		n <- dim(full)[1]						    
		dimention <-dim(full)[2]					

		# society
		nbhood <- society$nbhood
		repres.indices <-society$representatives
		community <- society$communities
		num.of.cummnities <- length(community)


		# Weightening 
		conductance.matrix <- conductance$conductance.matrix
		sigma <- conductance$sigma
		#if(talk) message(paste("sigma = ", sigma))



	
	max.conductance <- matrix( ,number.of.clusters.from.spectral, number.of.clusters.from.spectral)
	for (i in 1:number.of.clusters.from.spectral)
		for (j in 1:number.of.clusters.from.spectral)
			max.conductance[i,j] <- max( conductance.matrix[repres.cluster[[i]],repres.cluster[[j]]])

	connection <- max.conductance

	# Components
	components.list<-list()
	for (i in 1:dim(connection)[2])
		components.list[[i]] <-c(i)				
	
	concat.two.sided <- function(input.graph.matrix){

		# input check
		if(class(input.graph.matrix) != "matrix"){ 	
			result <- list (input.graph.matrix, components.list)
			return(result)	
		}
		output.graph.matrix <- input.graph.matrix
		within.connection <- diag(output.graph.matrix)


		# In the following lines, we will substitue the zero entries with a small value to avoid deviding by zero or producing NaN.
			absolutly.poor.clusters.index	<- which(within.connection==0)						# The clusters with zero within connection. 
		
		minumum.positive.connection <-  min(output.graph.matrix[output.graph.matrix!=0])
			# This is the smallest positive connection in the whole graph. 			
	
		substitute.small.value <- minumum.positive.connection * min(separation.factor, 1)/2	
			# We reduce it even more in case separation factor is so small and we
			# still would like to connect the two nodes.			
																										
		within.connection[absolutly.poor.clusters.index] <- substitute.small.value			# Replacing the zeros with the small value.
		#End of substituting.


		normlizer.matrix.for_connection <- diag(1/within.connection)
		normalized.output.graph.matrix <- output.graph.matrix %*% normlizer.matrix.for_connection	
		last.row <-input.graph.matrix[dim(input.graph.matrix)[1],]			
		for (i in 1:dim(output.graph.matrix)[2] ){	
			for (j in 1:dim(output.graph.matrix)[2] ){
				if ( (normalized.output.graph.matrix[i,j] >=1)  && (normalized.output.graph.matrix[j,i] >=1) && (i!=j) ) {	
					the.two.components <- c(last.row[i], last.row[j])
					host <- max(the.two.components) 				
					guest <- min(the.two.components) 				
					components.list[[host]] <- union( components.list[[host]], components.list[[guest]])
					components.list[[guest]] <- "vanished"
					output.graph.matrix[j,] <- pmax(output.graph.matrix[j,] , output.graph.matrix[i,])			
					output.graph.matrix[,j] <- pmax(output.graph.matrix[,j] , output.graph.matrix[,i])			
					output.graph.matrix <- output.graph.matrix[-i,]								
					output.graph.matrix <- output.graph.matrix[,-i]								
					result <- list (output.graph.matrix, components.list)
					return(result)
				}
			}
		}
		result <- list (output.graph.matrix, components.list)
		return(result)
	}


	
	concat <- function(input.graph.matrix){
    	## input check
    	if(class(input.graph.matrix) != "matrix"){	
    		result <- list (input.graph.matrix, components.list)
    		return(result)	
    	}
    
    	output.graph.matrix <- input.graph.matrix
    	number.of.columns <- dim(input.graph.matrix)[1]-1
    	last.row <-input.graph.matrix[dim(input.graph.matrix)[1],]	
    
    	for (i in 1:dim(output.graph.matrix)[2] ){	
    		matrix.without.last.row <- output.graph.matrix[1:number.of.columns,]		
			within.connection <- diag(matrix.without.last.row)									# Describes the connections inside the clusters.

			# In the following lines, we will substitue the zero entries with a small value to avoid deviding by zero or producing NaN.
				absolutly.poor.clusters.index	<- which(within.connection==0)						# The clusters with zero within connection. 
	
				minumum.positive.connection <-  min(matrix.without.last.row[matrix.without.last.row!=0])
					# This is the smallest positive connection in the whole graph. 			
	
				substitute.small.value <- minumum.positive.connection * min(separation.factor, 1)/2	
					# We reduce it even more in case separation factor is so small and we
					# still would like to connect the two nodes.			
																										
				within.connection[absolutly.poor.clusters.index] <- substitute.small.value			# Replacing the zeros with the small value.
			#End of substituting.

    		normlizer.matrix.for_connection <- diag(1/within.connection)
    		normalized.output.graph.matrix <- matrix.without.last.row %*% normlizer.matrix.for_connection			
    											
    		out.going.edges <- which(normalized.output.graph.matrix[,i] >= separation.factor )		
    																					
    		out.going.edges <- out.going.edges[ which( !(out.going.edges ==i)) ]		
    
    		## strategy
    		if (length(out.going.edges) >=1){
    			connection.to.others <- normalized.output.graph.matrix[,i]			
    			connection.to.others[i] <- -Inf
    			host.column <- which.max(connection.to.others)
    			host.row <- host.column	
    			host <- last.row[host.column]
    			guest.column <- i
    			guest <- last.row[i]		
    			components.list[[host]] <- union(components.list[[host]] , components.list[[guest]])
    			components.list[[guest]] <- "vanished"
    			ith.row <- output.graph.matrix[i,]
    			ith.column <- output.graph.matrix[,i]
    			output.graph.matrix[host.column,] <- pmax(output.graph.matrix[host.column,] , ith.row )		
    			output.graph.matrix[,host.row] <- pmax(output.graph.matrix[,host.row] , ith.column )		
    			output.graph.matrix[dim(output.graph.matrix)[1] , host.column] <- host		
    			output.graph.matrix	<-output.graph.matrix[-i,]			
    			output.graph.matrix	<-output.graph.matrix[,-i]			
       			result <- list (output.graph.matrix, components.list)
    			return(result)
    		}
        }
	result <- list (output.graph.matrix, components.list)
	return(result)
    }



	connection.with.indices	<- rbind(connection, 1:dim(connection)[2])	
	to.be.concatenated <- connection.with.indices	
	if(TRUE){
		repeat {	
			repeat {
				concat.result <-concat.two.sided(to.be.concatenated)	
				concatenated.graph.matrix <- concat.result[[1]]
				components.list <- concat.result[[2]]
				if ( sum(to.be.concatenated)  == sum(concatenated.graph.matrix) ) break
				to.be.concatenated <- concatenated.graph.matrix
			} 

			## separation.factor	
			concat.result <-concat(to.be.concatenated)	
			concatenated.graph.matrix <- concat.result[[1]]
			components.list <- concat.result[[2]]
			if ( sum(to.be.concatenated)  == sum(concatenated.graph.matrix) ) break		
			to.be.concatenated <- concatenated.graph.matrix
		}
	}

	## building components 
	components <- list()
	number.of.components <- 0
	for(i in 1:length(components.list))
		if(components.list[[i]][1] != "vanished"){			
			number.of.components <- number.of.components +1
			components[[number.of.components]] <- components.list[[i]]
		}

	# points in the components
	components.points <-list()					
	component.of <- rep(NA, times=dim(full)[1] )							
	number.of.points.in_component <-c()			

	for (i in 1:length(components)){			
		components.points[[i]] <- components[[i]][1]
		for (j in 1:length(components[[i]])){	
			jth.cluster <- cluster[[ components[[i]][j] ]] 
			components.points[[i]] <- 	union( components.points[[i]] , jth.cluster)
			component.of[ components.points[[i]] ]<- i 
		}
		number.of.points.in_component[i] <-length(components.points[[i]])
	}


	# Reporting number of components
	if(talk) message("Number of components after connecting:")
	if(talk) message(paste("------------------------------------->", length(components)))

	################################################################# E N D ########################################################################
	
    return(list(label=component.of,clusters.graph=max.conductance))
    	#clusters.graph is the graph of clusters before connnecting is done.
}















