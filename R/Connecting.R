

Connecting <- function(full, society, conductance, number.of.clusters, labels.for_num.of.clusters, separation.factor, talk=TRUE)
{
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
# data structure;
# cluster[[]]	:	the list of indices for each cluster.
#					e.g. cluster[[3]] is the vector of indices
#					in cluster 3.
#
# components	:	a list, each element is a vector that includes 
#					the indices of clusters in each component 
#
# separation.factor:
# 				# This is a positive factor typically less than 1.
				# The larger it is, it is more probable 
				# that we divide the data into two or more components.
				# If it is 0, we may not devide at all. 
				# It it is +Inf, we just return the input clusters without 
				# connecting any together.	
###########################################################################

	# Compute indices of points in cluster j from the output of civilized spcetral clustering
	number.of.clusters.from.spectral <- number.of.clusters			
		# We get this many clusters from spectral clustering.
	### The clusters
	labels <-labels.for_num.of.clusters[[number.of.clusters.from.spectral]]
	repres.cluster	<-list()
	cluster<- list()
		# cluster[[3]] should be the indices of points in cluter 3.
		# If this is not given, we compute it from the out put of 
		# civilized spcetral clustering as follows.
	
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

	### ###connection graph
		#This matrix represents a directed graph. There is a node for each cluster.
		# If a cluster i tends to attach to cluster j, then we put a directed edge from i to j.
		connection <- max.conductance

	### ### concatenatation
	# visualized in picture 4
	# If a node of the graph points to just one other node, we omit the first node.
	# This is done technically by deleting its column and adding its row to the row for the 
	# host node. We do this in order to remember which nodes had pointed to the deleted node.

	# Components
	components.list<-list()
	for (i in 1:dim(connection)[2])
		components.list[[i]] <-c(i)				
		# For now, each component contains exactly one node.		

	### concat functions
	# We will delete nodes from this graph step by step by the following function.
	# input.graph.matrix has an extra row at the bottom 
	# which contains the indices of the vertices in the original graph of clusters.

	# If both of two nodes have tendency to combine, we combine them.
	concat.two.sided <- function(input.graph.matrix){

		# input check
			# If input matrix is 2*1 then it will not be  of class "matrix".
			# Because the last row just keeps the indices, 
			# it means there is only one node that can not be concatinated more.
			if(class(input.graph.matrix) != "matrix"){ 	
				result <- list (input.graph.matrix, components.list)
				return(result)	
			}#End if(class(input.graph.matrix) != "matrix").
		
		output.graph.matrix <- input.graph.matrix
		# We would like to devide the label of outgoing edges from a node 
		# by the label of the edge from the node to itself.

		within.connection <- diag(output.graph.matrix)
		# Describes the connections inside the clusters.


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
			# This row just keeps the components.
		for (i in 1:dim(output.graph.matrix)[2] ){	# for all comlumns,
			for (j in 1:dim(output.graph.matrix)[2] ){ # for all comlumns,
				if ( (normalized.output.graph.matrix[i,j] >=1)  && (normalized.output.graph.matrix[j,i] >=1) && (i!=j) ) {	# high tendency
					#updating components
					the.two.components <- c(last.row[i], last.row[j])
					host <- max(the.two.components) 				
					guest <- min(the.two.components) 				
					components.list[[host]] <- union( components.list[[host]], components.list[[guest]])
					components.list[[guest]] <- "vanished"
					#message("host"); message(host)
					#message("guest"); message(guest)

					# combine node i to node j.
					output.graph.matrix[j,] <- pmax(output.graph.matrix[j,] , output.graph.matrix[i,])			# combining the rows.
					output.graph.matrix[,j] <- pmax(output.graph.matrix[,j] , output.graph.matrix[,i])			# combining the columns.
					output.graph.matrix <- output.graph.matrix[-i,]								
					output.graph.matrix <- output.graph.matrix[,-i]								
					# concatinating of just one pair is enough for each call of this function.
					result <- list (output.graph.matrix, components.list)
					return(result)
				}
			}
		}
		result <- list (output.graph.matrix, components.list)
		return(result)
	}###End of concat.two.sided.


	
	concat <- function(input.graph.matrix){
    	## input check
    	if(class(input.graph.matrix) != "matrix"){	
	    	# If input matrix is 2*1 then it will not be  of class "matrix".
				# Because the last row just keeps the indices, 
				# it means there is only one node that can not be concatinated more.
    		result <- list (input.graph.matrix, components.list)
    		return(result)	
    	}
    
    	output.graph.matrix <- input.graph.matrix
    	number.of.columns <- dim(input.graph.matrix)[1]-1
    	last.row <-input.graph.matrix[dim(input.graph.matrix)[1],]
	    	# This row just keeps the indices.	
    
    	for (i in 1:dim(output.graph.matrix)[2] ){	# for all comlumns,
    		#message(">>>>i : ");message(i)
			#message(last.row[i])
			
			### To figure out if this node should be omitted

			# We would like to devide the label of outgoing edges from a node 
			# by the label of the edge from the node to itself.
    		matrix.without.last.row <- output.graph.matrix[1:number.of.columns,]	
	    		# The last row contains the index of the node 
				# which is not good for us here.	

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
				# This may include the i, itself.
    																					
    		out.going.edges <- out.going.edges[ which( !(out.going.edges ==i)) ]		
	    		# All out.going.edges excluding the point i , itself.
    
    		## strategy
	   			# We could have two different approaches here:
				#	1) if there is more than one highly connected node, stop combining for now.
				#	2) combine with the node with maximum connection
	
			# first strategy:
				 #if (length(out.going.edges) ==1)	# This node is connected to exactly one node other than itself.
					#host.column <- out.going.edges[ which(out.going.edges >=1)  ]			
						# This is the column of the node that i is connected to.			
						# Then, we should concatinate i to the appropriate node.
						
			# second strategty:		
    		if (length(out.going.edges) >=1){
    			connection.to.others <- normalized.output.graph.matrix[,i]			
	    			# We would like to pick the most connected node other 
					# than the guest itself.
    			connection.to.others[i] <- -Inf
				### Once the host.column is selected, the rest is the same.	###
    			host.column <- which.max(connection.to.others)
    			host.row <- host.column	
    			host <- last.row[host.column]
		
				# gust will be omitted and its compoent will be added to host.
				# We will do omitting by putting 0 in its row and column.
				# so, if an entry in the lst row is 0, it means that the corespondig compunent is 
    			guest.column <- i
    			guest <- last.row[i]		
				#message("guest: "); message(guest);message( "host"); message(host)
				
				# adding the list of component of i to the list of host
    			components.list[[host]] <- union(components.list[[host]] , components.list[[guest]])
    			components.list[[guest]] <- "vanished"
				#message("components.list:     >>");message(components.list)

				### omiting the node guest.
    			ith.row <- output.graph.matrix[i,]
    			ith.column <- output.graph.matrix[,i]
    			
				# combining the omitted to the host.
    			output.graph.matrix[host.column,] <- pmax(output.graph.matrix[host.column,] , ith.row )		# combining ith.row and host.column
    			output.graph.matrix[,host.row] <- pmax(output.graph.matrix[,host.row] , ith.column )		# combining ith.column and host.row
				# last row is an exception for which we should not take max.				
    			output.graph.matrix[dim(output.graph.matrix)[1] , host.column] <- host		
    			output.graph.matrix	<-output.graph.matrix[-i,]			# deletes the i'th row.
    			output.graph.matrix	<-output.graph.matrix[,-i]			# deletes the i'th column.
       			result <- list (output.graph.matrix, components.list)	
    			return(result)
    		}#End if (length(out.going.edges) >=1).
        }#End for (i in 1:dim(output.graph.matrix)[2] ).
	result <- list (output.graph.matrix, components.list)
	return(result)
    }#End concat <- function(input.graph.matrix).



	connection.with.indices	<- rbind(connection, 1:dim(connection)[2])	
		# The extra row at the bottom contains the indices of 
		# the vertices in the original graph of clusters.
		
		
	# We concat edges untill there is no change.
	to.be.concatenated <- connection.with.indices	
	#Main loop.
	repeat {	
		#two sided
		repeat {
			concat.result <-concat.two.sided(to.be.concatenated)	
			concatenated.graph.matrix <- concat.result[[1]]
			components.list <- concat.result[[2]]
			if ( sum(to.be.concatenated)  == sum(concatenated.graph.matrix) ) break
				# means we don't need to continue, since the matrix hasn't change.
			to.be.concatenated <- concatenated.graph.matrix
		}#End of two sided.
		
		### separation.factor	
		concat.result <-concat(to.be.concatenated)	
		concatenated.graph.matrix <- concat.result[[1]]
		components.list <- concat.result[[2]]
		if ( sum(to.be.concatenated)  == sum(concatenated.graph.matrix) ) break		
			# means we don't need to continue, since the matrix hasn't change.
		to.be.concatenated <- concatenated.graph.matrix
	}

	## building components 
	# which is a list, each element is a vector that includes the indices of clusters in that component 
	components <- list()
	number.of.components <- 0
	for(i in 1:length(components.list))
		if(components.list[[i]][1] != "vanished"){	# this is a none-empty component.			
			number.of.components <- number.of.components +1
			components[[number.of.components]] <- components.list[[i]]
		}
	# points in the components
	components.points <-list()		# Each element is a vector of indices of points in a component.			
	component.of <- rep(NA, times=dim(full)[1] )	# For each point, we assign it to a component.							
	number.of.points.in_component <-c()				# How many point are there in each components?
	for (i in 1:length(components)){	# for each component,				
		components.points[[i]] <- components[[i]][1]
		for (j in 1:length(components[[i]])){	#each cluster in the component.
			jth.cluster <- cluster[[ components[[i]][j] ]] 
			components.points[[i]] <- 	union( components.points[[i]] , jth.cluster)
			component.of[ components.points[[i]] ]<- i 
		}
		number.of.points.in_component[i] <-length(components.points[[i]])
	}
	# plot(as.data.frame(full),pch='.',col=component.of)

	# Reporting number of components
	if(talk) message("Number of components after connecting:")
	if(talk) message(paste("------------------------------------->", length(components)))

	################################################################# E N D ########################################################################
    return(list(label=component.of,clusters.graph=max.conductance))
    	#clusters.graph is the graph of clusters before connnecting is done.
}















