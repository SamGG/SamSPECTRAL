###	INTENTION:
###########################################################################
# Here, we get a matrix full, and omit prune it to get the matrix selected.
# We would like to try this idea:
#
# 	- 1) nhbood: We compute h=nbhood based on the effective volume of the data and the final number of points we want to select.
#	- 2) density: We compute the density of all points in the data as the number of points which are in distance less than h.
#	- 3) pickup randomly: with probablity 1/density (We hope to get a space uniform sampling.)
#	Now, we are done with the sampling part. Next, we try to extract and keep information about the cummunities and represantatives.
#	Later, we'll use this information for assigning weight to the edges.
#
#	- register: For each cummunity, assign the members. This must be a PARTITIONING, each point exactly one cummunity.
#		    So, after regitering the members of each community we'll use [-index] to omit the registered points 
#			from the list.
#
#	Later on, our plan is to compute the conductance between two communities as the some of the conductance of their members.
#	We can do this, because we have list of members.	
#
# Data Structure:
# each community = (the representative, set of members) 
#
#
# So in R, we can implement the data structure as follows:
# 
#  community = list()	# This is the list of communities.
#
#  c1 <- list(repres = 346, members = c(237,8769,23))  		
# This is a typical cummunity where 346 is the index of the representative and the indeces of 
# the members is stored in the vector:  c(237,8769,23)
# 
# ...c1$repres
# ...c1$members
#
# c1 might be the 17th community:
# cummunity[[17]] <- c1
#
# To see all cummnities:
# for (i in 1:length(community))  message(community[[i]])
#
# Good luck with implementing this idea!
###########################################################################

Building_Communities <- function(full, m=3000, space.length=1, community.weakness.threshold=1,talk=TRUE, do.sampling=TRUE)
{
	

	t1<-Sys.time()
	if(talk) message(t1)
	########################### S T A R T ##########################


	### Initialization
	n <- dim(full)[1]			# number of points.
	dimention <-dim(full)[2]	# the number of measured parameters.

	# Call functions
	#source("Compute_Close_Poits.R")
	################################ packaging
    	# This function gets the index of a points and compute and return all close points to it which is basically the indeces of all points closer than h.
        # We use the metric: MAX (|X_1 -X_2| , |Y_1-Y_2|)
        Compute_Close_Points <- function(point.ind, h, n, dimention, epsilon){											
	
    	    point <-full[point.ind,]
    	    
    	    # computing the distance
    	    total.number <-	n							
    	    # This is the total number of points in the space.			
        
    	    pmatrix<- t(matrix(point, dimention, total.number)) 			
    	    # A matrix with height equal to total number of points which has copies of p on its     rows.
    	    difference <- abs(full - pmatrix)
        
            #packaging
	        #square.distance <- .Call("maximum_of_rows",difference)			# This vector, contains the distance between the point and all other points.
	        square.distance <- .Call("maximum_of_rows",difference, PACKAGE = "SamSPECTRAL")			
	        close.points <- which(square.distance < h )
        
			# To compute epsilon which is almost the minimum distance between points in the data,
				square.distance.without.THE.point <- square.distance[-point.ind]
				if (min(square.distance.without.THE.point) != 0){ #We no not consider points with the same coordinates.
    				epsilon <- min(union(epsilon,square.distance.without.THE.point))	
		        }
	        result <- list(close.points= close.points, epsilon = epsilon)											
    	    return(result)
        }


	################################ End of ("Compute_Close_Poits.R")

	###### ###### ###### Sampling...

	if(talk) message("We assume that the points are in a totaly random order.")
	if(talk) message(" Any position related order may cause problem.")

	### Computing h=nbhood
	nbhood <- (space.length/ (m ^ (1/dimention)) )/2 	# Points closer than this threshould are considered as neighbours, This derives from the formula: 
								# average neighbours of a pixel = n/m = (2 * nbhood)^d * (n/l^d) where l is the length of the universe
								# and d is its dimention. 

	if(talk) message(paste("neighbourhood = ", nbhood))


	repeat{									# If nbhood is not set properly, we may build less than enough number of communities.
											# So we repeat building with smaller nbhood to get more communities.
											# We will break this loop at the end if we make sure that enough number of communities
											# are built (repres.num > m/2) or if all points are selected(m=n).


								
		### Initializtion before building the communities
		epsilon <- Inf					# This will be almost the minimum distance between the points in the data.

		registered <- rep(FALSE, times = n)		# No point is registered yet.

		community <- list()				# List of all communities.

		repres.indeces <-c()				# This is the indeces of all representatives.

		repres.num <- 0					# Number of representatives so far which is equal to the number of communities.

		close.points <- c()				# For each point p_i, this is a vector that contains the indeces of close points to p_i.

		repres.density <- c()				# Defined to be the number of close points to the represantative in each community.



		### To compute the densities and picking up representatives. Then building communities:
		for (i in 1:n){

			if (! registered[i]){			# If this point is not already registered, it establishes a community.
				repres.num <- repres.num +1
				repres.indeces[repres.num] <- i
		
				if (do.sampling){
					ccp <- Compute_Close_Points(repres.indeces[repres.num],nbhood, n, dimention, epsilon)
					close.points <- ccp$close.points
					epsilon <- ccp$epsilon
					repres.density[repres.num] <- length(close.points)		
				} else {	# No sampling will be done, every poit will make a cummunity.
					close.points <- i
					epsilon <- "not_sampled"
					repres.density[repres.num] <- as.integer(1)
				}#End if (do.sampling).
				# Registering members; All of the close points to this representative who have not registered to a community before
				# must register as members of this community.
		
				not.registered.neighbours <- close.points[which(registered[close.points]== FALSE)]	# All neighbours which have not alreay registered.

				# Building the community
				community[[repres.num]] <- list (repres = i, members = not.registered.neighbours)
		
				# Updating the registered points.
				registered[not.registered.neighbours] <- TRUE		
			}
		}
		if(talk) message(repres.num)
		if(talk) message(paste("^--------^ communities have been build up."))

		### breaking the repeat loop.bbbbbbbbbbbbbbbbbbbb
		if( ( (repres.num <= m*(1.1) ) & (repres.num > 0.95*(m/2) ))  | repres.num >= n | !do.sampling)	break   			
    		                                                    # We break this loop only if we make sure that enough number of communities
																# have been built (repres.num > m/2).
																# The condision (repres.num < m) is to make sure the number of communities is not so high.
		
		
	
		k<- (m/repres.num)^(1/dimention)
		
		if (nbhood < epsilon) break                             # If this happens, then we will not get more communities by decreasing nbhood!
		#if (nbhood < nbhood / k) break                        # This happens when nbhood is not reduced by dividing by k!
	
    		                                                    # These situation may happen if many points have the same coordinates,
    		                                                    # then nbhood will converge to the persision level of the machine (10^(-324))
    		                                                    
		nbhood <- nbhood / k									# <= We can reduce nbhood to go to higher resolution. 
																# By dividing by 2^(1/dimention), we expect the number of cummunities to double
																# as the volumne of each community will be devided by two..
																
		if(talk) message(paste("neighbourhood is changed to: ", nbhood))
		#End breaking.bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb

	}#End repeat.	



	# Discarding weak communities
	ith.community <- 1
	while(ith.community <= length(community)){		#from all the ecommunities,
		if(length( community[[ith.community]]$members ) < community.weakness.threshold ){	#This is a low populated community.
			# omitting the community.
			community[[ith.community]] <- NULL
			repres.indeces <- repres.indeces[-ith.community]
			repres.density <- repres.density[-ith.community]
		}
		else 
			ith.community <- ith.community +1
	}
	if(talk) message("But just"); if(talk) message(paste("......>>",length(community), "are left as high populated ones."))

	if(talk) message(paste("epsilon= ",epsilon))
	society <- list(nbhood =nbhood, representatives=repres.indeces, communities= community, densities= repres.density, epsilon= epsilon)


	########################### E N D ###########################
	if(talk) message(Sys.time()-t1)
	
	society
}
