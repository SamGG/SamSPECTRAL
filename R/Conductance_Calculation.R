###########################################################################
# Here, we assume we have some commuties. Each community has a representative and a list of members.
# We intend to calcuate an estimation for  conductance between each two communites c1 and c2 according to th efollowing formula:
#
#	conductance(c1,c2) = sum( conductance(m1,m2) ), where the sumation is over all members of c1 and c2.
#
#	We can do this, because we have list of members.	
#
# Data Structure:
# each community = (the representative, set of members) 
#
#
# So in R, we have implement the data structure as follows:
# 
#  community = list()	# This is the list of communities.
#
# c1 <- list(repres = 346, members = c(237,8769,23))  		
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
# for (i in 1:length(community))  if(talk) message(community[[i]])
#######################################################



### Computing sigma...
	# (1) It takes too much time, so we estimate it.
	# (2) We may adjust locally it according to the density.
	# (3) We may learn it from some clustered samples. 
##########################################################################




Conductance_Calculation <- function(full, normal.sigma, space.length, society, precision, talk=TRUE, beta=4)
{
	t1<-Sys.time()
	if(talk) message(t1)
	########################### S T A R T ###########################

	# setting parameter sigma, determinds how the conductance decays by distance in th formula: 
	#point to point condustance = exp(- sigma * distance^2)
	sigma= normal.sigma / (space.length)^2	# since the formula is exp(-sigma*distance^2)
	if(talk) message(paste("sigma = ", sigma))


	# Loading society
	nbhood <- society$nbhood
	repres.indeces <-society$representatives
	community <- society$communities
	num.of.cummnities <- length(community)

	### Initialization
	n <- dim(full)[1]			# Total number of points.
	dimention <-dim(full)[2]		# The number of measured parameters.


	### Computing conductance between communites:
	#dyn.load("src/SamSPECTRAL.so")  # C function, packaging
	#conductance.matrix <-.Call("conductance_computation",society,full, sigma)
	conductance.matrix <-.Call("conductance_computation",society,full, sigma, beta, PACKAGE = "SamSPECTRAL") #package
	conductance.matrix <- round.conductance.matrix <- round(conductance.matrix,precision)	# solving precision problem.
	conductance <- list(conductance.matrix=conductance.matrix, sigma = sigma)

	### Outputing:
	#outfile = paste("communities/",num.of.cummnities, " communities with resistance,sigma", sigma, sep="")
	#save(conductance,file= paste(outfile,".Cnd",sep=""))       
	#pdf(file=paste(outfile,".pdf",sep=""))
	#plot(full[repres.indeces,])
	#dev.off()


	########################### E N D ###########################
	if(talk) message(Sys.time()-t1)

	conductance
}
