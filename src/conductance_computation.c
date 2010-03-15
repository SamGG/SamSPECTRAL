#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


static double betaa = 4;

// MACROs; Faster way than calling a function: 
// The index (i,j) in matrix M having n rows. i and j start from 0. 
#define INDEX(M,n,i,j) (M)[(i) + (j)*(n)]

#define RESISTANCE(coordinates, pNum, dimention, ind1, ind2, sigma) 	exp((2)*(sigma)* (Distance_power2)(coordinates, pNum, dimention, ind1, ind2)) 
	// That is: RESISTANCE =1/ CONDUCTANCE, so there is no negative sign in exp.

// Other option, without exp.
//#define RESISTANCE(coordinates, pNum, dimention, ind1, ind2, sigma) 	sqrt(Distance_power2(coordinates, pNum, dimention, ind1, ind2))




// Surprisingly, it is not built-in in C!
inline double max(double a, double b){
	if (a>b) return a;
	else return b;
}

// This is the square of distance between two points.
// c=coordinates, d=dimention
// inline, a tiny function that runs in register memory very fast.
inline double Distance_power2(double *c, int pNum, int d, int ind1, int ind2)
{
	int i;
	double sum;
	
	sum=0;
	// Remind: X1(i) = INDEX(c,pNum,ind1-1,i)
	for (i=0; i< d; i++)
		sum+= ( INDEX(c,pNum,ind1-1,i) - INDEX(c,pNum,ind2-1,i) ) 
				*  ( INDEX(c,pNum,ind1-1,i) - INDEX(c,pNum,ind2-1,i) );	
		// The above is: ( X1(i) - X2(i) ) ^2

	return sum; 
}







// The following function computes the conductance between the two communities.
//Currently,we are just adding the conductance of all of the edges between them.
double compute_sum_of_edges(double *coordinates, int pointNum, int dimention, 
			double nbhood, SEXP community, int ind1, int ind2, double sigma)
{
	int i, j, length1, length2;
	double kAprrox, sum, result;
	SEXP community_ind1, community_ind2;
	int *member_ind1, *member_ind2;								
	// These arrays will contain the indeces of memebers of the two communities.

	int repres_ind1 ,repres_ind2;


																
	// These two lists are two communites. 
	//In R, each is a list containing repres and members.
	PROTECT( community_ind1 = VECTOR_ELT(community, ind1));				
	PROTECT( community_ind2 = VECTOR_ELT(community, ind2));

	// Each of these will be the index of the represantative
	// of one of the two communities.
	repres_ind1 = INTEGER(VECTOR_ELT(community_ind1, 0))[0];	
	repres_ind2 = INTEGER(VECTOR_ELT(community_ind2, 0))[0];

	// These arrays will contain the indeces of memebers of the two communities.
	member_ind1 = INTEGER(VECTOR_ELT(community_ind1, 1));		
	member_ind2 = INTEGER(VECTOR_ELT(community_ind2, 1));
	length1 = length(VECTOR_ELT(community_ind1, 1));
	length2 = length(VECTOR_ELT(community_ind2, 1));


	

 	if ( sqrt(Distance_power2(coordinates, pointNum, dimention, 
			repres_ind1, repres_ind2)) > betaa * nbhood){						
		//This is NOT square distance,what we used for building the communities.
		// According to picture 3:	

		kAprrox = ((double)length1 * length2 ) / 
					RESISTANCE(coordinates, pointNum, dimention, 
						repres_ind1, repres_ind2, sigma);
		
		/*//To test:
		if (ind1 ==2 -1 && ind2==1322 -1){
			Rprintf("\n length1 = %d\n", length1);
			Rprintf("\n length2 = %d\n", length2);
			Rprintf("\n kAprrox = %e\n", kAprrox);
			Rprintf("\n 1/2 = %e\n", ((double)1/2) );			
		}*/

		result = kAprrox;
	}
		
	else{// compute sum;										
	// Now, we have to consider all the members 
	//because the communities are relatively close.

																
		//Picture 3 explains the following formulas 
		//for calculating the normal cut.
		sum =0;
		for (i=0; i<length1; i++)								
		// To add the conductance of all edges between the two communities. 
		// Loot at picture 2 for the mening of the formula 
		//used to compute conductance.
			for (j=0; j<length2; j++){

				sum += 1/RESISTANCE(coordinates, pointNum, dimention, 
							member_ind1[i], member_ind2[j], sigma);	
						// coordinates is the full matrix of all points of data.

				// To debug,
				//Rprintf("\n in1+1 = %d \n",ind1+1);
				//Rprintf("\n in2+1 = %d \n",ind2+1);
				//Rprintf("\n member_ind1[i] = %d \n",member_ind1[i]);
				//Rprintf("\n member_ind2[j] = %d \n",member_ind2[j]);
			}
 
		result= sum;
	}
	UNPROTECT(2);
	return result;			
}// End of compute_sum_of_edges.






// To compute the contuctance between all communites.
// "coordinates" contains the coordinates of all points iin the data.
SEXP conductance_computation(SEXP society, SEXP coordinates, SEXP sigmaVal)
{
	SEXP community, returnVal, coordDim;	
	int communityNum, memberLen, i, j, pointNum, dimention;
	double * coords, nbhood, di, dj, sigma, 
		sumOfEdges, weightBetweenCommunities; //epsilon
	int *density, *repres_ind;

	/* VARIATIONS: */
	//Rprintf("\n Strategy: sum of edges using exp formula.\n	");
	//Rprintf("       betaa = %f\n",betaa);	

	// To set the type from R to C.
	PROTECT( society = coerceVector(society, VECSXP));	
	PROTECT( coordinates = coerceVector(coordinates, REALSXP));	
	PROTECT( sigmaVal = coerceVector(sigmaVal, REALSXP));
	//// Computing the values of the parameters.
	sigma = REAL(sigmaVal)[0];
	PROTECT( coordDim = getAttrib(coordinates, R_DimSymbol));
    pointNum = INTEGER(coordDim)[0];
    dimention = INTEGER(coordDim)[1];
	coords = REAL(coordinates);
	// Reading information from society
	nbhood = REAL(VECTOR_ELT(society, 0))[0];		// nbhood

	// The indeces 	of the representatives.
	repres_ind = INTEGER(VECTOR_ELT(society, 1));	
	// society is a list which has the communities as its third argumant.
	PROTECT( community = VECTOR_ELT(society, 2));	
	
	// This is the total number of cummnities.
	communityNum = length(community);				

	// An estimate for the density of the represantative 
	//which is equal to the population of th ecommunity.
	density = INTEGER(VECTOR_ELT(society, 3));
	// Almost the minimum distance between the points in the data.
	//epsilon = REAL(VECTOR_ELT(society, 4))[0]; not needed for now.		



	// Computing the conductances:
	/************************************************************************/

	// Allocate memory for returning value to R which is kept in a matrix..
    PROTECT(returnVal = allocMatrix(REALSXP, communityNum, communityNum));	
				
	// Used and explained  in function: CommunityConductance().

	// For any pair of communities,
	// we should compute the normalized cut between them.
	for (i=0; i< communityNum; i++)
    	for (j=0; j< communityNum; j++)
    		if (i==j)
	    		REAL(returnVal)[i+j*communityNum] = 0;	// Not important. 
	    	else{

				// We are using the sum of edge conductance 
				//as an estimate for community conductance.	
				sumOfEdges= compute_sum_of_edges(coords, pointNum, dimention, 
								nbhood, community, i, j, sigma);	    		
				weightBetweenCommunities = sumOfEdges;	
				// This is JUST the enumarator of the normal cut.
				
				// Probably produced by divide by 0.
				if (isnan(weightBetweenCommunities)==TRUE){ 
					weightBetweenCommunities = 0;
				}
				// To debug:
					//Rprintf("\n i= %d - j=%d \n",i+1,j+1);
					//Rprintf("\n%f\n", sum_of_edges);
					//Rprintf("\n%f\n", weightBetweenCommunities);
					
	
				REAL(returnVal)[i+j*communityNum] =	weightBetweenCommunities;
			}

  	UNPROTECT(6);
	return returnVal;	
}










