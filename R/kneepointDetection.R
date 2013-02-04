########################################## A FUNCTIONS used in Civilized_Spectral_Clustering():##################
kneepointDetection <- function(vect, PlotFlag=FALSE){
	# OrthagonalResiduals=FALSE
	#Finds the change point using piece wise linear regressions. See the Rd document for details. 
	  n<-length(vect);
	  Vect<-vect;


	  
	  a<-as.data.frame(cbind(1:n, Vect[1:n]));
	  l<-lm(a[,2]~a[,1], data=a);
	  MinError=100000000;
	  MinIndex=1;
	  for (i in 2:(n-2)){
	  	#first line
	  	# We fix the first line to be always horizontal (Y=1).
		a<-as.data.frame(cbind(1:i, Vect[1:i]))
                l1<-lm(a[,2]~a[,1], data=a)
		#error 1
		e1<-sum(abs(1-a[,2]))

		#second line  
		a<-as.data.frame(cbind((i+1):n, Vect[(i+1):n]))
		l<-lm(a[,2]~a[,1], data=a)
		l2<-l
		#error 2
		e2<-sum(abs(l$residuals))
		Error=e1+e2;
		if (MinError>Error){
		  MinError=Error;
		  MinIndex=i;
		}

		if (PlotFlag){
		  if (!.Platform$OS.type=='unix')
		    stop('The plotting feature only works on unix systems')
		  if (!file.exists('tmpfigs'))
		    system('mkdir tmpfigs');      
		  png(sprintf('tmpfigs/%.3d.png', i))
		  a<-as.data.frame(cbind(1:length(Vect), Vect))
		  plot(a, xlim=c(0,n), ylim=c(0,max(a)), axes=FALSE, xlab='Iteration', ylab='Distance');
		  par(new=TRUE);
		  plot(a[MinIndex,], pch=19, col='green', xlim=c(0,n), ylim=c(0,max(a)), axes=FALSE, xlab='Iteration', ylab='Distance')
		  par(new=TRUE);
		  plot(a[i,], pch=19, col='red', xlim=c(0,n), ylim=c(0,max(a)), axes=FALSE, xlab='Iteration', ylab='Distance')
		  axis(1);axis(2);
		  title(sprintf('SSD: %.3f',Error))
		  title(sub='Red: Current Change Point       Green: Best Change Point');
		  abline(l1)
		  abline(l2)
		  dev.off();
		}
	  }
	  if (PlotFlag){
		system('convert -delay 40 tmpfigs/* animation.gif')
	  }

	  i<-MinIndex
	  a<-as.data.frame(cbind(1:i, Vect[1:i]))
	  l1<-lm(a[,2]~a[,1], data=a)
	  Error=mean(abs(l$residuals));
	  a<-as.data.frame(cbind((i+1):n, Vect[(i+1):n]))
	  l2<-lm(a[,2]~a[,1], data=a)
	  a<-as.data.frame(cbind(1:n, Vect))

	  return(list(MinIndex=MinIndex, l1=l1 ,l2=l2));
}# End kneepointDetection <-function(vect, PlotFlag=FALSE)	

















	
