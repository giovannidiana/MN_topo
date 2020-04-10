
# This function calculates the drift in z returning a list of two array corresponding to the first n samples (dropping the first "skip" samples) in the variable "draw". You can access those lists to make histograms like this:
# > t = drift_z(1000,100) 
# > hist(t$drift_z_wt)
# or, to calculate the difference between them you can do
# > t_diff=t$drift_z_wt-t$drift_z_mu
# and plot the histogram of the difference. You can then do statistics exactly in the same way you did for the KL but now you can get negative values, which then allows you to make statistical statements about the positivity or negativity of the shift.
# note: you want to choose n to be equal to the total number of draws that you had in your simulated draws.
drift_z <- function(n=500,skip=100){
	drift_wt <- rep(NA,n-skip)
	drift_mu <- rep(NA,n-skip)
	zval <- rep(0:14,4)
	counter=0
	for(i in (skip+1):n){
		counter=counter+1
		drift_wt[counter]=sum(zval*draw$alpha[1:60,1,i,1]/sum(draw$alpha[1:60,1,i,1]))
		drift_mu[counter]=sum(zval*draw$alpha[1:60,2,i,1]/sum(draw$alpha[1:60,2,i,1]))
	}
	return(list(drift_z_wt=drift_wt,drift_z_mu=drift_mu))
}

# Same as before but calculated with the lateral dimension. 
# note the different definition used for lval compared to zval in the previous function
drift_lateral <- function(n=1000,skip=100){
	drift_wt <- rep(NA,n-skip)
	drift_mu <- rep(NA,n-skip)
	lval <- rep(1:4,rep(15,4))
	counter=0;
	for(i in (skip+1):n){
		counter=counter+1
		drift_wt[counter]=sum(lval*draw$alpha[1:60,1,i,1]/sum(draw$alpha[1:60,1,i,1]))
		drift_mu[counter]=sum(lval*draw$alpha[1:60,2,i,1]/sum(draw$alpha[1:60,2,i,1]))
	}
	return(list(drift_l_wt=drift_wt,drift_l_mu=drift_mu))
}

# This function plots the two shifts (lateral and over z) as a scatter plot 
plot_drifts <- function(){
	tz=drift_z(1000,100)
	tx=drift_lateral(1000,100)
	plot(tz$drift_z_wt,tx$drift_l_wt,xlim=c(5,7.5),ylim=c(1.8,2.5),pch=19,cex=0.2)
	points(tz$drift_z_mu,tx$drift_l_mu,col="red",pch=19,cex=0.2)
	legend("topright",fill=c("black","red"),legend=c("WT","mutant"))
}

