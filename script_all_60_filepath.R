model="model {

    for( i in 1:60) {
		alpha[i,1] ~ dunif(0.5,50)
		alpha[i,2] ~ dunif(0.5,50)
	}

	for( s in 1:41 ){
		p[s,1:60] ~ ddirch(alpha[1:60,gbys[s]])
	}

	for( i in 1:N ){
		ind[i] ~ dcat(p[S[i],1:60])
	} 
  
}"

tab<-read.table("data_IO.csv",header=T)
gbys<-rep(NA,41)
for(i in 1:41) gbys[i] = tab$G[tab$S==i][1]
data=list(ind=tab$Z+15*(tab$X-1)+1,S=tab$S,gbys=gbys,N=nrow(tab));
varnames=c("alpha","p")
burn_in=10;
steps=20000;
thin=1;

library(rjags)
fileConn=file("~/Documents/Giovanni_bayesian_inference/model.tmp")
writeLines(model,fileConn);
close(fileConn)

m=jags.model(file="~/Documents/Giovanni_bayesian_inference/model.tmp",data=data);
update(m,burn_in)
draw=jags.samples(m,steps,thin=thin,variable.names=varnames)

KL <- function(n=2000){
  kl <- rep(NA,n)
  for(i in 1:n){
    p1=draw$alpha[1:60,1,i,1]/sum(draw$alpha[1:60,1,i,1])
    p2=draw$alpha[1:60,2,i,1]/sum(draw$alpha[1:60,2,i,1])
    kl[i] = sum(p1 * log(p1/p2) )
  }
  return(kl)
}



