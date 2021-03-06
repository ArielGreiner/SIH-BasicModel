function(dispersal=0.001,species=9,patches=30){
#Constants####
N<- matrix(10,ncol=species,nrow=patches) # Community x Species abundance matrix
R<-rep(10*(species/10),patches) #Initial resources

rInput<-150 #resource input
rLoss<-10 #resource loss 
eff<-0.2 #conversion efficiency
mort<-0.2 #mortality
Ext<- 0.1 #extinction Threshold

ePeriod<-40000 #period of env sinusoidal fluctuations
eAMP<-1 #amplitude of envrionment sinusoidal fluctuations

Tmax<-140000 #number of time steps in Sim
DT<- 0.08 # % size of discrete "time steps" - this is the Euler value

#vectors####
eOptimum<-1-seq(0,eAMP, by=eAMP/(species-1)) #species environmental optima

#dispersal conditions####
dispersal_matrix<-matrix(1/(patches-1),patches,patches)
diag(dispersal_matrix)<-0

calc.immigration <- function(N,a,dispersal_matrix) dispersal_matrix%*%N*rep(a,each=patches)

Prod<-matrix(NA,species*patches,40000)
Abund<-Prod

for(TS in 1:Tmax){
  Immigrants<-calc.immigration(N,dispersal,dispersal_matrix)
  envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(1:patches)*2*pi/patches)+1)
  consume <- 0.1*(1.5-abs(sapply(eOptimum,'-',envt.v)))
  Nt <- N*(1+DT*(eff*R*consume - dispersal - mort)) + DT*Immigrants #abundance step
  Rt <- DT*rInput+R*(1-DT*(rLoss + rowSums(consume*N))) #resource step    
  N <- Nt * (Nt>Ext) # set to 0 if below extinction threshold
  R <- Rt
  
  if(TS>=100000){
    Prod[,(TS-100000)] <- c(t(eff*consume*R*N))
    Abund[,(TS-100000)] <- c(t(N))
  }
} 

Prod<-array(t(Prod),dim=c(40000,species,patches))
Prod<-Prod[seq(1,40000,100),,]

Abund<-array(t(Abund),dim=c(40000,species,patches))
Abund<-Abund[seq(1,40000,100),,]
return(list(Prod,Abund))
}


