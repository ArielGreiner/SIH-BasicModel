#Initial

#Note: rep(thingtoreplicate,numberoftimesitreplicatesit)

#Values
nspecies = 9
npatches = 30
nfunctions = 7
continuity = 0.08
e = 0.2
m = 0.2
I = 150
l = 10
binom_p = 0.5
timesteps = 140000
Period = 40000
dispersalrate = 0.01 
sumforR = 0
sumforN = 0
#Random iterands...
x=0
y=0

#Tip: Use NA instead of NaN, NaN can cause some problems (what comes up when there's a mistake), seq() is super useful, NaN = not a number

#initialize abundance matrix
N<-matrix(10,ncol = nspecies,nrow = npatches)

#initialize resource vector
R<-rep(9,npatches)

#initialize environment matrix 
#E <- matrix(NA, nrow = npatches, ncol = timesteps)
#E[,1] <-seq(0,1,length=npatches)
E <-seq(0,1,length=npatches)

##initialize species array
##someData <- rep(NA, timesteps*npatches*nspecies) 
#ar_species <- array(NA, c(nspecies, npatches, timesteps))
#ar_species[,,1]<- N

#better species array
ar_species <- array(NA, dim=c(timesteps, nspecies, npatches))

#set species to environment relationship  
H<-vector(mode='numeric',length=nspecies)
H<-seq(0,1,length=nspecies)

#initialize consumption matrix
cons <- matrix(NA, nrow = nspecies, ncol = npatches)

#initialize productivity matrix
p <- matrix(NA, nrow = nspecies, ncol = npatches)

#set dispersal matrix
dispersal_matrix<-matrix(1/(npatches-1),npatches,npatches)
diag(dispersal_matrix)<-0

#function matrix
#A<-matrix(rbinom(n=nspecies*nfunctions,size=1,prob=binom_p),nspecies,nfunctions)*runif(nspecies*nfunctions)

for(t in 1:timesteps){
		
#E[,t]<-0.5*1*(sin((2*pi/Period)*t+1+(1:npatches)*2*pi/npatches)+1)
E <-0.5*1*(sin((2*pi/Period)*t+1+(1:npatches)*2*pi/npatches)+1)


cons <- 0.1*(1.5-abs(sapply(H,'-',E)))


#dN/dt = N_ij(t)*[e*c_ij(t)*R_j(t) - m] + (a/(M-1))All patches(N(t)-aN(t))
Nt <- N*(1 + continuity*(e*R*cons - dispersalrate - m)) + continuity*(dispersal_matrix%*%N*rep(dispersalrate,npatches))
N <- Nt
	
#dR/dt
Rt <- continuity*I+R*(1-continuity*(l + rowSums(cons*N)))
R <- Rt


		
		
if(t >= 100000)	{
#Add in the species abundances into the array
for(i in 1:nspecies){
	for(j in 1:npatches){
		ar_species[t-100000,i,j] = N[j,i]
		if(N[j,i] < 0.1){
			ar_species[t-100000,i,j] = 0
		}
	}
}
}	

}
ar_species<-ar_species[seq(1,40000,100),,] #these numbers need to change as the timestep changes
#%*% - this is standard matrix multiplication, if not R multiplies matrices like it adds