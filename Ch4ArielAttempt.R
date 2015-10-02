#Initial

#run one in R and the other in R studio

#Values
nspecies = 9
npatches = 30
nfunctions = 7
euler = 0.08
e = 0.2
m = 0.2
I = 150
l = 10
binom_p = 0.5
timesteps = 140
Period = 40
dispersalrate = 0.01 
sumforR = 0
sumforN = 0
#Random iterands...
x=0
y=0

#Tip: Use NA instead of NaN, NaN can cause some problems, seq() is super useful

#initialize abundance matrix
N<-matrix(10,nrow = nspecies,ncol = npatches)

#initialize resource matrix
R<-matrix(NA,nrow = npatches, ncol = timesteps)
R[,1]<-9

#initialize environment matrix 
E <- matrix(NA, nrow = npatches, ncol = timesteps)
E[,1] <-seq(0,1,length=npatches)

#initialize species array
someData <- rep(NA, timesteps*npatches*nspecies) 
ar_species <- array(someData, c(nspecies, npatches, timesteps))
ar_species[,,1]<- N

#initialize species to environment relationship  
H<-vector(mode='numeric',length=nspecies)
H<-seq(0,1,length=nspecies)

#initialize consumption matrix
c <- matrix(NA, nrow = nspecies, ncol = npatches)

#initialize productivity matrix
p <- matrix(NA, nrow = nspecies, ncol = npatches)

#function matrix
A<-matrix(rbinom(n=nspecies*nfunctions,size=1,prob=binom_p),nspecies,nfunctions)*runif(nspecies*nfunctions)

for(t in 1:timesteps){
	

for(j in 1:npatches){
	E[j,t]=(0.5)*(sin(E[j,1] + (2*pi*t)/Period) + 1) #changed the '-1' to '+1' -> + values 
	}


for(i in 1:nspecies) {
	for(j in 1:npatches) {
		c[i,j]=((1.5) - abs(H[i] - E[j,t]))/10
	}}


#productivity_ij(t) = ec_ij(t)R_j(t)N_ij(t)
e = 0.2 #safety
for(i in 1:nspecies) {
	for(j in 1:npatches) {
		p[i,j]= e*c[i,j]*R[j,t]*N[i,j]
	}
}

#dR/dt=I-l*R_j(t)sumofallspecies(c_ij(t)*N_ij(t))
I = 150
l = 10
for(j in 1:npatches){
	for(i in 1:nspecies){
	sumforR = sumforR + c[i,j]*10
	}
	if(t == 1){
	R[j,t+1] = R[j,t] + euler*(I  - l*R[j,t] - R[j,t]*sumforR)
	}else{
		R[j,t] = R[j,t-1] + euler*(I  - l*R[j,t-1] - R[j,t-1]*sumforR)
		}
	}

#dN/dt = N_ij(t)*[e*c_ij(t)*R_j(t) - m] + (a/(M-1))All patches(N(t)-aN(t))
e = 0.2 #safety
m = 0.2 #safety
dispersalrate = 0.01 #safety
#N_new <- matrix(NA, nrow = nspecies, ncol = npatches) 
for(i in 1:nspecies){
	for(j in 1:npatches){
		for(k in 1:npatches){
			if(k != j){
			sumforN = sumforN + (N[i,k]) 	
			}
		}
		
	N[i,j] = N[i,j] + euler*N[i,j]*(e*c[i,j]*R[j,t] - dispersalrate - m) + euler*(dispersalrate/(npatches - 1))*sumforN 
	#N[i,j] = N_new[i,j]
}
}	
	
		
		
if(t > 1)	{
#Add in the species abundances into the array
for(i in 1:nspecies){
	for(j in 1:npatches){
		ar_species[i,j,t] = N[i,j]
		#if(N[i,j] < 0.1){
		#ar_species[i,j,t] = 0
		#}
	}
}
}	

}
#%*% - this is standard matrix multiplication, if not R multiplies matrices like it adds