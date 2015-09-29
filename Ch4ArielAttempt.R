#Initial

#Values
nspecies = 9
npatches = 30
nfunctions = 7
e = 0.2
m = 0.2
I = 150
l = 10
binom_p = 0.5
timesteps = 400
Period = 2
dispersalrate = 1
sumforR = 0
sumforN = 0
#Random iterands...
x=0
y=0
#??need to ask Patrick about the euler approximation
#??what order should the abundance and resource equations go?
#thinking about things: to make the species respond to the environment in such a way that follows something of a trait structure, should probably just add another environmental factor into the consumption equation? like maybe make the environmental factor be a summation of lots of environmental factors and add that into the consumption equation (same with resources eventually) <- that would be for the response traits
#for the functional traits - need to futz with the production of function part of things



#Define a set of patches and their initial environments
#initialize a vector with evenly spaced initial environment variables
x = 0 #probably not needed
Enviro<- vector(mode='numeric',length=(npatches*timesteps))
for(i in 1:npatches){
	x = x + (1/npatches)
	Enviro[i] = x
}
for(i in (npatches + 1):(npatches*timesteps)){ #don't need this for-loop, auto initialized with 0's
	Enviro[i] = 0
}
E <- matrix(data = Enviro, nrow = npatches, ncol = timesteps) #might have to flip rows w columns later...


#Make a species abundance matrix, assign them all an abundance of 10 in all 30 of the patches
i_species<- vector(mode='numeric',length=(npatches*nspecies)) 
for(i in 1:(npatches*nspecies)) { 
i_species[i]=10
}
N <- matrix(data = i_species, nrow = nspecies, ncol = npatches)

#3D array for species abundances
#put t = 1 in
someData <- rep(NaN, timesteps*npatches*nspecies); 
ar_species <- array(someData, c(nspecies, npatches, timesteps));
for(i in 1:nspecies){
	for(j in 1:npatches){
		ar_species[i,j,1] = N[i,j]
	}
}


#Make a resource abundance matrix, assign them all an abundance of 9 in the t = 1 column
x=0
Resource<- vector(mode='numeric',length=(npatches*timesteps))
for(i in 1:npatches){
	Resource[i] = 9
}
for(i in (npatches + 1):(npatches*timesteps)){ #don't need this for-loop, auto initialized with 0's
	Resource[i] = 0
}
R <- matrix(data = Resource, nrow = npatches, ncol = timesteps) #might have to flip rows w columns later...

#Set H_i
#may need to change this later, not sure what sort of values he set these to
y=0 #probably not needed
H <- vector(mode='numeric',length=nspecies)
for(i in 1:nspecies){
	y = y + (1/nspecies)
	H[i] = y
}

#Set c_ij(1)
consumption<- vector(mode='numeric',length=(npatches*nspecies)) 
for(i in 1:(npatches*nspecies)) { 
consumption[i]=0
}
c <- matrix(data = consumption, nrow = nspecies, ncol = npatches)
for(i in 1:nspecies) {
	for(j in 1:npatches) {
		c[i,j]=((1.5)*abs(H[i]-E[j,1]))/10
	}
}

#Initialize productivity matrix for t = 1
e = 0.2
productivity<- vector(mode='numeric',length=(npatches*nspecies)) 
for(i in 1:(npatches*nspecies)) {
productivity[i]=0
}
p <- matrix(data = productivity, nrow = nspecies, ncol = npatches)
for(i in 1:nspecies) {
	for(j in 1:npatches) {
		p[i,j]= e*c[i,j]*R[j,1]*N[i,j] #check the indices
	}
}

#Set A_id 
binom_p = 0.5 #sets the functional overlap, repeat from above
func<- vector(mode='numeric',length=(nspecies*nfunctions)) 
for(i in 1:(nfunctions*nspecies)) {
func[i]=0
}
A <- matrix(data = func, nrow = nspecies, ncol = nfunctions)
set.seed(2)
for(d in 1:nspecies){
	for(j in 1:nfunctions){
		x = rbinom(n=1,size=1,prob=binom_p) #check the indices, not sure what the size should be...
		A[d,j]= x
	}
}

for(d in 1:nspecies){
	for(j in 1:nfunctions){
		x = runif(n=1,min = 0, max = 1) #check the indices, not sure what the size should be...
		A[d,j]= A[d,j]*x
	}
}

for(t in 2:timesteps){
	
#Environment 
	#E_j(t)=(0.5)(sin(E_initj + 2pit/T) - 1)
	for(j in 1:npatches){
		E[j,t]=(0.5)*(sin(E[j,1] + (2*pi*t)/Period) + 1) #changed the '-1' to '+1' -> + values 
	}
#Consumption Rate
	#c_ij(t) = (1.5)(abs(H_i - Ej(t)))/(10)
	#abs() is the built-in absolute value function

for(i in 1:(npatches*nspecies)) { #just to be safe
consumption[i]=0
}
c <- matrix(data = consumption, nrow = nspecies, ncol = npatches)
for(i in 1:nspecies) {
	for(j in 1:npatches) {
		c[i,j]=((1.5)*abs(H[i]-E[j,t]))/10
	}
}


#Productivity
	#productivity_ij(t) = ec_ij(t)R_j(t)N_ij(t)
e = 0.2 #for safety
for(i in 1:(npatches*nspecies)) {#for safety
productivity[i]=0
}
p <- matrix(data = productivity, nrow = nspecies, ncol = npatches)
for(i in 1:nspecies) {
	for(j in 1:npatches) {
		p[i,j]= e*c[i,j]*R[j,t]*N[i,j]
	}
}

#And the bigger functions...
#??what order should the abundance and resource equations go?
#dR/dt=I-l*R_j(t)sumofallspecies(c_ij(t)*N_ij(t))
I = 150
l = 10
#R <- matrix(data = Resource, nrow = 30, ncol = 10)
#c <- matrix(data = consumption, nrow = 9, ncol = 30)
#N <- matrix(data = i_species, nrow = 9, ncol = 30)
for(j in 1:npatches){
	for(i in 1:nspecies){
	sumforR = sumforR + c[i,j]*N[i,j]
	}
	R[j,t] = I - l*R[j,t-1] - R[j,t-1]*sumforR 
}


#dN/dt = N_ij(t)*[e*c_ij(t)*R_j(t) - m] + (a/(M-1))All patches(N(t)-aN(t))
e = 0.2
m = 0.2
dispersalrate = 1
someData <- rep(NaN, nspecies*npatches); 
N_new <- matrix(data = someData, nrow = nspecies, ncol = npatches) 
for(i in 1:nspecies){
	for(j in 1:npatches){
		for(k in 1:npatches){
			if(i != k){
			sumforN = sumforN + (N[i,k] - dispersalrate*N[i,k])	
			}
		}
		N_new[i,j] = N[i,j]*(e*c[i,j]*R[j,t] - m) + (dispersalrate/(30 - 1))*sumforN
	}
}


#Add in the species abundances into the array
for(i in 1:nspecies){
	for(j in 1:npatches){
		ar_species[i,j,t] = N_new[i,j]
		if(N_new[i,j] < 0.1){
		ar_species[i,j,t] = 0
		}
	}
}


}


