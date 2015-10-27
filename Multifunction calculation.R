source("SIH code_PLT.r") #the SIH model code
require(vegan)

nspecies<-9 #the number of species
npatches<-30 #the number of patches
nfunctions<-2 #the number of functions
function_overlap<-0.5 #the amount of functional overlap

DispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1) #the dispersal rates 

#runs the SIH model at all dispersal rates in DispV and saves the abundances and productivity in a list
SIH_data<-sapply(DispV,SIH_function,species=nspecies,patches=npatches)

#define traits#### - I have added this to make some but you can define your traits any way you like
repeat{
  MTraits<-t(matrix(runif(nspecies*nfunctions),nspecies,nfunctions))
  MTraits[rbinom(nspecies*nfunctions, 1, function_overlap)==0]<-0 #0.265 is the mean proportion of species that contribute to each process in Hector and Bagchi
  MTraits<-decostand(MTraits,"total",2)
  if(sum(colSums(MTraits))==nspecies){break}
}
print(MTraits) #just shows the traits

Function_rates<-data.frame(Rate=NA,Function=factor(1:nfunctions),Dispersal=rep(DispV,each=nfunctions),Scale=rep(c("Local","Regional"),each=nfunctions*length(DispV))) #storage dataframe for functional rates
for(i in 1:length(DispV)){ #loops through all dispersal rates
  LFunc_rate1<-apply(SIH_data[["Abund",i]],1,function(M){MTraits%*%M}) #calculates the amount of each function based on the abundances from the SIH model
  LFunc_rate<-array(t(LFunc_rate1),dim=c(nrow(SIH_data[["Abund",1]]),nfunctions,npatches)) #arranges the functions into an array that is easier to use
  Function_rates[Function_rates$Scale=="Local" & Function_rates$Dispersal==DispV[i],"Rate"]<-rowMeans(apply(LFunc_rate,3,colMeans)) #calculates and saves the mean local rate for each function
  Function_rates[Function_rates$Scale=="Regional" & Function_rates$Dispersal==DispV[i],"Rate"]<-colMeans(apply(LFunc_rate,2,rowSums)) #calculates and saves the mean regoinal rate for each function
}

require(ggplot2) #plotting package
ggplot(Function_rates,aes(x=Dispersal,y=Rate,color=Function,group=Function))+ #defines the variables that you are plotting
  geom_line()+ #plots data as lines
  facet_grid(Scale~.,scale="free")+ #separates out local and regional into two panels
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
  