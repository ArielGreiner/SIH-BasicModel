#source("SIH code_PLT.r") #the SIH model code
require(vegan)
#SIH_data<-sapply(DispV,SIH_function,species=nspecies,patches=npatches,eff_vary=T,eff_values=eff_values)

nspecies<-9 #the number of species
npatches<-30 #the number of patches
nfunctions<-1 #the number of functions #CURRENTLY ONLY LOOKING AT PRODUCTION OF 1 FUNCTION
function_overlap<-0.5 #the amount of functional overlap

DispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1) #the dispersal rates 

#runs the SIH model at all dispersal rates in DispV and saves the abundances and productivity in a list
eff_values<-rnorm(nspecies,mean=0.2,sd=0.005)
SIH_data<-sapply(DispV,SIH_function,species=nspecies,patches=npatches)

#define traits#### - I have added this to make some but you can define your traits any way you like
MTraits<-t(matrix(runif(nspecies*nfunctions)*rbinom(nspecies*nfunctions,1,function_overlap),nspecies,nfunctions)) #no structure 
#MTraits[1,]<-c(1,1,1,0,0,0,0,0,0)
#MTraits[2,]<-c(0,0,0,1,1,1,0,0,0)
#MTraits[3,]<-c(0,0,0,0,0,0,1,1,1)
#print(MTraits) #just shows the traits


Function_rates<-data.frame(Rate=NA,Function=factor(1:nfunctions),Dispersal=rep(DispV,each=nfunctions),Scale=rep(c("Local","Regional"),each=nfunctions*length(DispV))) #storage dataframe for functional rates
Function_number<-data.frame(Number=NA,Dispersal=DispV,Scale=rep(c("Local","Regional"),each=length(DispV))) #storage dataframe for functional rates

for(i in 1:length(DispV)){ #loops through all dispersal rates
  LFunc_rate1<-apply(SIH_data[["Abund",i]],1,function(M){MTraits%*%M}) #calculates the amount of each function in each patch based on the abundances from the SIH model
  LFunc_rate<-array(t(LFunc_rate1),dim=c(nrow(SIH_data[["Abund",1]]),nfunctions,npatches)) #arranges the functions into an array that is easier to use
  Function_number[Function_number$Scale=="Local" & Function_number$Dispersal==DispV[i],"Number"]<-mean(apply(LFunc_rate>0,3,rowSums))
  Function_number[Function_number$Scale=="Regional" & Function_number$Dispersal==DispV[i],"Number"]<-mean(rowSums(apply(LFunc_rate,2,rowSums)>0))
  #Single function variant:
  Function_rates[Function_rates$Scale=="Local" & Function_rates$Dispersal==DispV[i],"Rate"]<-mean(apply(LFunc_rate,3,colMeans))
  #multifunction variant:
  #Function_rates[Function_rates$Scale=="Local" & Function_rates$Dispersal==DispV[i],"Rate"]<-rowMeans(apply(LFunc_rate,3,colMeans)) #calculates and saves the mean local rate for each function
  Function_rates[Function_rates$Scale=="Regional" & Function_rates$Dispersal==DispV[i],"Rate"]<-colMeans(apply(LFunc_rate,2,rowSums)) #calculates and saves the mean regoinal rate for each function
}

require(ggplot2) #plotting package
ggplot(Function_rates,aes(x=Dispersal,y=Rate,color=Function,group=Function,linetype=Function))+ #defines the variables that you are plotting
  geom_line(size=2)+ #plots data as lines
  facet_grid(Scale~.,scale="free")+ #separates out local and regional into two panels
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(Function_number,aes(x=Dispersal,y=Number,group=Scale,color=Scale,linetype=Scale))+ #defines the variables that you are plotting
  geom_line(size=2)+ #plots data as lines
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
  