source("SIH code_PLT.r")
require(vegan)

nspecies<-9
npatches<-30
nfunctions<-2
function_overlap<-0.5

DispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)

#runs the SIH model at all dispersal rates in DispV and saves the abundances and productivity in a list
SIH_data<-sapply(DispV,SIH_function,species=nspecies,patches=npatches)

#define traits#### - I have added this to make some but you can define your traits any way you like
repeat{
  MTraits<-t(matrix(runif(nspecies*nfunctions),nspecies,nfunctions))
  MTraits[rbinom(nspecies*nfunctions, 1, function_overlap)==0]<-0 #0.265 is the mean proportion of species that contribute to each process in Hector and Bagchi
  MTraits<-decostand(MTraits,"total",2)
  if(sum(colSums(MTraits))==nspecies){break}
}
print(MTraits)

Function_rates<-data.frame(Rate=NA,Function=factor(1:nfunctions),Dispersal=rep(DispV,each=nfunctions),Scale=rep(c("Local","Regional"),each=nfunctions*length(DispV)))
for(i in 1:length(DispV)){
  LFunc_rate1<-apply(SIH_data[["Abund",i]],1,function(M){MTraits%*%M})
  LFunc_rate<-array(t(LFunc_rate1),dim=c(nrow(SIH_data[["Abund",1]]),nfunctions,npatches))
  Function_rates[Function_rates$Scale=="Local" & Function_rates$Dispersal==DispV[i],"Rate"]<-rowMeans(apply(LFunc_rate,3,colMeans))
  Function_rates[Function_rates$Scale=="Regional" & Function_rates$Dispersal==DispV[i],"Rate"]<-colMeans(apply(LFunc_rate,2,rowSums))
}

require(ggplot2)
ggplot(Function_rates,aes(x=Dispersal,y=Rate,color=Function,group=Function))+ #defines the variables that you are plotting
  geom_line()+ #plots data as lines
  facet_grid(Scale~.,scale="free")+ #separates out local and regional into two panels
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
  