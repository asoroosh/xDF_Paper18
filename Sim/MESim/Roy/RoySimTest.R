#load necessary libraries
library(RcppArmadillo)
library(Rcpp)
library(R.matlab)
#library(MTS)
#library(Matrix)
#library(MASS)
#library(matrixcalc)
#library(statmod)
#library(R.matlab)

setwd('/home/wmrnaq/ACAnal/Sim/MESim/Roy/')

.libPaths( c( .libPaths(), "/home/wmrnaq/R/x86_64-pc-linux-gnu-library/3.3") )

itr = as.numeric(Sys.getenv("SGE_TASK_ID"))

itr     = 1; 

r_cnt   = 2;
t_cnt   = 1;

TList 	= c(100,200,600,1200,1500,2000);
RhoList = c(0,0.2,0.5,0.7,0.9); 

T       = TList[t_cnt];
rho 	= RhoList[r_cnt]; 

print(paste('For T: ',as.character(T),'  RHO: ',as.character(rho),'  ON ITR:', as.character(itr),sep=''))
 
ACn     = 6
ACnN    = ACn*(ACn-1)/2+ACn #Total number of available AC structures

print('Libs are on!')

#bring in Rcpp functions
sourceCpp('bin/LongitInclude.cpp')

print('RoysEstimator is ON, ready for take off!')

Path2Storage = '/storage/essicd/data/HCP/Soroosh/ACAnal/MESim/'

RoyVar = matrix(0, ACn, ACn)
Idx = which(upper.tri(matrix(1, ACn, ACn), diag = 1))

for (AC_cnt in 1:ACnN){
  print(paste('ACStructure: ',as.character(AC_cnt),sep=''))
  
  WhereFrom = paste(Path2Storage,'R_tsSim_t', as.character(t_cnt) ,'_r', as.character(r_cnt) ,'/ts_AC', as.character(AC_cnt),'_',as.character(T),'_r',as.character(rho),'_', as.character(itr) ,'.txt' ,sep = '')
  print(WhereFrom)

	#Read the time series
  ts = t(as.matrix(read.csv(WhereFrom,sep = ',',header = FALSE)))
  
	#Do the Roy variance estimator on them 
  RoyFunc = Roy(ts,lag = 20,bw = 5,SigmaType = 'Diagonal') #BW was copied from Mark's paper, number of lags copied from his code -- couldn't find anything about it in the paper --
  RoyVar[Idx[AC_cnt]]  = RoyFunc[[3]]
}

#Save them as a matlab variable 
WhereOnStorage = paste(Path2Storage,'R_MESim_t', as.character(t_cnt) ,'_r', as.character(r_cnt) ,'/',sep='')

print(WhereOnStorage)

Where2Save = paste(WhereOnStorage,'MESim_t',as.character(T),'_r',as.character(rho),'_',as.character(itr),'_Roy_TVOff.mat' ,sep = '');
writeMat(RoyVar = RoyVar, Where2Save)

