#include <iostream>
#include <string>
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;
using namespace std;


//////////////////////////////    Functions used in cross sectional and Longitudinal Models //////////

//calculate cross covariance for an individual visit
// [[Rcpp::export]]
cube ccv(const mat& xt, int lag){ 
  //xt is matrix for a given visit moving across ROIs and down time
  //lag is maximum number of lags to use in the calculations
  
  int P = xt.n_cols; //number of ROIs
  int T = xt.n_rows; //length of time series
  cube covMat(P, P, 2*lag+1, fill::zeros);
  vec xi(T, fill::zeros);
  double xibar;
  vec xj(T, fill::zeros);
  double xjbar;
  
  for(int i = 0; i < P; i++){
    for(int j = 0; j < P; j++){
      xi = xt.col(i);
      xj = xt.col(j);
      xibar=mean(xi);
      xjbar=mean(xj);
      for(int l = 0; l <= lag; l++){
        vec xil = xi.head(T-l);
        vec xjl = xj.tail(T-l);
        double sum = 0;
        for(int t=0; t<(T-l); t++){
          sum += (xil(t)-xibar)*(xjl(t)-xjbar);
        }
        covMat(i,j,lag+l) = covMat(j,i,lag-l) = sum/T;//cross covariance for ROI i and j at lag l
      }
    }
  }
  return covMat; //return cube with ROIs in rows and columns and lags in slices (lag 0 is middle slice)
}


//function to calculate the lag 0 cross correlation matrix for a given visit
// [[Rcpp::export]]
mat ccr(const cube& C){
  //C is a PxPxL cross covariance cube output by the ccv function
  
  int P = C.n_cols; //Number of ROIs
  int L = C.n_slices; //2 times the number of lags plus 1
  mat R(P,P,fill::zeros);
  mat lag0cov = C.slice((L-1)/2); //take the lag 0 slice from covariance cube
  for(int i = 0; i < P; i++){
    for(int j = 0; j < P; j++){
      R(i,j)=lag0cov(i,j)/std::sqrt(lag0cov(i,i)*lag0cov(j,j)); //convert to correlation
    }
  }
  return R; //return correlation matrix for lag 0
}

 
//calculate weights to be used in the estimation of delta
// [[Rcpp::export]]
arma::vec MB(double lag, double b){
  //lag is maximum number of lags to use in calculation
  //b is the bandwidth of the windowing function times sqrt(T) where T is the length of the time series
  
  vec weights(2*lag+1, fill::zeros);
  for(int i = -lag; i < lag+1; i++){
    if(std::abs(i/b)<=1){weights[i+lag] = 1-std::abs(i/b);}
    else{weights[i+lag] = 0;}
  }
  return weights; //return vector of weights of length 2 times lag + 1 (lag 0 in center)
}

  
//function to calculate delta which are used in calculation of Roy covariance cov(r_ij, r_lm)
// [[Rcpp::export]]
double delta(int i, int j, int l, int m, const vec& weights, const cube& C, int lag, int lagln){
    //i, j, l, and m are the ROIs for the desired delta value
    //weights is a vector of output weights from the MB function
    //C is cross covariance cube from ccv function
    //lag is maximum lag to be used in calculations
    //lagln is 2*lag+1
    
    double sum=0;
    for(int t=0; t<lagln; t++){
      sum += weights(t)*C(i-1,j-1,t)*C(l-1,m-1,t); //estimate theta
    }
    double result = sum/std::sqrt(C(i-1,i-1,lag)*C(j-1,j-1,lag)*C(l-1,l-1,lag)*C(m-1,m-1,lag)); //convert to correlation
    return result;
}


//function calculate Roy covariance cov(r_ij, r_lm)
// [[Rcpp::export]]
double roycov(int i, int j, int l, int m, const cube& C, const mat& R, double T, 
              const vec& weights, int lag, int lagln){
    //i,j,l,m are ROIs of the desired covariance
    //C is cross covariance cube from ccv function
    //R is lag 0 cross correlation from ccr function
    //T is the length of the time series
    //weights is a vector of weights from the MB function(or other windowing function)
    //lag is maximum lag to be used in calculations
    //lagln is 2*lag+1
    
    //formula for Roy covariance
    double result=(0.5*R(i-1,j-1)*R(l-1,m-1)*(delta(i,l,i,l,weights,C,lag,lagln)+
                                              delta(i,m,i,m,weights,C,lag,lagln)+
                                              delta(j,l,j,l,weights,C,lag,lagln)+
                                              delta(j,m,j,m,weights,C,lag,lagln)) -
      R(i-1,j-1)*(delta(i,l,i,m,weights,C,lag,lagln)+delta(j,l,j,m,weights,C,lag,lagln)) -
      R(l-1,m-1)*(delta(j,l,i,l,weights,C,lag,lagln)+delta(j,m,i,m,weights,C,lag,lagln)) +
      delta(i,l,j,m,weights,C,lag,lagln) + delta(j,l,i,m,weights,C,lag,lagln))/T;
    return result;
}


//construct within subject covariance matrix using Roy covariance function, unstructured
// [[Rcpp::export]]
mat RoyMat(const cube& C, const mat& R, double T, const vec& weights, int lag, string SigmaType){
  //C is cross covariance cube from ccv function
  //R is lag 0 cross correlation matrix from ccr function
  //T is the length of the time  series
  //weights is a vector of weights from the MB function(or other windowing function)
  //lag is maximum lag to be used in calculations
  //SigmaType is a string to signify if SigmaRoy variance should be "Zero", "Unstructured", or "Diagonal"
  
  int P=C.n_rows; //number of ROIs
  int Q=P*(P-1)/2; //number of ROI pairs
  int lagln=2*lag+1; 
  
  mat result(Q,Q,fill::zeros);
  
  if(SigmaType=="Unstructured"){
    vec resultvec(Q*Q, fill::zeros);
  
    //create vector of roy covariances
    int count=0;
    for(int i=0; i<(P-1); i++){
      for(int j=i+1; j<P; j++){
        for(int l=0; l<(P-1); l++){
          for(int m=l+1; m<P; m++){
            if((l > i) || (l == i && m >= j)){
            resultvec(count) = roycov(i+1,j+1,l+1,m+1,C,R,T,weights,lag,lagln);
            count += 1;
            }
          }
        }
      }
    }
  
    //populate within subject Roy variance matrix from the vector of Roy variances
    for(int i=0; i<Q; i++){
      for(int j=i; j<Q; j++){
        result(i,j)=resultvec[0];
        result(j,i)=resultvec[0];
        resultvec=resultvec.tail(resultvec.n_elem-1);
      }
    }
  }
  
  if(SigmaType=="Diagonal"){
    int count=0;
    for(int i=0; i<(P-1); i++){
      for(int j=i+1; j<P; j++){
        result(count,count) = roycov(i+1,j+1,i+1,j+1,C,R,T,weights,lag,lagln);
        count += 1;
      }
    }
  }
  
  if(SigmaType=="Zero"){
    result.fill(0);
  }
  
  return result;
}


//Apply previous functions to get within subject Roy variance for a single visit
// [[Rcpp::export]]
List Roy(const mat& x, int lag, double bw, string SigmaType){ 
  //x is matrix for a given subject moving across ROIs and down time
  //lag is maximum number of lags to use in the calculations
  //SigmaType is a string to signify if SigmaRoy variance should be "Zero", "Unstructured", or "Diagonal"
  //bw is the bandwidth to be used by the MB function
  
  int T=x.n_rows; //length of time series
  lag=min(lag, T);
  cube C=ccv(x,lag); //cross covariance cube
  mat R=ccr(C); //lag 0 correlation matrix
  
  double b = bw*sqrt(T); //bandwidth used in windowing function
  vec weights=pow(MB(lag,b),2); 
  
  List out(2);
  out("Sigma")=RoyMat(C,R,T,weights,lag,SigmaType);
  
  out("R")=R;
  return out;
}