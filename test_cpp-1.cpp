/*
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.

*/
#include <RcppArmadillo.h>
#include <cmath>
#include<algorithm> 
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <math.h>

using namespace Rcpp ;
using namespace arma ;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double UFVhd_cpp(const arma::colvec& v, const arma::colvec& bet,const arma::mat& X)
  {
  int n = X.n_rows, p = X.n_cols;
  arma::mat w(n,p), v1(n,1);
  arma::colvec Nvect, Dvect, quotient;
  w.cols(1,p-1)=X.cols(0,p-2);
  v1.fill(1.0);
  w.col(0)=v1;
  Nvect=X.col(p-1)-w*bet;
  Dvect=w*v;
  for (int i=0; i<n; ++i){
    if(Dvect(i)==0){Dvect(i)=pow(10, -20);} 
    }  ;
  quotient=Nvect/Dvect;
  //Rprintf("median value: %f \n", median(quotient));
  double ufv=std::abs(median(quotient));
  return ufv;          
  }

// [[Rcpp::export()]]
double compute_UF_cpp(const arma::mat X, const arma::colvec bet,
                      int ND, double UF_min)
{
  //std::cout<<"test50:\n"<<std::endl;
  int n = X.n_rows, p = X.n_cols;
  double ufold=0.0;
  
  //std::cout<<"test5:\n"<<std::endl;
  
//generate T_i's    
  arma::mat w(n,p), TM(n,p), v1(n,1), DI(n,n);
  arma::colvec Res;
  
  //std::cout<<"test6:\n"<<std::endl;
  
  w.cols(1,p-1)=X.cols(0,p-2);
  v1.fill(1.0);
  w.col(0)=v1;
  Res=X.col(p-1)-w*bet;
  
  //std::cout<<"test7:\n"<<std::endl;
  
  for (int i=0; i<n; ++i){
    if(Res(i)==0){Res(i)=pow(10, -20);} 
  };
  Res=v1/Res;
  
  DI=diagmat(Res);
  TM=DI*w;
//  TM.rows(0,9).print();
//----------------------------------------------------------------
//take care of directions that are perpendicular to T_i's
for (int i=0; i<n; ++i)
 {
  arma::colvec uu(p); 
  arma::uvec id;
  uu.fill(0.0);
 // if(i==0) {TM.row(i).print();}
  id=find(TM.row(i));// TM_i is a 1xp rowvec, id is its non-zero indeices
  //if (i==0){id.t().print();} std::cout<<"id.n_rows"<<id.n_rows<<std::endl;}
  if(id.n_rows==1)    // contruct a vector uu that is perpendicular to TM_i
  { if (id(0)==0)
    { 
    uu(0)= TM.row(i)(1)  ;
    uu(1)= -TM.row(i)(0) ;
    }
    else
    { uu(0)=-TM.row(i)(id(1));
      uu(id(1))=TM.row(i)(0);
    }  
  }  
  else  
  {
   uu(id(0))=-TM.row(i)(id(1));
   uu(id(1))=TM.row(i)(id(0));  
  }
  //if(i==0) {uu.t().print(); TM.row(i).print();}
  double res=arma::sum(uu.t()*uu);
  arma::vec denom(p);
  arma::vec v=uu/denom.fill(std::sqrt(res));
  //if(i==0) {v.t().print();}
  double ufnew=UFVhd_cpp(v, bet, X);
  
  //std::cout<<"test8:\n"<<std::endl;
  if (ufnew>=UF_min){ ufold=pow(10,10); return(ufold);}
    else {ufold=std::max(ufold, ufnew);}
        
 }//for end brace for the TM_i

  //std::cout<<"test9:\n"<<std::endl;
 //now take care of p axis directions  
   arma::mat D(p,p); 
   D.eye();
   for (int i=0; i<p; ++i)
   { 
     arma::vec v=D.col(i);
     double ufnew=UFVhd_cpp(v,bet, X);
   
       
     if (ufnew>=UF_min){ufold=pow(10,10);return(ufold);}
     else{ufold=std::max(ufold, ufnew);}
   } 
 //std::cout<<"test10:\n"<<std::endl;
//now considere the N Normal directions of the hyperplanes formed by p points   
 
   for (int J=0; J<ND; ++J)
   {
    Rcpp::IntegerVector sequence, id(p), idx, index_2;
    arma::mat X_temp(p,p), y_temp(p,p), v1(p,1), yvect(p,1);
    sequence=Rcpp::seq_len(n);
    arma::vec vv_1=randu<vec>(n);
    index_2=sort_index(vv_1);

    //std::cout<<"test11:\n"<<std::endl;
    for (int i=0; i<p; ++i) { id(i)=index_2(i);}
//  id=RcppArmadillo::sample(sequence, p, false);
    for (int j=0; j<p; ++j){ 
    y_temp.row(j)=X.row(id(j)); 
    yvect.row(j)=X(id(j),p-1); 
    }
    X_temp.cols(1,p-1)=y_temp.cols(0,p-2);  
    v1.fill(1.0);
    X_temp.col(0)=v1;
    
    //std::cout<<"test12:\n"<<std::endl;
    arma::vec uu(p);
    arma::mat vv(p,1), beta1(p-1,1);
    
    vv=arma::inv(X_temp)*yvect; //this is (beta_0, beta_1) in y=(1,x')beta
    beta1=vv.rows(1, p-1);
    uu.rows(0,p-2)=-beta1; //here could use arma::join_rows(A,B) function
    uu(p-1)=1;         //but A and B must have the same number of rows
   
    double res=arma::sum(uu.t()*uu);
    arma::vec denom(p);
    arma::vec v=uu/denom.fill(std::sqrt(res));
    double ufnew=UFVhd_cpp(v, bet, X);
   
   if (ufnew>=UF_min){ufold=pow(10,10);return(ufold);}
   else      {ufold=std::max(ufold, ufnew);} 
   //std::cout<<"J:~"<<J<<std::endl;
 }  
   
   return(ufold); // return a arma:vec ret(2); ret(0)=ufold; ret(1)=UF_min;
}//function end brace


// [[Rcpp::export()]]
arma::vec compute_deepest_PRD2_cpp(const arma::mat& X, const arma::mat& B,
 const arma::vec& UFbeta, int ND, int Nbet, double UF_min, int RN)
{
 int nn = B.n_rows, p = B.n_cols, N2=nn-RN;
 arma::uvec index=arma::sort_index(UFbeta);
 int id_min=index[0];
 arma::uvec J=index.rows(0,p);
 
 int nrow=Nbet+N2+1; //plus one make sure that zeros(p) as a candidate of beta is evaluated
 arma::mat beta(nrow,p); //mat beta(nrow,p,fill::zeros)
 beta.zeros();
 arma::vec UFbeta_final(nrow);
 
 //first treat beta=zeros(p), the last row of beta 
 int k=nrow-1;
 double temp0=compute_UF_cpp(X, beta.row(k).t(), ND, UF_min);
 UFbeta_final(k)= temp0; 
 UF_min=std::min(UF_min, temp0); 
 
 for (int k=0; k<Nbet; ++k)
 {
  arma::vec b=randu<vec>(p+1);
  arma::vec dem= Rcpp::rep(sum(b), p+1);
  arma::vec a=b/dem; //repmat(sum(b), p+1);
  
  beta.row(k)=a.t()*B.rows(J); 
  double temp=compute_UF_cpp(X, beta.row(k).t(), ND, UF_min);
  UFbeta_final(k)= temp; 
   UF_min=std::min(UF_min, temp); 
 }
 
 arma::uvec JJ(N2);
 for (int ii=RN; ii<nn; ++ii) {JJ[ii-RN]=ii;}
// Rcpp::IntegerVector JJ=Rcpp::seq(RN, (n-1));
 for (int l=0; l<N2; ++l)
 {
   arma::vec b1=randu<vec>(N2);
   arma::vec dem1=Rcpp::rep(sum(b1), N2);
   arma::vec a=b1/dem1;
  
   int k=l+Nbet;
   beta.row(k)=a.t()*B.rows(JJ);
//   if(arma::sum(beta.row(k))==0.0)//(beta.row(k).has_nan())
//   {std::cout<<
//     "beta.row(k) is a zero vector:\n"<<(beta.row(k))
//     //<<"a.t(): "<<a.t()
//       <<std::endl;}
   UFbeta_final(k)=compute_UF_cpp(X, beta.row(k).t(), ND, UF_min);
   UF_min=std::min(UF_min, UFbeta_final(k)); 
 }
 
 double UF_min_final=min(UFbeta_final);
 uvec index_fin=arma::find(UFbeta_final==UF_min_final);
 int id_fin=index_fin(0);

// if (B.has_nan()){std::cout<<"B contains NAN and needed to be cleaned"<<std::endl;}
// if (beta.has_nan()){std::cout<<"beta contains NAN and needed to be cleaned"
//                    <<std::endl;}
/*
 std::cout<<"\n first 5 rows of beta:\n"<<std::endl;
 beta.rows(0,4).print();
 std::cout<<"\n"<<std::endl;
 std::cout<<"\n first 5 rows of B2:\n"<<std::endl;
 B.rows(0,4).print();
 std::cout<<"\n"<<std::endl;

 for (int m=0; m<nrow; ++m)
   {if(sum(beta.row(m)==0)) 
   {std::cout<<"\n nrow=:\n"<<nrow<<"Nbet=:\n"<<Nbet<<"N2=:\n"<<N2<<
     "beta mth row is zero and m=:\n"
             <<beta.row(m)<<m<<std::endl;
     }
   }
*/
  vec  beta_vec(p);
 if (UFbeta_final[id_fin]< UFbeta[id_min] )
 {beta_vec=beta.row(id_fin).t(); //  double UF=UFbeta_final[id_fin];
 }
 else
 {beta_vec=B.row(id_min).t();//  double UF=UFbeta[id_min];
 }
// std::cout<<"beta.row(id_fin).t() and B.row(id_min).t():\n" << beta.row(id_fin)
//          <<"\n"<<B.row(id_min)<<std::endl;
// std::cout<<"beta_vec=:\n"<<beta_vec.t()<<std::endl;
return  beta_vec; //Rcpp::list::crate(Rcpp::Named("beta")=beta,
              //           Rcpp::Named("UF")=UF);  
}  //end of the big function


// [[Rcpp::export]]
List get_Beta_UFbeta_cpp(const mat X, const mat Z,const mat dpbeta, 
                     int ND, int N2, int RN)
{
  int p=X.n_cols, nrow=RN+N2;  
  double  UF_min=pow(10,10);
  arma::vec UFbeta2(nrow); 
  arma::rowvec beta;
  
  for (int j=0; j<RN; ++j)
  { 
    beta=Z.row(j); 
    UFbeta2(j)= compute_UF_cpp(X,beta.t(),ND,UF_min);
    if (j>=(p+1)){UF_min=std::min(UF_min, UFbeta2(j));}
  }  

  for (int kk=0; kk<N2; ++kk)
  {
    int ll=kk+RN; 
    UFbeta2(ll)=compute_UF_cpp(X,dpbeta.row(kk).t(), ND, UF_min);
    UF_min=std::min(UF_min, UFbeta2(ll));
  } 
  
  List res;
  res["UFbeta2"]=UFbeta2;
  res["UF_min"]=UF_min;
  
  return(res);
} //end of the big

//added on 2/24/20

// [[Rcpp::export]]
List get_dprd_wprd(//const arma::mat& X,
                   const arma::mat& B,
                   const arma::vec& UFbeta, //int ND, 
                   double K)
{
 int //n=B.n_rows, 
   p=B.n_cols;
// int m=floor(alpha*n);  
 arma::uvec index=arma::sort_index(UFbeta);
 arma::vec dprd1;
 
 //std::cout<<"T1:"<<std::endl;
 arma::rowvec dpofB=B.row(index[0]);
 arma::uvec J=index.rows(0,p);
 
 //std::cout<<"T2:"<<std::endl;
 arma::mat temp=B.rows(J);
 dprd1=sum(temp).t()/(p+1);
 
 //std::cout<<"T3:"<<std::endl;
 arma::vec weight=zeros(p+1);
 double r0=UFbeta(index(p-1));
 
 for (int i=0;i<p+1;++i)
 { double r=UFbeta(index(i));
   weight(i)=(r<=r0)+
     (r>r0)*(exp(K*(2*r0/r-pow(r0/r,2)))-1)/(exp(K)-1);
 }
 //std::cout<<"T4:"<<std::endl;
 arma::rowvec temp2=weight.t()*B.rows(J);
 arma::rowvec dprd2=temp2/sum(weight);

/* the following attempt is not very successful  
// std::cout<<"T4:"<<std::endl;
 //arma::mat deepmat=randu<mat>(p,3);
// std::cout<<"T44:"<<std::endl;
 List deepmat;
 deepmat["0"]=dpofB;
 deepmat["1"]=dprd1;
 deepmat["2"]=dprd2;
// deepmat=join_cols(dpofB, dprd1, dprd2);
// std::cout<<"T5:"<<std::endl;
 
 double UF_min=pow(10,10);
 arma::vec zer0=zeros(p);
 UF_min=compute_UF_cpp(X, zer0, ND, UF_min);
 
// std::cout<<"T6:"<<std::endl;
 arma::vec UFdeep(3);
// std::cout<<"T66:"<<std::endl;
 UFdeep(0)= UFbeta(index(0));    //compute_UF_cpp(X,dpofB.t(), ND, UF_min);
// std::cout<<"T666:"<<std::endl;
 double t2, t3;
 arma::vec b1=dprd1;
 arma::vec b2=dprd2.t();
// std::cout<<"T6666:"<<std::endl;
 t2=compute_UF_cpp(X, b1, ND, UF_min);
 t3=compute_UF_cpp(X, b2, ND, UF_min);
// std::cout<<"T66666:"<<std::endl;
 UFdeep(1)= t2;//compute_UF_cpp(X,dprd1.t(), ND, UF_min);
 UFdeep(2)= t3;//compute_UF_cpp(X,dprd2.t(), ND, UF_min);
// std::cout<<"T6666:"<<std::endl;
 arma::uvec index1=arma::sort_index(UFdeep);
 int k=index1(0);
 arma::rowvec dprd22=deepmat[k];
// std::cout<<"T7:"<<std::endl;
*/
 List res;
 res["dprd"]=dpofB; //deepest of among B
 res["wprd1"]=dprd1; //weighted based on deepest (p+1) beta's (mean)
 res["wprd2"]=dprd2; //weighted based on all beta's in B (UF weighted)
 // res["wprd3"]=dprd3; //wieughted based on all beta's in Beta generated 
 //                     //in compute_deepest_PRD2
 
 return(res);
}  
