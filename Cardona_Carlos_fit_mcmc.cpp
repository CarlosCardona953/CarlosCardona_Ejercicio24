#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <cmath>

using namespace std; 

double fun_deg(int poly_degree, double x_obs, double* coeficientes);
double bayes(int poly_degree, double* coeficientes);
double MCMC_polynomial(int n_steps, int poly_degree);
int read_file(string filename, double* &dir);

  double *x=NULL;
  double *y=NULL;
  int n_x=0;
  

int main(){


  n_x = read_file("valores_x.txt", x);  
  n_x = read_file("valores_x.txt", y); 
  srand48(time(0));
  MCMC_polynomial(1000000, 3);
  
  return 0;
}

int read_file(string filename, double* &dir){
  ifstream infile; 
  string line;
  int n = 0;

  infile.open(filename);
  while(infile){
    
    n+=1;
    getline(infile,line);
    
  }
  infile.close();

  infile.open(filename);
  dir=new double[n];
  for (int i=0; i<n;i++){
    getline(infile,line);
    dir[i]=atof(line.c_str());
  }
  infile.close();

  return n;
}
       
double fun_deg(int poly_degree, double x_obs, double* coeficientes){

double y_fun =0;
for (int n =0; n<poly_degree+1;n++){
y_fun+=coeficientes[n]*pow(x_obs,n);
}
return y_fun;
}
            
double bayes(int poly_degree, double* coeficientes){
    double d=0;
    for (int i=0;i<n_x;i++){
    d += (-0.5*pow((y[i]-fun_deg(poly_degree, x[i], coeficientes)),2));
    }
    return d;      
}
double MCMC_polynomial(int n_steps, int poly_degree){


    double sigma = 0.1;      

    int N = n_steps;
    
     
    double* coeficientes=new double[poly_degree+1];
     
   
    for(int i=0; i<N; i++) {
        
        double* propuesta_coeficientes=new double[poly_degree+1];
        
               
        for (int j=0; j<poly_degree +1; j++){
                   
            propuesta_coeficientes[j]  = coeficientes[j] + (drand48()-0.5)*2*sigma;
           
        }

        double evaluar_viejo = bayes(poly_degree, coeficientes);
        double evaluar_nuevo = bayes(poly_degree,propuesta_coeficientes) ;
        double evaluado= evaluar_nuevo-evaluar_viejo;
        double r = min(1.0,exp(evaluado));
        double alpha = drand48();
        if(alpha<r){
             coeficientes=propuesta_coeficientes;
        } 
        cout << coeficientes[0];
        for (int k=1;k<poly_degree+1;k++){
        cout << " " << coeficientes[k];
        }
    cout << endl;
  }
}
    
