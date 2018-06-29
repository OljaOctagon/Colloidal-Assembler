#include "sphere.h"


	sphere::sphere(){

		 
		 trans_periodic = new double[3];
		
  
        }  


    sphere::~ sphere(){
  
         delete[] trans_periodic;
		
        }  
