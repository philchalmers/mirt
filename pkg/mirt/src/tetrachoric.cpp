#include "TetCorr.h"
#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>

extern "C"{
SEXP tetrachoric(SEXP Rdatavec) {
		
	int *datavec, npairs, Table[2][2], k = 0;    
	
	//Make pointers and protect variables	
	PROTECT(Rdatavec = AS_INTEGER(Rdatavec));	
	datavec = INTEGER_POINTER(Rdatavec);	
	npairs = LENGTH(Rdatavec) / 4;	
	
	SEXP Rreturn;	
	double *Rs;	
	PROTECT(Rreturn = NEW_NUMERIC(npairs));
	Rs = NUMERIC_POINTER(Rreturn);		
	
	double R = 0.0; // Resulting tetrachoric r.
	double SE = 0.0; //Standard error of r.
	double SE0 = 0.0; //Standard error of r for test that rho = 0.
	double Z1 = 0.0;
	double Z2 = 0.0;
	double PHI_R = 0.0;
	int ITYPE = 0;	

	/*        IFAULT:  Termination code:
	*                  0:  normal termination
	*                  1:  iteration failed to converge
	*                  2:  table has at most one nonzero row/column or at
	*                     least one negative frequency.
	*/
	int IFAULT = 0, IFAULTVEC[npairs];	
	TetCorr tetCorr;
	
	for(int i = 0; i < npairs; i++) {	  
	//int i = 0;
	  Table[0][0] = datavec[0 + k];
	  Table[1][0] = datavec[1 + k];
	  Table[0][1] = datavec[2 + k];
	  Table[1][1] = datavec[3 + k];	  
	  tetCorr.Tetra(Table,R,SE,SE0,Z1,Z2,PHI_R,ITYPE,IFAULT);
	  Rs[i] = R;	  
	  IFAULTVEC[i] = IFAULT;
	  k += 4;
	}    

	UNPROTECT(2);	
	return(Rreturn);	

}	
}
