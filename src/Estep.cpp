#include <Rcpp.h>

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

//Estep for mirt
RcppExport SEXP Estep(SEXP Ritemtrace, SEXP Rprior, SEXP RX,  
	SEXP Rnfact, SEXP Rr) 
{
    BEGIN_RCPP

    /*
        Ritemtrace = numeric matrix. Probability traces
        Rprior = numeric vector. Normalized prior
        RX = integer matrix. Response vector with 9's as missing
        Rnfact = integer. Number of factors
        Rr = integer vector.  Counts of same response vector
    */
    
    Rcpp::NumericVector prior(Rprior);
    Rcpp::IntegerVector nfact(Rnfact);
    Rcpp::IntegerVector r(Rr);
    Rcpp::IntegerMatrix data(RX);
    Rcpp::NumericMatrix itemtrace(Ritemtrace);
    int nquad = prior.length();
    int nitems = data.ncol();
    int npat = r.length();      
    double expd=0.0, posterior[nquad];
    int i=0, k=0, item=0;	
    Rcpp::NumericMatrix r1(nitems, nquad);    
    Rcpp::NumericMatrix r0(nitems, nquad);
    Rcpp::NumericVector expected(npat);
	
    // Begin main function body 				
	for (int pat = 0; pat < npat; pat++){		
  	    for (k = 0; k < nquad; k++)
		    posterior[k] = prior[k];
			
		for (item = 0; item < nitems; item++){
		    if (data(pat,item) == 9) 
			    continue;
		    if (data(pat,item)){
                for (k = 0; k < nquad; k++)
			        posterior[k] = posterior[k] * itemtrace(item,k);
    		} else {
	    	    for (k = 0; k < nquad; k++)
		    	    posterior[k] = posterior[k] * (1.0 - itemtrace(item,k));
		    }			
		}    
	    expd = 0.0;
	    for (i = 0; i < nquad; i++)
	        expd += posterior[i];	
	    expected[pat] = expd;		
		
	    for (i = 0; i < nquad; i++)
	        posterior[i] = r[pat] * posterior[i] / expd;	
			
        for (item = 0; item < nitems; item++){
	        if (data(pat,item) == 9) 
	    	    continue;
    	    if (data(pat,item)){
	            for (k = 0; k < nquad; k++)
	                r1(item,k) += posterior[k];
	    	} else {
		        for (k = 0; k < nquad; k++)
    		        r0(item,k) += posterior[k];
		    }			
	    }
	} //end main 		
 
    //return R list of length 3 with list("r1","r0","expected") 
    Rcpp::List ret;
    ret["r1"] = r1;
    ret["r0"] = r0;
    ret["expected"] = expected;
    return(ret);
    
    END_RCPP
}


//Estep for bfactor
RcppExport SEXP Estepbfactor(SEXP Ritemtrace, SEXP Rprior, 
    SEXP RX, SEXP Rnfact, SEXP Rr, SEXP Rsitems) 
{
    BEGIN_RCPP

    /*
        Ritemtrace = numeric matrix. Probability traces
        Rprior = numeric vector. Normalized prior
        RX = integer matrix. Response vector with 9's as missing
        Rnfact = integer. Number of factors
        Rr = integer vector. Counts of same response vector
        Rsitems = integer matrix. Specific factor indicator
    */

    Rcpp::NumericVector prior(Rprior);
    Rcpp::IntegerVector tmpnfact(Rnfact);
    Rcpp::IntegerVector r(Rr);
    Rcpp::IntegerMatrix data(RX);
    Rcpp::NumericMatrix itemtrace(Ritemtrace);
    Rcpp::IntegerMatrix sitems(Rsitems);
    int nfact = tmpnfact[0];
    int sfact = nfact - 1;
    int nitems = data.ncol();
    int npquad = prior.length();
    int nquad = npquad * npquad;  
    int npat = r.length();      
    int i=0, j=0, k=0, item=0, fact=0;	
    Rcpp::NumericMatrix r1(nitems*sfact,nquad);    
    Rcpp::NumericMatrix r0(nitems*sfact,nquad);
    Rcpp::NumericVector expected(npat);
        
	//declare dependent arrays 
	double likelihoods[sfact][nquad], L[npquad][npquad], tempsum[npquad],
		Plk[sfact][npquad], Elk[sfact][npquad], Pls[npquad], posterior[sfact][nquad];	
				
	// Begin main function body here				
	for (int pat = 0; pat < npat; pat++){		
		for (fact = 0; fact < sfact; fact++){ 	
			for (k = 0; k < nquad; k++)
				likelihoods[fact][k] = 1;				
			for (item = 0; item < nitems; item++){
				if (data(pat,item) == 9) 
					continue;
				if (data(pat,item)){	
					if(sitems(fact,item))
					    for (k = 0; k < nquad; k++)
					  	    likelihoods[fact][k] = likelihoods[fact][k] * itemtrace(item,k);					
				} else {
					if (sitems(fact,item))
					    for (k = 0; k < nquad; k++)
						    likelihoods[fact][k] = likelihoods[fact][k] * (1.0 - itemtrace(item,k));
				}			
			}
		}			
		for (fact = 0; fact < sfact; fact++){
			k = 0;
			for (j = 0; j < npquad; j++){
				tempsum[j] = 0.0;
			    for (i = 0; i < npquad; i++){
			  	    L[i][j] = likelihoods[fact][k];
			  	    k++;
			    }
			}
			for (j = 0; j < npquad; j++)				
			    for (i = 0; i < npquad; i++)
			  	    L[i][j] = L[i][j] * prior[j];
			for (j = 0; j < npquad; j++)				
		        for (i = 0; i < npquad; i++)
			        tempsum[j] += L[j][i];
			for (i = 0; i < npquad; i++)
			    Plk[fact][i] = tempsum[i];    			
		}				
		expected[pat] = 0.0;
		for (i = 0; i < npquad; i++){
		    Pls[i] = 1.0; 		  		
			for(fact = 0; fact < sfact; fact++)
			    Pls[i] = Pls[i] * Plk[fact][i];			
			expected[pat] += Pls[i] * prior[i];  
		}				
		for (fact = 0; fact < sfact; fact++)
		    for (i = 0; i < npquad; i++)
		  	    Elk[fact][i] = Pls[i] / Plk[fact][i];		  	
		for (fact = 0; fact < sfact; fact++)
		    for (i = 0; i < nquad; i++)  			  	
		        posterior[fact][i] = likelihoods[fact][i] * r[i] * Elk[fact][i % npquad] / expected[pat];
		for (fact = 0; fact < sfact; fact++){			
			for (item = 0; item < nitems; item++){
				if (data(pat,item) == 9) 
					continue;
				if (data(pat,item)){
					for (k = 0; k < nquad; k++)
						r1(item + nitems*fact,k) += posterior[fact][k];
				} else {
					for (k = 0; k < nquad; k++)
						r0(item + nitems*fact,k) += posterior[fact][k];
				}			
			}
		}	
	}	//end main 
	
    //return R list of length 3 with list("r1","r0","expected") 
    Rcpp::List ret;
    ret["r1"] = r1;
    ret["r0"] = r0;
    ret["expected"] = expected;
    return(ret);

    END_RCPP
}

