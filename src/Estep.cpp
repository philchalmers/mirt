#include <Rcpp.h>
using namespace Rcpp;

//Estep for mirt
RcppExport SEXP Estep(SEXP Ritemtrace, SEXP Rprior, SEXP RX,  
	SEXP Rnfact, SEXP Rr) 
{
    BEGIN_RCPP

    /*
        Ritemtrace = numeric matrix. Probability traces
        Rprior = numeric vector. Normalized prior
        RX = integer matrix.
        Rnfact = integer. Number of factors
        Rr = integer vector.  Counts of same response vector
    */
    
    NumericVector prior(Rprior);
    NumericVector log_prior = prior;
    IntegerVector nfact(Rnfact);
    IntegerVector r(Rr);
    IntegerMatrix data(RX);
    NumericMatrix itemtrace(Ritemtrace);
    NumericMatrix log_itemtrace = itemtrace;    
    int nquad = prior.length();
    int nitems = data.ncol();
    int npat = r.length();      
    double expd = 0.0;
    int i = 0, k = 0, item = 0;    
    NumericMatrix r1(nquad, nitems);    
    NumericVector expected(npat), posterior(nquad);
    List ret;
    for (item = 0; item < nitems; item++)
        for (k = 0; k < nquad; k++)
            log_itemtrace(k,item) = log(itemtrace(k,item));
    for (k = 0; k < nquad; k++)
        log_prior(k) = log(prior(k));
	
    // Begin main function body 				
	for (int pat = 0; pat < npat; pat++){		
  	    for (k = 0; k < nquad; k++)
		    posterior(k) = 0.0;
			
		for (item = 0; item < nitems; item++){
            if (data(pat,item))
                for (k = 0; k < nquad; k++)
			        posterior(k) += log_itemtrace(k,item);
		}    
        for (i = 0; i < nquad; i++)
            posterior(i) = exp(log_prior(i) + posterior(i));
	    expd = 0.0;
	    for (i = 0; i < nquad; i++)
	        expd += posterior(i);	
	    expected(pat) = expd;		
		
	    for (i = 0; i < nquad; i++)
	        posterior(i) = r(pat) * posterior(i) / expd;	
        for (item = 0; item < nitems; item++){              
            if (data(pat,item))
	            for (k = 0; k < nquad; k++)
	                r1(k,item) += posterior(k);
		    			
	    }
	} //end main  		
 
    //return R list of length 3 with list("r1","r0","expected") 
    ret["r1"] = r1;
    ret["expected"] = expected;
    return(ret);
    
    END_RCPP
}


//Estep for bfactor
RcppExport SEXP Estepbfactor(SEXP Ritemtrace, SEXP Rprior, SEXP RX, SEXP Rr, SEXP Rsitems) 
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
    NumericMatrix itemtrace(Ritemtrace);
    NumericMatrix log_itemtrace = itemtrace; 
    NumericVector prior(Rprior);
    IntegerVector r(Rr);
    IntegerMatrix data(RX);
    IntegerMatrix sitems(Rsitems);
    int nfact = sitems.ncol() + 1;
    int sfact = nfact - 1;
    int nitems = data.ncol();
    int npquad = prior.length();
    int nquad = npquad * npquad;  
    int npat = r.length();      
    int i=0, j=0, k=0, item=0, fact=0;
    for (item = 0; item < nitems; item++)
        for (k = 0; k < nquad; k++)
            log_itemtrace(k,item) = log(itemtrace(k,item));
        
	//declare dependent arrays 
	NumericVector tempsum(npquad), expected(npat), Pls(npquad);
	NumericMatrix likelihoods(nquad,sfact), L(npquad,npquad), r1(nquad,nitems*sfact),
		Plk(npquad,sfact), Elk(npquad,sfact), posterior(nquad,sfact);	
				
	// Begin main function body here				
	for (int pat = 0; pat < npat; pat++){		
		for (fact = 0; fact < sfact; fact++){ 	
			for (k = 0; k < nquad; k++)
				likelihoods(k,fact) = 0;				
			for (item = 0; item < nitems; item++){
				if (data(pat,item))	
					if(sitems(item,fact))
					    for (k = 0; k < nquad; k++)
					  	    likelihoods(k,fact) += log_itemtrace(k,item);					
			}
		} 
        for (fact = 0; fact < sfact; fact++)
            for (k = 0; k < nquad; k++)
		        likelihoods(k,fact) = exp(likelihoods(k,fact));					
		for (fact = 0; fact < sfact; fact++){
			k = 0;
			for (j = 0; j < npquad; j++){
				tempsum(j) = 0.0;
			    for (i = 0; i < npquad; i++){
			  	    L(i,j) = likelihoods(k,fact);
			  	    k++;
			    }
			}
			for (j = 0; j < npquad; j++)				
			    for (i = 0; i < npquad; i++)
			  	    L(i,j) = L(i,j) * prior(j);
			for (j = 0; j < npquad; j++)				
		        for (i = 0; i < npquad; i++)
			        tempsum(j) += L(j,i);
			for (i = 0; i < npquad; i++)
			    Plk(i,fact) = tempsum(i);    			
		}				
		expected(pat) = 0.0;
		for (i = 0; i < npquad; i++){
		    Pls(i) = 1.0; 		  		
			for(fact = 0; fact < sfact; fact++)
			    Pls(i) = Pls(i) * Plk(i,fact);			
			expected(pat) += Pls(i) * prior(i);  
		}				
		for (fact = 0; fact < sfact; fact++)
		    for (i = 0; i < npquad; i++)
		  	    Elk(i,fact) = Pls(i) / Plk(i,fact);		  	
		for (fact = 0; fact < sfact; fact++)
		    for (i = 0; i < nquad; i++)  			  	
		        posterior(i,fact) = likelihoods(i,fact) * r(pat) * Elk(i % npquad,fact) / expected(pat);
		// ordered specific factor packets, each the size of itemtrace
		for (fact = 0; fact < sfact; fact++){			
			for (item = 0; item < nitems; item++)
				if (data(pat,item))
					for (k = 0; k < nquad; k++)
						r1(k,item + nitems*fact) += posterior(k,fact);
			
		}	
	}	//end main 
	
    //return R list of length 3 with list("r1","r0","expected") 
    List ret;
    ret["r1"] = r1;
    ret["expected"] = expected;
    return(ret);

    END_RCPP
}

