#include <Rcpp.h>
using namespace Rcpp;

//Estep for mirt
RcppExport SEXP Estep(SEXP Ritemtrace, SEXP Rprior, SEXP RX,  
	SEXP Rnfact, SEXP Rr) 
{
    BEGIN_RCPP

    NumericVector prior(Rprior);
    NumericVector log_prior(prior.length());
    IntegerVector nfact(Rnfact);
    IntegerVector r(Rr);
    IntegerMatrix data(RX);
    NumericMatrix itemtrace(Ritemtrace);
    NumericMatrix log_itemtrace(itemtrace.nrow(), itemtrace.ncol());    
    const int nquad = prior.length();
    const int nitems = data.ncol();
    const int npat = r.length();          
    NumericMatrix r1(nquad, nitems);    
    NumericVector expected(npat);
    List ret;
    for (int item = 0; item < nitems; item++)        
        log_itemtrace(_,item) = log(itemtrace(_,item));
    log_prior = log(prior);
	
    // Begin main function body 				
	for (int pat = 0; pat < npat; pat++){		
  	    NumericVector posterior(nquad);
		for (int item = 0; item < nitems; item++)
            posterior += data(pat,item) * log_itemtrace(_,item);
        posterior = exp(log_prior + posterior);
	    double expd = sum(posterior);
	    expected(pat) = expd;
	    posterior = r(pat) * posterior / expd;	
        for (int item = 0; item < nitems; item++){              
            if (data(pat,item))
	            for (int k = 0; k < nquad; k++)
	                r1(k,item) += posterior(k);		    			
	    }
	} //end main  		
     
    ret["r1"] = r1;
    ret["expected"] = expected;
    return(ret);
    
    END_RCPP
}


//Estep for bfactor
RcppExport SEXP Estepbfactor(SEXP Ritemtrace, SEXP Rprior, SEXP RX, SEXP Rr, SEXP Rsitems) 
{
    BEGIN_RCPP

    NumericMatrix itemtrace(Ritemtrace);
    NumericMatrix log_itemtrace(itemtrace.nrow(), itemtrace.ncol()); 
    NumericVector prior(Rprior);
    NumericVector log_prior(prior.length());
    IntegerVector r(Rr);
    IntegerMatrix data(RX);
    IntegerMatrix sitems(Rsitems);
    const int nfact = sitems.ncol() + 1;
    const int sfact = nfact - 1;
    const int nitems = data.ncol();
    const int npquad = prior.length();
    const int nquad = npquad * npquad;  
    const int npat = r.length();      
    int i=0, j=0, k=0, item=0, fact=0;
    for (int item = 0; item < nitems; item++)        
        log_itemtrace(_,item) = log(itemtrace(_,item));
    log_prior = log(prior);
        
	NumericVector tempsum(npquad), expected(npat), Pls(npquad);
	NumericMatrix likelihoods(nquad,sfact), L(npquad,npquad), r1(nquad,nitems*sfact),
		Plk(npquad,sfact), Elk(npquad,sfact), posterior(nquad,sfact);	
				
	// Begin main function body here				
	for (int pat = 0; pat < npat; pat++){		
		for (fact = 0; fact < sfact; fact++){ 	
			for (k = 0; k < nquad; k++)
    			likelihoods(k,fact) = 0.0;
			for (item = 0; item < nitems; item++){
				if (data(pat,item) && sitems(item,fact))										    
				    for (k = 0; k < nquad; k++)
    				    likelihoods(k,fact) += log_itemtrace(k,item);					
			}
		}         					
		for (fact = 0; fact < sfact; fact++){
			k = 0;
			for (j = 0; j < npquad; j++){
				tempsum(j) = 0.0;
			    for (i = 0; i < npquad; i++){
			  	    L(i,j) = likelihoods(k,fact);
			  	    k++;
			    }
			}
			for (i = 0; i < npquad; i++)
			    L(i,_) = exp(L(i,_) + log_prior);
			for (i = 0; i < npquad; i++)
			        tempsum += L(_,i);
			for (i = 0; i < npquad; i++)
			    Plk(i,fact) = tempsum(i);    			
		}		
        for (fact = 0; fact < sfact; fact++)
            for (k = 0; k < nquad; k++)
    	        likelihoods(k,fact) = exp(likelihoods(k,fact));
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
	    
    List ret;
    ret["r1"] = r1;
    ret["expected"] = expected;
    return(ret);

    END_RCPP
}

//EAP estimates used in multipleGroup
RcppExport SEXP EAPgroup(SEXP Rlog_itemtrace, SEXP Rtabdata, SEXP RTheta, SEXP Rprior, SEXP Rmu) 
{
    BEGIN_RCPP

    NumericMatrix log_itemtrace(Rlog_itemtrace); 
    NumericMatrix tabdata(Rtabdata); 
    NumericMatrix Theta(RTheta); 
    NumericVector prior(Rprior);
    NumericVector mu(Rmu);
    int i, j, k, ind;
    const int n = prior.length(); //nquad
    const int nitems = tabdata.ncol();
    const int nfact = Theta.ncol();

    NumericVector L(n), thetas(nfact), thetas2(nfact*(nfact+1)/2); 
    NumericMatrix scores(tabdata.nrow(), nfact), scores2(tabdata.nrow(), nfact*(nfact + 1)/2);
    double LL, denom;

    for(int pat = 0; pat < tabdata.nrow(); pat++){
        
        for(j = 0; j < n; j++){
            LL = 0.0;
            for(i = 0; i < nitems; i++)
               LL += tabdata(pat, i) * log_itemtrace(j, i); 
            L(j) = exp(LL);
        }
        
        thetas.fill(0.0);
        denom = 0.0;
        for(j = 0; j < n; j++)
            denom += (L(j)*prior(j));
        
        for(k = 0; k < nfact; k++){            
            for(j = 0; j < n; j++)
                thetas(k) += Theta(j, k) * L(j) * prior(j) / denom;
            scores(pat, k) = thetas(k);
        }

        ind = 0;
        thetas2.fill(0.0);
        for(i = 0; i < nfact; i++){
            for(k = 0; k < nfact; k++){
                if(i <= k){
                    for(j = 0; j < n; j++)
                        thetas2(ind) += (Theta(j, i) - thetas(i)) * (Theta(j, k) - thetas(k)) *
                            L(j) * prior(j) / denom;
                    thetas2(ind) += (thetas(i) - mu(i)) * (thetas(k) - mu(k));
                    scores2(pat, ind) = thetas2(ind);
                    ind += 1;
                }
            }
        }
    }

    List ret;    
    ret["scores"] = scores;
    ret["scores2"] = scores2;
    return(ret);

    END_RCPP
}



