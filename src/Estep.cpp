#include <Rcpp.h>
using namespace Rcpp;

//Estep for mirt
RcppExport SEXP Estep(SEXP Ritemtrace, SEXP Rprior, SEXP RX,  
	SEXP Rnfact, SEXP Rr) 
{
    BEGIN_RCPP

    NumericVector prior(Rprior);
    IntegerVector nfact(Rnfact);
    IntegerVector r(Rr);
    IntegerMatrix data(RX);
    NumericMatrix itemtrace(Ritemtrace);
    const int nquad = prior.length();
    const int nitems = data.ncol();
    const int npat = r.length();          
    NumericMatrix r1(nquad, nitems);    
    NumericVector expected(npat);
    NumericVector posterior(nquad); 
    List ret;
	
    // Begin main function body 				
	for (int pat = 0; pat < npat; ++pat){		  	           
        for(int q = 0; q < nquad; ++q)
            posterior(q) = prior(q);  
		for (int item = 0; item < nitems; ++item)
            if(data(pat,item))
                for(int q = 0; q < nquad; ++q)
                    posterior(q) = posterior(q) * itemtrace(q,item);
	    double expd = sum(posterior);
	    expected(pat) = expd;
        for(int q = 0; q < nquad; ++q)
	        posterior(q) = r(pat) * posterior(q) / expd;	
        for (int item = 0; item < nitems; ++item){              
            if (data(pat,item))	            
                for(int q = 0; q < nquad; ++q)
	                r1(q,item) = r1(q,item) + posterior(q);		    			
	    }
	} //end main  		
     
    ret["r1"] = r1;
    ret["expected"] = expected;
    return(ret);
    
    END_RCPP
}


//Estep for bfactor
RcppExport SEXP Estepbfactor(SEXP Ritemtrace, SEXP Rprior, SEXP RPriorbetween, SEXP RX, 
    SEXP Rr, SEXP Rsitems) 
{
    BEGIN_RCPP

    List ret;
    NumericMatrix itemtrace(Ritemtrace);    
    NumericVector prior(Rprior);
    NumericVector Priorbetween(RPriorbetween);
    IntegerVector r(Rr);
    IntegerMatrix data(RX);
    IntegerMatrix sitems(Rsitems);    
    const int sfact = sitems.ncol();
    const int nitems = data.ncol();
    const int npquad = prior.length();
    const int nbquad = Priorbetween.length();
    const int nquad = nbquad * npquad;  
    const int npat = r.length();
        
	NumericVector expected(npat), Pls(nbquad);
	NumericMatrix likelihoods(nquad,sfact), L(nbquad,npquad), r1(nquad,nitems*sfact),
		Plk(nbquad,sfact), Elk(nbquad,sfact), posterior(nquad,sfact);	
				
	// Begin main function body here				
	for (int pat = 0; pat < npat; ++pat){
        likelihoods.fill(1.0);
		for (int fact = 0; fact < sfact; ++fact){ 	
			for (int item = 0; item < nitems; ++item){
				if (data(pat,item) && sitems(item,fact))										    
				    for (int k = 0; k < nquad; ++k)
    				    likelihoods(k,fact) = likelihoods(k,fact)*itemtrace(k,item);					
			}
		}         					
		for (int fact = 0; fact < sfact; ++fact){
			int k = 0;
			for (int j = 0; j < npquad; ++j){				
			    for (int i = 0; i < nbquad; ++i){
			  	    L(i,j) = likelihoods(k,fact);
			  	    ++k;
			    }
			}
			for (int i = 0; i < nbquad; ++i)
                for (int q = 0; q < npquad; ++q)
			        L(i,q) = L(i,q) * prior(q);
            NumericVector tempsum(nbquad);
			for (int i = 0; i < npquad; ++i)
                for (int q = 0; q < nbquad; ++q)
			        tempsum(q) += L(q,i);
			for (int i = 0; i < nbquad; ++i)
			    Plk(i,fact) = tempsum(i);    			
		}
        Pls.fill(1.0);
		for (int i = 0; i < nbquad; ++i){		     		  		
			for(int fact = 0; fact < sfact; ++fact)
			    Pls(i) = Pls(i) * Plk(i,fact);			
			expected(pat) += Pls(i) * Priorbetween(i);  
		}				
		for (int fact = 0; fact < sfact; ++fact)
		    for (int i = 0; i < nbquad; ++i)
		  	    Elk(i,fact) = Pls(i) / Plk(i,fact);		  	
		for (int fact = 0; fact < sfact; ++fact)
		    for (int i = 0; i < nquad; ++i)  			  	
		        posterior(i,fact) = likelihoods(i,fact) * r(pat) * Elk(i % nbquad,fact) / expected(pat);		
        for (int item = 0; item < nitems; ++item)
    		if (data(pat,item))
		        for (int fact = 0; fact < sfact; ++fact)
                    for(int q = 0; q < nquad; ++q)
					    r1(q,item + nitems*fact) = r1(q,item + nitems*fact) + posterior(q,fact);
	}	//end main 
	    
    ret["r1"] = r1;
    ret["expected"] = expected;
    return(ret);

    END_RCPP
}

//EAP estimates used in multipleGroup
RcppExport SEXP EAPgroup(SEXP Ritemtrace, SEXP Rtabdata, SEXP RTheta, SEXP Rprior, SEXP Rmu) 
{
    BEGIN_RCPP

    NumericMatrix itemtrace(Ritemtrace); 
    NumericMatrix tabdata(Rtabdata); 
    NumericMatrix Theta(RTheta); 
    NumericVector prior(Rprior);
    NumericVector mu(Rmu);
    const int n = prior.length(); //nquad
    const int nitems = tabdata.ncol();
    const int nfact = Theta.ncol();

    NumericVector L(n), thetas(nfact), thetas2(nfact*(nfact+1)/2); 
    NumericMatrix scores(tabdata.nrow(), nfact), scores2(tabdata.nrow(), nfact*(nfact + 1)/2);
    double denom;

    for(int pat = 0; pat < tabdata.nrow(); ++pat){
        
        L.fill(1.0);
        for(int j = 0; j < n; ++j){            
            for(int i = 0; i < nitems; ++i)
                if(tabdata(pat, i))
                    L(j) = L(j) * itemtrace(j, i);             
        }
        
        thetas.fill(0.0);
        denom = 0.0;
        for(int j = 0; j < n; ++j)
            denom += (L(j)*prior(j));
        
        for(int k = 0; k < nfact; ++k){            
            for(int j = 0; j < n; ++j)
                thetas(k) += Theta(j, k) * L(j) * prior(j) / denom;
            scores(pat, k) = thetas(k);
        }

        int ind = 0;
        thetas2.fill(0.0);
        for(int i = 0; i < nfact; ++i){
            for(int k = 0; k < nfact; ++k){
                if(i <= k){
                    for(int j = 0; j < n; ++j)
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
