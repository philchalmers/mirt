#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP fullScores(SEXP Rfulldata, SEXP Rtabdata, SEXP Rscores)
{	
    BEGIN_RCPP
	int i, j, pat, nf, nfact, ncol, N;
    NumericMatrix fulldata(Rfulldata);    
    NumericMatrix tabdata(Rtabdata);    
    NumericMatrix scores(Rscores);    
    N = fulldata.nrow();
    ncol = fulldata.ncol();
    nfact = scores.ncol();
	NumericMatrix ret(N,nfact);			
	IntegerVector val(N);

    for(pat = 0; pat < tabdata.nrow(); pat++){
	    for(i = 0; i < N; i++){
            val(i) = 0; 
            for(j = 0; j < ncol; j++)
                val(i) += (fulldata(i,j) == tabdata(pat,j));
            if(val(i) == ncol)
                for(nf = 0; nf < nfact; nf++)
                    ret(i, nf) = scores(pat, nf); 
        } 
    }
	return(ret);
	END_RCPP
}


