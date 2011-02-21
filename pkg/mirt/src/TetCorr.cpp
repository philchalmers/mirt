/*
 *  TetCorr.cpp
 *
 *  Created by J.S.Fleming translated by Nic Jackson
 *
 */

/*
 *     Based on Morton Brown's algorithm to calculate tetrachoric
 *     correlations from a four-fold table T with cells:
 *
 *         _____________
 *        |      |      |
 *        |   A  |  B   | A + B
 *        |______|______|
 *        |      |      |
 *        |   C  |  D   | C + D
 *        |______|______|
 *
 *          A + C  B + D    N
 *
 *     The program was re-coded in structured and modularized code, and
 *     all calculations were made double precision. The Number of 
 *     iterations for both the series and quadrature was increased.
 *
 *     As another modification, Divgi's initial estimate replaces Yules'
 *     a the start value in all cases except when an observed cell 
 *     frequency is quite small. 
 *
 *     Subroutines used:
 *
 *        SUBROUTINE TETINIT:    Determines start value for correlation.
 *        SUBROUTINE TSERIES:    Computes tetrachoric series.
 *        SUBROUTINE GQ:         Performs Gaussian quadrature.
 *        REAL*8 FUNCTION Z2P:   Normal curve function.
 *        REAL*8 FUNCTION P2Z:   Inverse normal curve function.
 *
 *     Arguments:
 *
 *        T:       Four-fold (2 x 2) frequency table.
 *        R:       Resulting tetrachoric r.
 *        SE:      Standard error of r.
 *        SE0:     Standard error of r for test that rho = 0.
 *        ITYPE:   Method actually used:
 *                 1:  Gaussian quadrature
 *                 2:  cosine function (for z1 = z2 = 0)
 *                 3:  r is 1 or -1
 *                 4:  r is zero
 *                >5:  number of terms in the tetrachoric series
 *                <0:  as above, but all cells were modified by +/- 0.5
 *        IFAULT:  Termination code:
 *                 0:  normal termination
 *                 1:  iteration failed to converge
 *                 2:  table has at most one nonzero row/column or at
 *                     least one negative frequency.
 *
 *     Constants (defaults in parentheses):
 *
 *        CITER:   Convergence constant for iterations (10**-6).
 *        CONV:    Convergence consant for series convergence (10**-8)
 *        NITR:    Max iterations (99; modified from 25)
 *        BIGR:    Value of |r| that determines whether of not to use
 *                 quadrature (.95)
 *        UPLIM:   Upper limit c of the integral for quadrature (5)
 *        RLIMIT   Upper limit for  TrU when tetrachoric series or 
 *                 quadrature is used (.9999)
 *        CONST:   Constant used to rescale terms in series to avoid
 *                 underflow (10**-60)
 *        CHALF:   Used with CONST to rescale terms and must equal 
 *                 SQRT(CONST)
 *
 *     Sources:
 *
 *        Brown, M. B. (1977). Algorithm AS 116: The tetrachoric 
 *     correlation and its asymptotic standard error. Applied 
 *     statistics, 26, 343-351.
 *
 *        Fleming, J. S. (in press). TETCORR: A computer program 
 *     to compute smoothed tetrachoric correlation matrices. Behavior
 *     Research Methods, Instruments, and Computers.
 *
 */

#include "TetCorr.h"
#include <cmath>

TetCorr::TetCorr()  
{
	ZERO = 0.0;
	ONE = 1.0;
	TWO = 2.0;
	FOUR = 4.0;
	SIX = 6.0;
	HALF = 0.5;
	TWOPI = 6.28318531;
	SQT2PI = 2.50662827;
	RLIMIT = 0.9999;
	BIGR = 0.95;
	UPLIM = 5.0;
	CONST = 1.0e-60;
	CHALF = 1.0e-30;
	CONV = 1.0e-8;
	CITER = 1.0e-6;
	NITER = 99;	
}

TetCorr::~TetCorr() { }

void TetCorr::Tetra(int T[2][2],double & R, double & SE, double & SE0, double & Z1, double & Z2, double & PHI_R, int & ITYPE, int & IFAULT) {
	
	bool ERR = false;
	bool LOFREQ = false;
	bool FINIS = false;
	bool QUADR = false;
	bool RNEG = false;

	double POINT01 = 0.01;

	//Initializations.
	R = ZERO;
	SE = ZERO;
	SE0 = ZERO;
	
	IFAULT = 0;
	ITYPE = 0;
	int JTYPE = 1;
	double XINCR = ZERO;

	A = T[0][0];
	B = T[0][1];
	C = T[1][0];
	D = T[1][1];
	
	/*     Compute maginal and grand totals, and product terms for later 
	 *     computations.
	 */
	
	double R1 = A + B;
	double R2 = C + D;
	double C1 = A + C;
	double C2 = B + D;
	TOT = R1 + R2;
	
	double AD = A*D;
	double BC = B*C;
	
	/*
	 *     No marginal can be 0.
	 */
	if(R1 == ZERO || R2 == ZERO || C1 == ZERO || C2 == ZERO) {
		IFAULT = 2;
		return;
	}
	
	
	/*
	 *     Detect low cell frequencies.
	 */
	double FRAC = POINT01*TOT;
	
	if(A < FRAC || B < FRAC || C < FRAC || D < FRAC) {
		LOFREQ = true;		
	}else {
		LOFREQ = false;
	}
	
	/*
	 *     Determine sign of r.
	 */
	
	RNEG = (AD < BC);
	
	/*
	 *
	 *     Search for special cases.
	 *
	 */
	if(A == ZERO && D == ZERO) {
		// Case: r = -1.
		R = -ONE; // not sure if this works with C++
		ITYPE = 3;
		
	} else if(B == ZERO && D == ZERO) {
		// Case: r = +1.
		R = ONE;
		ITYPE = 3;
	} else if(AD == BC) {
		// Case: r=0.
		R = ZERO;
		ITYPE = 4;
	} else if(A == ZERO || D == ZERO) {
		// Case: Just one cell is zero (set increment for continuity correction
		XINCR = HALF;
		JTYPE = -1;
	} else if(B == ZERO || C == ZERO) {
		// Case: Just one cell is zero (set increment for continuity correction
		XINCR = HALF;
		JTYPE = -1;	
	}

	/*
	 *     Add increment to cells, recompute products. (Note: marginal
	 *     totals -- R1, R2, C1, C2 -- are not affected.)
	 */
	
	if(JTYPE == -1) {
		A = A + XINCR;
		B = B - XINCR;
		C = C - XINCR;
		D = D + XINCR;
		AD = A*D;
		BC = B*C;
	}
			
	/*
	 *     Compute phi coefficient. (Note - "phi" and "phi0" represent
	 *     bivariate PDFs and *     are not the same as the conventional
	 *     phi coefficient, here denoted "phi_r".)
     */
	double XNUM = AD - BC;
	double XDEN = sqrt(R1*R2*C1*C2);
	PHI_R = XNUM/XDEN;
	
	/*
	 *     Compute Z scores for variables 1 and 2, and the bivariate
	 *     function phi for r = 0: phi(z1, z2, r|r = 0).
	 */
	
	
	double PRA = 0;
	double PRAC = 0;
	
	if(!RNEG) {
		PRA  = A/TOT;
		PRAC = (A + C)/TOT;
	}else {
		PRA  = B/TOT;
		PRAC = (B + D)/TOT;
	}
	
	double PRAB = (A + B)/TOT;
	
	/*
	 *     [Note: Z1 corresponds to the column variable and Z2 is the row 
	 *     variable, keeping consistent with Brown's notation. This 
     *     more conventional matrix notation should be switched when 
	 *     printing to be consistent with conventional notation (rows = 
	 *     variable 1, cols = variable 2 for first table, etc.).]
	 */
	
	Z1 = P2Z (PRAC, ERR);
	Z2 = P2Z (PRAB, ERR);
	
	double PHI0 = exp(-HALF * (pow(Z1,2) + pow(Z2,2)) )/ TWOPI; // looks dodgy
	
	/*
	 *
	 *     Case: Equal marginals (use cosine formula). Use original table
	 *     values (T) because these have not been adjusted.
	 */
	if(ITYPE == 0) {
	
		if( (T[0][0] == T[1][1]) && (T[0][1] == T[1][0]) )
			ITYPE = 2;
		
	}
	
	/*
	 *     If no special cases were identified (ITYPE = 0), proceed
	 *     with tetrachoric series--unless estimated r > .95, in
	 *     which case use Gaussian quadrature.
	 */
	
	double R0;
	double RX;
	
	if(ITYPE == 0) {
		// Get start value (R0) for tetrachoric r.
		
		TETINIT (R0, LOFREQ, AD, BC, Z1, Z2); // fugly should probably be changed to return value R0
		R0 = abs(R0);
		RX = R0;
		
	}
						   
	/*
	 *     Do estimation iterations here.
	 *
	 *     Note that we can switch from tetrachoric series to quadra-
	 *     ture if R gets bigger than BIGR, but the reverse isn't true
	 *     -- once in quadrature we stay there.
	 */
	double ITER = 0;
    double JTER = 0;
	double RPREV = 0;
	double SUM = 0;
	
	FINIS  = false;
	
    QUADR  = (R0 > BIGR);
	
	while(!FINIS) {
		
		if (QUADR){
			
			GQ (R0, RPREV, SUM, PRAB, PRAC, Z1, Z2, RNEG,ITYPE, JTER, IFAULT, PHI0, FINIS);
		
		} else {
						   
			TSERIES (R0, PRA, PRAB, PRAC, Z1, Z2,ITYPE, ITER, IFAULT, PHI0, FINIS);
		   
			/*
		     *           If no convergence, or if r > BIGR, or r = 0 after
		     *           calling series, try quadrature.
		     */
			
			if(IFAULT > 0 || R0 > BIGR || R0 == ZERO) {
				QUADR = true;
				IFAULT = 0;
				ITYPE = 1;
				FINIS = false;
				ITER = 0;
				R0 = RX;
			}
			
		}
		
		if(IFAULT != 0)
			return;
		
	}
	
	R = R0;
	
	//Correct sign of r. Compute standard errors.
	if(ITYPE != 3 && ITYPE != 4){
		
		if(RNEG) {
			R  = -R;
			Z1 = -Z1;
		}
		
		ITYPE    = ITYPE * JTYPE;
		double ZZ1      = Z1*Z1;
		double ZZ2      = Z2*Z2;
		double Z12      = Z1*Z2;
		double RTRM     = ONE - R*R;
		double ROOT     = sqrt(RTRM);
		double PHI      = exp(-HALF*(ZZ1 - TWO*R*Z12 + ZZ2)/RTRM)/ (TWOPI*ROOT);
		double PR1      = Z2P ((Z1 - R*Z2)/ROOT, ERR) - HALF;
		double PR2      = Z2P ((Z2 - R*Z1)/ROOT, ERR) - HALF;
		double AB       = A*B;
		double CD       = C*D;
		double AC       = A*C;
		double BD       = B*D;
		SE = (A + D) * (B + C) / FOUR + PR2 * PR2*C1*C2 + PR1*PR1*R1*R2 + TWO*PR1*PR2*(AD - BC)- PR2*(AB - CD) - PR1*(AC - BD);
		
		if(SE < ZERO) 
			SE = ZERO;
		SE = sqrt(SE) / (sqrt(TOT)* TOT * PHI);
	}
	
	SE0      = sqrt((R1*R2*C1*C2)/TOT)/(TOT*TOT*PHI0);
	
	if(R == ZERO)
		SE = ZERO;
	
}
							   
double TetCorr::P2Z(double PROB, bool & ERR) {
							   
	/*
   REAL*8 FUNCTION P2Z (PROB, ERR)
   *
   *     Computes inverse gaussian value (given probability PROB, return
   *     Z-value. Algorithm described in Kennedy & Gentle (1980), based on
   *     Odeh & Evans (1974).
   *
   *     Because the Z value associated with zero or one is - or + ,
   *     infinity these values of PROB should obviously be avoided. An 
   *     error is raised (ERR = .TRUE.) if this occurs, or if very small 
   *     (close to zero) orvery large (close to one) values of PROB are 
   *     input.
   *
   *     Sources: 
   *
   *        Kennedy, W. J., & Gentle, J. W. (1980). Statistical 
   *     computing. New York, Dekker.
   *
   *        Odeh, R. E., & Evans, J. O. (1974). Algorithm AS 70: 
   *     Percentage points of the normal distribution. Applied Statistics, 
   *     23, 96-97.
   *
   */
	
	double ReturnValue;
	
	double PLIM;
	double TEMP;
	double Y;
	
	double P[5] = {-0.322232431088,-1.0,-0.342242088547,-0.0204231210245,-0.453642210148e-4};
	double Q[5] = {0.0993484626060,0.588581570495,0.531103462366,0.103537752850,0.38560700634e-2};
	
	PLIM = 10.0e-20;

   
	TEMP = PROB;
	double PR = PROB;
	ReturnValue = ZERO;
	
	if(PR > HALF)
		PR = ONE - PR;
	
	if(PR < PLIM) {
		ERR = true; // this is not good C++ transcribed almost literally from FORTRAN
		return 0.0; // we have raised an error however we need to return something
	}else {
		ERR = false;
	}
	
	if(PR != HALF) {
		Y = sqrt(log(ONE / (PR * PR)));
		ReturnValue = Y + ( ( ( (Y * P[4] + P[3]) * Y + P[2]) * Y + P[1]) * Y + P[0]) / ( ( ( (Y * Q[4] + Q[3]) * Y + Q[2]) * Y + Q[1]) * Y + Q[0]);
	}
	
	if(TEMP < HALF) 
		ReturnValue = -ReturnValue;
	
	return ReturnValue;
}

void TetCorr::TETINIT (double & R0, bool LOFREQ, double AD, double BC, double Z1, double Z2) {
	
	/*
	 *     Compute starting value of tetrachoric correlation. The tetra-
	 *     choric series assumes R0 > 0, so the initial estimate's sign may 
	 *     have to be reversed.
	 */
	
	double C1 = 0.12454;
	double C2 = 0.27102;
	double C3 = 0.82281;
	double C4 = 1.03514;
	double C5 = 0.07557;
	double C6 = 0.51141;
	double C7 = 2.05793;
	double C8 = 0.79289;
	double C9 = 4.28981;
	double C10 = 3.30231;
	
	/*
	 *     Compute Yule's estimate only if at least one cell has a low 
	 *     frequency; otherwise comute Divgi's estimate (default).
	 */
	
	
	double Q1, Q2,RR,SGN;
	
	
	if(LOFREQ) {
		R0 = pow((sqrt(AD) - sqrt(BC) ),2) / abs(AD - BC);
	} else {
		
		if(abs(Z1) > abs(Z2)) {
			Q1 = abs(Z1);
			Q2 = abs(Z2);
		} else {
			Q1 = abs(Z2);
			Q2 = abs(Z1);
		}
		
		RR = AD/BC;
		SGN = SIGN(1.0, Z1) * SIGN(1.0, Z2); //Not sure this is correct
		double T1 = Q1 * Q1 + Q2 * Q2; // need to check this is not some dodgy FORTRAN global variable
		double T2 = (Q1 - Q2) * (Q1 - Q2);
		double S1 = sqrt(T1);
		double PI = TWOPI / TWO;
		double AA = HALF / (ONE + T1 * (C1 - C2 * ( ONE - (Q1/S1) ) ) );
		double BB = HALF / (ONE + T1 * (C3 - C4 * (Q2 / S1) ) );
		double CC = C5 * Q1 + T2 * ( (C6 / (Q1 + C7) ) - (C6 / Q1) );
		double DD = SGN * Q2 * (C8 + C9/ (ONE + C10 * Q1) );
		double ALPHA = AA + BB * ( -ONE + ONE / ( ONE + CC * pow(log(RR) - DD,2 ) ) ); // double check this line could have misplaced bracket
		R0 = abs(cos(PI / (ONE + pow(RR,ALPHA) ) ) );
		
		/*        If Divgi's estimate is too big, convergence may be a problem. 
		          Revert to Yules. */
		if(R0 > BIGR) {
			// Below two lines were commented out not sure why
			// R0 = (SQRT(AD) - SQRT(BC))**2/(ABS(AD - BC))
			// LOFREQ = .TRUE.
			R0 = BIGR;
		}
			
	}
		
	
	
}

void TetCorr::GQ (double & R, double & RPREV, double & SUM, double PRAB, double PRAC, double Z1, double Z2, bool RNEG, int & ITYPE,double & ITER, int &IFAULT, double PHI0, bool & FINIS) {
	
/*
*     Use Gaussian quadrature to estimate tetrachoric r based.
*     This program is called by TETRA.
*/
	double X[16] = {0.9972638618, 0.9856115115, 0.9647622556, 0.9349060759,0.8963211558, 0.8493676137, 0.7944837960, 0.7321821187,
			0.6630442669, 0.5877157572, 0.5068999089, 0.4213512761,0.3318686023, 0.2392873623, 0.1444719616, 0.0483076657};
	double W[16] = {0.0070186100, 0.0162743947, 0.0253920653, 0.0342738629,0.0428358980, 0.0509980593, 0.0586840935,0.0658222228,
			0.0723457941, 0.0781938958, 0.0833119242, 0.0876520930,0.0911738787, 0.0938443991, 0.0956387201, 0.0965400885};
	
	bool DONE = false;
	bool DUMVAR = false;
	if(ITER == 0) {
		SUM = PRAB  * PRAC; //possbile return var
		RPREV = ZERO; //possbile return var
	}
	
	double SUMPRV = PRAB - SUM;
	double PROB;
	
	if(RNEG)
		PROB = A/TOT;
	else
		PROB = B/TOT;
	
	ITYPE = 1;

	/*
	 *     Loop to find estimate of correlation. Computation of
	 *     integral (sum) by quadrature.
	 */
	
	while(!DONE) {
			
		double RRSQ = sqrt(pow(ONE - R,2));
		double AMID = HALF * (UPLIM + Z1);
		double XLEN = UPLIM - AMID;
		SUM = ZERO;
		
		for(int IQUAD=0; IQUAD < 16; IQUAD++) {
			double XLA = AMID + X[IQUAD] * XLEN;
			double XLB = AMID - X[IQUAD] * XLEN;
			
			/*
			 *           To avoid underflows, TEMPA and TEMPB are used.
			 */
			
			double TEMPA = (Z2 - R * XLA) / RRSQ;
			if (TEMPA >= -SIX) 
				SUM = SUM + W[IQUAD]* exp(-HALF* pow(XLA,2)) * Z2P(TEMPA, DUMVAR);
			
			double TEMPB = (Z2 - R * XLB) / RRSQ;
			
			if (TEMPB >= -SIX)
				SUM = SUM + W[IQUAD] * exp(-HALF * pow(XLB,2)) * Z2P(TEMPB, DUMVAR);
		}
					
		SUM = SUM * XLEN / SQT2PI;
					
		/*
		 *        Convergence check.
		 */ 
		ITER = ITER + 1;

		if (abs(PROB - SUM) <= CITER) {
			FINIS = true; // return var
			return;
		}else {
			ITER = ITER + 1;
			if (ITER > NITER){
				IFAULT = 1; // return var
				return;
			}
		}
					
		/*
   	     *        Estimate correlation for next iteration by linear 
		 *        interpolation.
		 */

		double REST = ((PROB - SUM) * RPREV - (PROB - SUMPRV) * R ) / (SUMPRV - SUM);
	
		//Is estimate positive and less than upper limit?

		if (REST > RLIMIT)
			REST = RLIMIT;
					
		if (REST < ZERO)
			REST = ZERO;
			
		RPREV = R;
		R = REST;
		SUMPRV = SUM;

		//        If estimate has same value on two iterations, stop.
		if (R == RPREV)
			DONE = true;
			
					
	}
	
	return;
}

double TetCorr::Z2P (double X, bool UPPER) {
	
	/*
	*     This algorithm is due to:
	*
	*        Hill, I. D. (1973). Algorithm AS 66: The normal integral.
	*     Applied Statistics, 22, 424-427.
	*
	*     Evaluates the tail area of the standard normal curve from x
	*     to infinity (if UPPER is true), or from -infinity to x (other-
																  *     wise).
	* 
	*     LTONE and UTZERO must be set to suit the needs of a particular
	*     computer (see article). These are for evaluting the function
	*     near the extremes for the lower and upper tails. (We go with the
	*     author's numbers here -- however, we expand to double precision,
	*     meaning that more accuracy should be possible for very extreme
	*     values.)
	*/
	
	bool UP;
	double LTONE = 7.0;
	double UTZERO = 18.66;

	double CON = 1.28;
	
	double ReturnVar;
	
	UP = UPPER;
	double Z = X;
	if(Z < 0) {
		Z = -Z; // ??????????
		UP = !UP;
	}
	
	if(Z < LTONE || ((UP < UTZERO) && (Z < UTZERO)) ) {
		
		double Y = HALF * Z * Z;
		
		if(Z < CON) {
			ReturnVar = HALF - Z * (0.398942280444 - 0.399903438504 * Y / (Y + 5.75885480458 - 29.8213557808 / 
					    (Y + 2.62433121679 + 48.6959930692 / (Y + 5.92885724438))));
		}else {
			ReturnVar = 0.398942280385 * exp(-Y) / (Z - 3.8052e-8 + 1.00000615302 / (Z + 3.98064794e-4 + 1.98615381364 /
					   (Z - 0.151679116635 + 5.29330324926 / (Z + 4.8485912808 - 15.1508972451 /
					   (Z + 0.742380924027 + 30.789933034 / (Z + 3.99019417011))))));
		}
		
	} else {
		ReturnVar = 0.0;
	}
	
	if(!UP)
		ReturnVar = ONE - ReturnVar;

	return ReturnVar;

}

void TetCorr::TSERIES(double & R, double PRA, double PRAB, double PRAC, double Z1, double Z2, int & ITYPE, double & ITER,int & IFAULT, double PHI0, bool & FINIS) {

/*
 *     Computes tetrachoric series. Called by SUBROUTINE TETRA.
 */
	

	bool DONE = false;	

	double VA = ONE;
	double VB = Z1;
	double WA = ONE;
	double WB = Z2;
	double TERM = ONE;
	double ITERM = 0;
	double SUM = PRAB*PRAC;
	double DERIV = ZERO;
	double SR = PHI0;

	// Begin tetrachoric series.

	while(!DONE) {
 
		if(abs(SR) <= CONST) {
			SR = SR / CONST;
			VA = VA * CHALF;
			VB = VB * CHALF;
			WA = WA * CHALF;
			WB = WB * CHALF;
		}
 
		//Form sum and derivative of series.

		double DR = SR * VA * WA;
		SR = SR * R / TERM;
		double COF = SR * VA * WA;

		//ITERM counts number of consecutive terms < CONV.
		
		ITERM = ITERM + 1;
	 
		if(abs(COF) > CONV)
			ITERM = 0;
		
		SUM = SUM + COF;
		DERIV = DERIV + DR;
		double VAA = VA;
		double WAA = WA;
		VA = VB;
		WA = WB;
		VB = (Z1 * VA) - (TERM * VAA);
		WB = (Z2 * WA) - (TERM * WAA);
		TERM = TERM + ONE;
		
		if(ITERM >= 2 && TERM > SIX)
			DONE = true;
	
	}
	
	//Check for convergence.

	if(abs(SUM - PRA) <= CITER) {

		// Iteration has converged, set ITYPE.
		ITYPE = TERM;
		FINIS = true;
		return;
	}
	

	// Calculate next estimate of correlation.

	ITER = ITER + 1;

    // If too many iterations, terminate.
	if(ITER >= NITER) {
		IFAULT = 1;
		return;
	}

	double DELTA = (SUM - PRA) / DERIV;
	//double RRPREV = R; THERE SEEMS TO BE NO POINT TO THIS LINE?
	R = R - DELTA;
	
	if(ITER == 1)
		R = R + HALF * DELTA;
	
	if(R > RLIMIT)
		R = RLIMIT;
	
	if(R < ZERO)
		R = ZERO;

	return;
}

double TetCorr::SIGN(double A,double B) {

	if(A >= 0 && B < 0)
		return -A;
	
	if(A < 0 && B >=0)
		return +A;
	
	return A;

}
