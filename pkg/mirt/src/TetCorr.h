/*
 *  TetCorr.h
 *
 *  Created by J.S.Fleming translated by Nic Jackson
 *
 */

#ifndef TetCorr_H
#define TetCorr_H

using namespace std;

extern "C" {
class TetCorr {
	
public:
	TetCorr(); // constuctor
	~TetCorr(); // destructor

	void Tetra(int T[2][2],double & R, double & SE, double & SE0, double & Z1, double & Z2, double & PHI_R, int & ITYPE, int & IFAULT);
	
	
private:
	
	/*
	 COMMON /TC/ A, B, C, D, TOT, ZERO, ONE, TWO, FOUR, SIX, HALF,
     &   TWOPI, SQT2PI, RLIMIT, BIGR, UPLIM, CONST, CHALF, CITER, CONV,
     &   NITER
	 */
	
	// need to change all these common things to const
	double A;
	double B;
	double C;
	double D;
	double TOT;
	
	double ZERO;
	double ONE;
	double TWO;
	double FOUR;
	double SIX;
	double HALF;
	double TWOPI;
	double SQT2PI;
	double RLIMIT;
	double BIGR;
	double UPLIM;
	
	double CONST;
	double CHALF;
	double CITER;
	double CONV;
	double NITER;
	
	double P2Z(double PROB, bool & ERR);
	double Z2P (double X, bool UPPER);
	
	void TETINIT (double & R0, bool LOFREQ, double AD, double BC, double Z1, double Z2);
	
	void GQ (double & R, double & RPREV, double & SUM, double PRAB, double PRAC, double Z1, double Z2, 
			 bool RNEG, int & ITYPE,double & ITER, int &IFAULT, double PHI0, bool & FINIS);
	
	void TSERIES (double & R, double PRA, double PRAB, double PRAC, double Z1, double Z2, int & ITYPE, double & ITER,int & IFAULT, double PHI0, bool & FINIS);
	
	/* compares the sign of A with B if A has a different sign to B it assigns B's sign to A
	 * i.e SIGN(1,-7) = -1;
	 *     SIGN(1,7) = 1;
	 */
	double SIGN(double A,double B); 
	
};
}

#endif

