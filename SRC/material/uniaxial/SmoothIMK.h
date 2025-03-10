/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.0 $                                                                     
// Written: 
// IMK material with exponential transition curves and multiple cyclic modes (Bilinear, poinched and peak-oriented)


#ifndef SmoothIMK_h
#define SmoothIMK_h

#include <UniaxialMaterial.h>

class SmoothIMK : public UniaxialMaterial
{
public:
	SmoothIMK(int tag,
		double pd1, double pf1,
		double pd2, double pf2,
		double pd3, double pf3,
		double nd1, double nf1,
		double nd2, double nf2,
		double nd3, double nf3,
		double gama, double _c,
		double r0, double r1, double r2,
		int cyclicRule,
		double sigInit = 0.0);

	SmoothIMK(void);
	virtual ~SmoothIMK();


	const char* getClassType(void) const { return "SmoothIMK"; };

	double getInitialTangent(void);
	UniaxialMaterial* getCopy(void);

	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag = 0);

	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);

	virtual double getEnergy() { return EnergyP; };

	double getInitYieldStrain() { return Fy / E0; }
	virtual void resetEnergy(void) { EnergyP = 0; }

protected:

private:
	double pf1, pd1, pf2, pd2, pf3, pd3;
	double nf1, nd1, nf2, nd2, nf3, nd3;
	double r0, r1, r2;
	int cyclicRule;  // 1:bilinear, 2:pinched, 3:peak-oriented
	double FailEnerg, c;			//damage parameters
	double FydP, FydN;		//Pos and Neg Fy's affected by damage
	double ExcurEnergy;
	void updateDamage();
	double EnergyP; //by SAJalali
	double sigini; // initial 
	
	//HISTORY VARIABLES
	double epsminP; //  = hstvP(1) : max eps in compression
	double epsmaxP; //  = hstvP(2) : max eps in tension
	double epsplP;  //  = hstvP(3) : plastic excursion
	double epss0P;  //  = hstvP(4) : eps at asymptotes intersection
	double sigs0P;  //  = hstvP(5) : sig at asymptotes intersection
	double epssrP;  //  = hstvP(6) : eps at last inversion point
	double sigsrP;  //  = hstvP(7) : sig at last inversion point
	int    statP;    //  = hstvP(8) : index for loading/unloading
	// hstv : STEEL HISTORY VARIABLES   
	double epsP;  //  = strain at previous converged step
	double sigP;  //  = stress at previous converged step
	double eP;    //   stiffness modulus at last converged step;

	double epsmin;
	double epsmax;
	double epspl;
	double epss0;
	double sigs0;
	double epsr;
	double sigr;
	int    stat;
	double sig;
	double e;
	double eps;   //  = strain at current step
};


#endif

