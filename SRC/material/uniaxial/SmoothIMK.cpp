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

// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SmoothIMK.cpp,v $

// Written: S. A. Jalali 10/2019
// Adding Cyclic and in-cycle deterioration modes to steel02 UniaxialMaterial

#include <math.h>

#include <stdlib.h>
#include <SmoothIMK.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <elementAPI.h>
#include <OPS_Globals.h>
#include "classTags.h"

static int numThisCall = 0;

void*
OPS_SmoothIMK()
{
	if (numThisCall == 0) {
		opserr << "------ SmoothIMK unaxialMaterial -------\n";
		opserr << "-------Syntax:\n";
		opserr << "-------UniaxialMaterial SmoothIMK $matTag \n";
		numThisCall = 1;
	}
	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;

	int numData = 1;
	int tag = 0;
	if (OPS_GetIntInput(&numData, &tag) != 0) {
		opserr << "WARNING invalid uniaxialMaterial SmoothIMK tag" << endln;
		return 0;
	}

	numData = OPS_GetNumRemainingInputArgs();

	if (numData != 8 && numData != 11 && numData != 15 && numData != 16) {
		opserr << "Invalid SmoothIMK #args for: " << tag << " see the syntax" << endln;
		return 0;
	}
	//default parameters:
	numData = 1;
	double pf1, pd1, pf2, pd2, pf3, pd3;
	double nf1, nd1, nf2, nd2, nf3, nd3;
	double gama, c, r0, r1, r2, sigInit;
	int rule;

	if (OPS_GetDoubleInput(&numData, &pd1) != 0) {
		opserr << "SmoothIMK:: invalid pd1 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &pf1) != 0) {
		opserr << "SmoothIMK:: invalid pf1 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &pd2) != 0) {
		opserr << "SmoothIMK:: invalid pd2 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &pf2) != 0) {
		opserr << "SmoothIMK:: invalid pf2 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &pd3) != 0) {
		opserr << "SmoothIMK:: invalid pd3 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &pf3) != 0) {
		opserr << "SmoothIMK:: invalid pf3 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &nd1) != 0) {
		opserr << "SmoothIMK:: invalid nd1 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &nf1) != 0) {
		opserr << "SmoothIMK:: invalid nf1 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &nd2) != 0) {
		opserr << "SmoothIMK:: invalid nd2 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &nf2) != 0) {
		opserr << "SmoothIMK:: invalid nf2 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &nd3) != 0) {
		opserr << "SmoothIMK:: invalid nd3 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &nf3) != 0) {
		opserr << "SmoothIMK:: invalid nf3 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &gama) != 0) {
		opserr << "SmoothIMK:: invalid gama for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &c) != 0) {
		opserr << "SmoothIMK:: invalid c for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &r0) != 0) {
		opserr << "SmoothIMK:: invalid r0 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &r1) != 0) {
		opserr << "SmoothIMK:: invalid r1 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &r2) != 0) {
		opserr << "SmoothIMK:: invalid r2 for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetIntInput(&numData, &rule) != 0) {
		opserr << "SmoothIMK:: invalid cyclicRule for material : " << tag << endln;
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &sigInit) != 0) {
		opserr << "SmoothIMK:: invalid sigInit for material : " << tag << endln;
		return 0;
	}

	// Parsing was successful, allocate the material
	theMaterial = new SmoothIMK(tag, pd1, pf1, pd2, pf2, pd3, pf3, nd1, nd2, nf1, nf2, nd3, nf3, gama, c, r0, r1, r2, rule, sigInit);


	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type SmoothIMK Material\n";
		return 0;
	}

	return theMaterial;
}


SmoothIMK::SmoothIMK(int tag,
	double _pd1, double _pf1,
	double _pd2, double _pf2,
	double _pd3, double _pf3,
	double _nd1, double _nf1,
	double _nd2, double _nf2,
	double _nd3, double _nf3,
	double _gama, double _c,
	double _r0, double _r1, double _r2,
	int _cyclicRule,
	double sigInit) :
	UniaxialMaterial(tag, MAT_TAG_SmoothIMK),
	pd1(_pd1), pf1(_pf1),
	pd2(_pd2), pf2(_pf2),
	pd3(_pd3), pf3(_pf3),
	nd1(_nd1), nf1(_nf1),
	nd2(_nd2), nf2(_nf2),
	nd3(_nd3), nf3(_nf3),
	c(_c), FailEnerg(_gama* pf1* pd1),
	r0(_r0), r1(_r1), r2(_r2),
	cyclicRule(_cyclicRule),
	sigini(sigInit)
{
	EnergyP = 0;
	eP = pf1 / pd1;
	epsP = sigini / (pf1 / pd1);
	sigP = sigini;
	sig = 0.0;
	eps = 0.0;
	e = pf1 / pd1;

	statP = 0;
	epsmaxP = pd1;
	epsminP = nd1;
	epsplP = 0.0;
	epss0P = 0.0;
	sigs0P = 0.0;
	epssrP = 0.0;
	sigsrP = 0.0;
	ExcurEnergy = 0;
	FydP = pf1;
	FydN = nf1;
}

SmoothIMK::SmoothIMK(void) :
	UniaxialMaterial(0, MAT_TAG_SmoothIMK),
	pd1(0), pf1(0),
	pd2(0), pf2(0),
	pd3(0), pf3(0),
	nd1(0), nf1(0),
	nd2(0), nf2(0),
	nd3(0), nf3(0),
	c(0), FailEnerg(0),
	r0(0), r1(0), r2(0),
	cyclicRule(1),
	sigini(0)
{
	EnergyP = 0;
	eP = 0;
	epsP = 0;
	sigP = 0;
	sig = 0.0;
	eps = 0.0;
	e = 0;

	statP = 0;
	epsmaxP = 0;
	epsminP = 0;
	epsplP = 0.0;
	epss0P = 0.0;
	sigs0P = 0.0;
	epssrP = 0.0;
	sigsrP = 0.0;

	ExcurEnergy = 0;
	FydP = FydN = 0;
}

SmoothIMK::~SmoothIMK(void)
{
	// Does nothing
}

UniaxialMaterial*
SmoothIMK::getCopy(void)
{
	SmoothIMK* theCopy = new SmoothIMK(this->getTag(), pd1, pf1, pd2, pf2, pd3,
		pf3, nd1, nd2, nf1, nf2, nd3, nf3, FailEnerg / pf1 / pd1, c, r0, r1, r2, cyclicRule, sigini);

	return theCopy;
}

double
SmoothIMK::getInitialTangent(void)
{
	return pf1/pd1;
}

int
SmoothIMK::setTrialStrain(double trialStrain, double strainRate)
{
	double EshP = (pf2-pf1)/(pd2-pd1);
	double EshN = (nf2-nf1)/(nd2-nd1);
	double E0p = pf1 / pd1;
	double E0n = nf1 / nd1;
	eps = trialStrain + sigini > 0 ? sigini/ E0p : sigini/E0n;
	const double deps = eps - epsP;
	if (fabs(deps) < 1.e-20)
		return 0;
	double b2p = EshP/E0p;
	double b2n = EshN/E0n;

	epsmax = epsmaxP;
	epsmin = epsminP;
	epspl = epsplP;
	epss0 = epss0P;
	sigs0 = sigs0P;
	epsr = epssrP;
	sigr = sigsrP;
	stat = statP;

	if (stat == 0) {
		epsmax = pd1;
		epsmin = nd1;
		if (deps < 0.0) {
			stat = 2;
			epss0 = epsmin;
			sigs0 = FydN;
			epspl = epsmin;
		}
		else
		{
			stat = 1;
			epss0 = epsmax;
			sigs0 = FydP;
			epspl = epsmax;
		}
	}

	if (stat == 2) {
		if (deps > 0.0)
		{
			// in case of load reversal from negative to positive strain increment, 
			// update the minimum previous strain, store the last load reversal 
			// point (epspl) and calculate the stress and strain (sigs0 and epss0) at the 
			// new intersection between elastic and strain hardening asymptote 
			stat = 1;
			epsr = epsP;
			sigr = sigP;
			if (epsP < epsmin)
				epsmin = epsP;
			const double d1 = (epsmax - epsmin) / (2.0 * pd1);
			epss0 = (FydP - EshP * pd1 - sigr + E0p * epsr) / (E0p - EshP);
			sigs0 = FydP  + EshP * (epss0 - pd1);
			epspl = epsmax;
		}
		else if (eps < (nd2 + nd1) / 2)
		{
			epss0 = nd2;
			sigs0 = FydN + EshN * (nd2 - nd1);
			epsr = epsP;
			sigr = sigP;
			stat = 4;
		}
	}
	else if (stat == 1) {
		if (deps < 0.0)
		{
			stat = 2;
			epsr = epsP;
			sigr = sigP;
			if (epsP > epsmax)
				epsmax = epsP;

			const double d1 = (epsmax - epsmin) / (2.0 * nd1);
			epss0 = (FydN + EshN * nd1 - sigr + E0n * epsr) / (E0n - EshN);
			sigs0 = FydN + EshN * (epss0 + nd1);
			epspl = epsmin;
		}
		else if (eps > (pd1 + pd2) / 2)
		{
			epss0 = pd2;
			sigs0 = FydP + EshP * (pd2 - pd1);
			epsr = epsP;
			sigr = sigP;
			stat = 3;
		}
	}
	if (stat == 4)
	{
		if (deps > 0.0)
		{
			stat = 1;
			epsr = epsP;
			sigr = sigP;
			if (epsP < epsmin)
				epsmin = epsP;
			const double d1 = (epsmax - epsmin) / (2.0 * (a4 * epsy));
			const double shft = 1.0 + a3 * pow(d1, 0.8);
			epss0 = (FydP * shft - Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
			sigs0 = FydP * shft + Esh * (epss0 - epsy * shft);
			epspl = epsmax;
		}
		else
			b2 = pstcpEFac / b;

	}
	if (stat == 3)
	{
		if (deps < 0.0)
		{
			stat = 2;
			epsr = epsP;
			sigr = sigP;
			if (epsP > epsmax)
				epsmax = epsP;

			const double d1 = (epsmax - epsmin) / (2.0 * (a2 * epsy));
			const double shft = 1.0 + a1 * pow(d1, 0.8);
			epss0 = (-FydN * shft + Esh * epsy * shft - sigr + E0 * epsr) / (E0 - Esh);
			sigs0 = -FydN * shft + Esh * (epss0 + epsy * shft);
			epspl = epsmin;
		}
		else
			b2 = pstcpEFac / b;

	}
	if (stat == 5)
	{
		sig = resFac * Fy;
		e = 1e-15 * E0;
	}
	else {

		double R = R0;
		// calculate current stress sig and tangent modulus E 
		if (cR1 != 0.0 && cR2 != 0.0)
		{
			const double xi = fabs((epspl - epss0) / epsy);
			R *= (1.0 - (cR1 * xi) / (cR2 + xi));
		}
		const double epsrat = (eps - epsr) / (epss0 - epsr);
		const double dum1 = 1.0 + pow(fabs(epsrat), R);
		const double dum2 = pow(dum1, (1 / R));

		sig = b2 * epsrat + (1.0 - b2) * epsrat / dum2;
		sig = sig * (sigs0 - sigr) + sigr;

		e = b2 + (1.0 - b2) / (dum1 * dum2);
		e *= (sigs0 - sigr) / (epss0 - epsr);
		if (stat == 3 && sig < resFac * Fy)
		{
			sig = resFac * Fy;
			e = 1e-15 * E0;
			stat = 5;
		}
		else if (stat == 4 && sig > -resFac * Fy)
		{
			sig = -resFac * Fy;
			e = 1e-15 * E0;
			stat = 5;
		}
	}
	return 0;
}

double
SmoothIMK::getStrain(void)
{
	return eps;
}

double
SmoothIMK::getStress(void)
{
	return sig;
}

double
SmoothIMK::getTangent(void)
{
	return e;
}

int
SmoothIMK::commitState(void)
{
	epsminP = epsmin;
	epsmaxP = epsmax;
	epsplP = epspl;
	epss0P = epss0;
	sigs0P = sigs0;
	epssrP = epsr;
	sigsrP = sigr;

	updateDamage();
	statP = stat;
	eP = e;
	sigP = sig;
	epsP = eps;

	return 0;
}

int
SmoothIMK::revertToLastCommit(void)
{
	epsmin = epsminP;
	epsmax = epsmaxP;
	epspl = epsplP;
	epss0 = epss0P;
	sigs0 = sigs0P;
	epsr = epssrP;
	sigr = sigsrP;
	stat = statP;

	e = eP;
	sig = sigP;
	eps = epsP;
	return 0;
}

int
SmoothIMK::revertToStart(void)
{
	EnergyP = 0;	//by SAJalali
	eP = E0;
	epsP = sigini / E0;
	sigP = sigini;
	sig = 0.0;
	eps = 0.0;
	e = E0;

	statP = 0;
	epsmaxP = Fy / E0;
	epsminP = -epsmaxP;
	epsplP = 0.0;
	epss0P = 0.0;
	sigs0P = 0.0;
	epssrP = 0.0;
	sigsrP = 0.0;


	ExcurEnergy = 0;
	FydP = FydN = Fy;
	return 0;
}

int
SmoothIMK::sendSelf(int commitTag, Channel& theChannel)
{
	static Vector data(33);
	data(0) = Fy;
	data(1) = E0;
	data(2) = b;
	data(3) = R0;
	data(4) = cR1;
	data(5) = cR2;
	data(6) = a1;
	data(7) = a2;
	data(8) = a3;
	data(9) = a4;
	data(10) = epsminP;
	data(11) = epsmaxP;
	data(12) = epsplP;
	data(13) = epss0P;
	data(14) = sigs0P;
	data(15) = epssrP;
	data(16) = sigsrP;
	data(17) = statP;
	data(18) = epsP;
	data(19) = sigP;
	data(20) = eP;
	data(21) = this->getTag();
	data(22) = sigini;
	data(23) = EnergyP;
	data(24) = epsPCFac;
	data(25) = pstcpEFac;
	data(26) = gama;
	data(27) = c;
	data(28) = resFac;
	data(29) = FydP;
	data(30) = FydN;
	data(31) = ExcurEnergy;
	data(32) = gama;
	if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "SmoothIMK::sendSelf() - failed to sendSelf\n";
		return -1;
	}
	return 0;
}

int
SmoothIMK::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	static Vector data(33);	//editted by SAJalali for energy

	if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "SmoothIMK::recvSelf() - failed to recvSelf\n";
		return -1;
	}

	Fy = data(0);
	E0 = data(1);
	b = data(2);
	R0 = data(3);
	cR1 = data(4);
	cR2 = data(5);
	a1 = data(6);
	a2 = data(7);
	a3 = data(8);
	a4 = data(9);
	epsminP = data(10);
	epsmaxP = data(11);
	epsplP = data(12);
	epss0P = data(13);
	sigs0P = data(14);
	epssrP = data(15);
	sigsrP = data(16);
	statP = int(data(17));
	epsP = data(18);
	sigP = data(19);
	eP = data(20);
	this->setTag(int(data(21)));
	sigini = data(22);
	EnergyP = data(23);
	epsPCFac = data(24);
	pstcpEFac = data(25);
	gama = data(26);
	c = data(27);
	resFac = data(28);
	FydP = data(29);
	FydN = data(30);
	ExcurEnergy = data(31);
	gama = data(32);
	FailEnerg = gama * Fy * Fy / E0;

	e = eP;
	sig = sigP;
	eps = epsP;

	return 0;
}

void
SmoothIMK::Print(OPS_Stream& s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		//    s << "SmoothIMK:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
		s << "SmoothIMK tag: " << this->getTag() << endln;
		s << "  fy: " << Fy << ", ";
		s << "  E0: " << E0 << ", ";
		s << "   b: " << b << ", ";
		s << "  R0: " << R0 << ", ";
		s << " cR1: " << cR1 << ", ";
		s << " cR2: " << cR2 << ", ";
		s << "  a1: " << a1 << ", ";
		s << "  a2: " << a2 << ", ";
		s << "  a3: " << a3 << ", ";
		s << "  a4: " << a4;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"SmoothIMK\", ";
		s << "\"E\": " << E0 << ", ";
		s << "\"fy\": " << Fy << ", ";
		s << "\"b\": " << b << ", ";
		s << "\"R0\": " << R0 << ", ";
		s << "\"cR1\": " << cR1 << ", ";
		s << "\"cR2\": " << cR2 << ", ";
		s << "\"a1\": " << a1 << ", ";
		s << "\"a2\": " << a2 << ", ";
		s << "\"a3\": " << a3 << ", ";
		s << "\"a4\": " << a4 << ", ";
		s << "\"sigini\": " << sigini << "}";
	}
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
SmoothIMK::setParameter(const char** argv, int argc, Parameter& param)
{

	if (strcmp(argv[0], "sigmaY") == 0 || strcmp(argv[0], "fy") == 0 || strcmp(argv[0], "Fy") == 0) {
		param.setValue(Fy);
		return param.addObject(1, this);
	}
	if (strcmp(argv[0], "E") == 0) {
		param.setValue(E0);
		return param.addObject(2, this);
	}
	if (strcmp(argv[0], "b") == 0) {
		param.setValue(b);
		return param.addObject(3, this);
	}
	if (strcmp(argv[0], "a1") == 0) {
		param.setValue(a1);
		return param.addObject(4, this);
	}
	if (strcmp(argv[0], "a2") == 0) {
		param.setValue(a2);
		return param.addObject(5, this);
	}
	if (strcmp(argv[0], "a3") == 0) {
		param.setValue(a3);
		return param.addObject(6, this);
	}
	if (strcmp(argv[0], "a4") == 0) {
		param.setValue(a4);
		return param.addObject(7, this);
	}

	return -1;
}



int
SmoothIMK::updateParameter(int parameterID, Information& info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->Fy = info.theDouble;
		break;
	case 2:
		this->E0 = info.theDouble;
		break;
	case 3:
		this->b = info.theDouble;
		break;
	case 4:
		this->a1 = info.theDouble;
		break;
	case 5:
		this->a2 = info.theDouble;
		break;
	case 6:
		this->a3 = info.theDouble;
		break;
	case 7:
		this->a4 = info.theDouble;
		break;
	default:
		return -1;
	}

	return 0;
}


void SmoothIMK::updateDamage()
{
	if (((statP == 1 || statP == 3) && sig < 0) || ((statP == 2 || statP == 4) && sig > 0))
	{
		//submit Pos damage and reset for new excursion
		double zeroSigEps = epsP - sigP / E0;
		double dE = 0.5 * sigP * (zeroSigEps - epsP);
		EnergyP += dE;
		if (EnergyP < 0) EnergyP = 0.;
		if (gama > 9999)
		{
			return;
		}
		double& Fyd = (statP == 2 || statP == 4) ? FydP : FydN;
		ExcurEnergy += dE;
		if (ExcurEnergy < 0) ExcurEnergy = 0.;
		double beta = pow(ExcurEnergy / (FailEnerg - EnergyP), c);
		if (beta > 0.999 || beta < 0)
		{
			opserr << "\nSmoothIMK:" << this->getTag() << " WARNING! Maximum Energy Absorbance Capacity Reached\n" << endln;
			beta = 0.999;
		}
		Fyd = (1. - beta) * Fyd + beta * resFac * Fyd;
		ExcurEnergy = 0.0;
	}
	else
	{
		double dE = 0.5 * (sig + sigP) * (eps - epsP);
		ExcurEnergy += dE;
		EnergyP += dE;
	}
}

