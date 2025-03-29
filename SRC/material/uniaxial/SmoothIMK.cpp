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
#include <MaterialResponse.h>

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

	if (numData < 20) {
		opserr << "Invalid SmoothIMK #args for: " << tag << " see the syntax" << endln;
		return 0;
	}
	//default parameters:
	numData = 1;
	double pf1, pd1, pf2, pd2, pf3, pd3, pdu;
	double nf1, nd1, nf2, nd2, nf3, nd3, ndu;
	double gama, c, r0, r1, r2, sigInit = 0;
	double pinchXPos = 0, pinchXNeg = 0, pinchYPos = 0, pinchYNeg = 0;
	int numPinchInput = 0;
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
	if (OPS_GetDoubleInput(&numData, &pdu) != 0) {
		opserr << "SmoothIMK:: invalid pdu for material : " << tag << endln;
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
	if (OPS_GetDoubleInput(&numData, &ndu) != 0) {
		opserr << "SmoothIMK:: invalid ndu for material : " << tag << endln;
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
	if (rule < 1 || rule > 3)
	{
		opserr << "SmoothIMK::invalid cyclicRule for material : " << tag
			<< "; can be 1,2 or 3 for bilinear, peak-oriented or piched rules, respectively" << endln;
		return 0;
	}
	if (rule == 3)
	{
		numPinchInput = 4;
		if (OPS_GetDoubleInput(&numData, &pinchXPos) != 0) {
			opserr << "SmoothIMK:: invalid pinchXPos for material : " << tag << endln;
			return 0;
		}
		if (OPS_GetDoubleInput(&numData, &pinchYPos) != 0) {
			opserr << "SmoothIMK:: invalid pinchYPos for material : " << tag << endln;
			return 0;
		}
		if (OPS_GetDoubleInput(&numData, &pinchXNeg) != 0) {
			opserr << "SmoothIMK:: invalid pinchXNeg for material : " << tag << endln;
			return 0;
		}
		if (OPS_GetDoubleInput(&numData, &pinchYNeg) != 0) {
			opserr << "SmoothIMK:: invalid pinchYNeg for material : " << tag << endln;
			return 0;
		}
	}
	if (OPS_GetNumRemainingInputArgs() > 20 + numPinchInput)
	{
		if (OPS_GetDoubleInput(&numData, &sigInit) != 0) {
			opserr << "SmoothIMK:: invalid sigInit for material : " << tag << endln;
			return 0;
		}
	}

	// Parsing was successful, allocate the material
	theMaterial = new SmoothIMK(tag, pd1, pf1, pd2, pf2, pd3, pf3, pdu, nd1, nf1, nd2, nf2, nd3, nf3, ndu, gama, c, r0, r1, r2, rule, pinchXPos, pinchYPos,
		pinchXNeg, pinchYNeg, sigInit);


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
	double _pdu,
	double _nd1, double _nf1,
	double _nd2, double _nf2,
	double _nd3, double _nf3,
	double _ndu,
	double _gama, double _c,
	double _r0, double _r1, double _r2,
	int _cyclicRule,
	double _pinchXPos, double _pinchYPos,
	double _pinchXNeg, double _pinchYNeg,
	double sigInit) :
	UniaxialMaterial(tag, MAT_TAG_SmoothIMK),
	pd1(_pd1), pf1(_pf1),
	pd2(_pd2), pf2(_pf2),
	pd3(_pd3), pf3(_pf3),
	nd1(_nd1), nf1(_nf1),
	nd2(_nd2), nf2(_nf2),
	nd3(_nd3), nf3(_nf3),
	pdu(_pdu), ndu(_ndu),
	c(_c), FailEnerg(_gama* pf1* pd1),
	r0(_r0), r1(_r1), r2(_r2),
	cyclicRule(_cyclicRule),
	pinchXPos(_pinchXPos), pinchYPos(_pinchYPos),
	pinchXNeg(_pinchXNeg), pinchYNeg(_pinchYNeg),
	sigini(sigInit)
{
	revertToStart();
}

SmoothIMK::SmoothIMK(void) :
	UniaxialMaterial(0, MAT_TAG_SmoothIMK),
	pd1(0), pf1(0),
	pd2(0), pf2(0),
	pd3(0), pf3(0),
	nd1(0), nf1(0),
	nd2(0), nf2(0),
	nd3(0), nf3(0),
	pdu(0), ndu(0),
	c(0), FailEnerg(0),
	r0(0), r1(0), r2(0),
	cyclicRule(1),
	sigini(0)
{
	revertToStart();
}

SmoothIMK::~SmoothIMK(void)
{
	// Does nothing
}

UniaxialMaterial*
SmoothIMK::getCopy(void)
{
	SmoothIMK* theCopy = new SmoothIMK(this->getTag(), pd1, pf1, pd2, pf2, pd3,
		pf3, pdu, nd1, nf1, nd2, nf2, nd3, nf3, ndu, FailEnerg / pf1 / pd1, c, r0, r1, r2, cyclicRule,
		pinchXPos, pinchYPos, pinchXNeg, pinchYNeg, sigini);
	theCopy->revertToStart();
	return theCopy;
}

double
SmoothIMK::getInitialTangent(void)
{
	return E0p;
}

double SmoothIMK::getStrain(void)
{
	return eps;
}

double SmoothIMK::getStress(void)
{
	return sig;
}

double SmoothIMK::getTangent(void)
{
	return e;
}

int
SmoothIMK::setParameter(const char** argv, int argc, Parameter& param)
{
	return -1;
}

int
SmoothIMK::updateParameter(int parameterID, Information& info)
{
	return -1;
}

int
SmoothIMK::commitState(void)
{
	epsminP = epsmin;
	epsmaxP = epsmax;
	epsLimitP = epsLimit;
	epss0P = epss0;
	sigs0P = sigs0;
	epssrP = epsr;
	sigsrP = sigr;
	slopeRatP = slopeRat;
	onEnvelopeP = onEnvelope;
	epsPlP = epsPl;
	updateDamage();
	isPosDirP = isPosDir;
	branchP = branch;
	eP = e;
	sigP = sig;
	epsP = eps;
	R0P = R0;
	return 0;
}

int
SmoothIMK::revertToLastCommit(void)
{
	epsmin = epsminP;
	epsmax = epsmaxP;
	epsLimit = epsLimitP;
	epss0 = epss0P;
	sigs0 = sigs0P;
	epsr = epssrP;
	sigr = sigsrP;
	branch = branchP;
	isPosDir = isPosDirP;
	slopeRat = slopeRatP;
	R0 = R0P;
	e = eP;
	sig = sigP;
	eps = epsP;
	onEnvelope = onEnvelopeP;
	epsPl = epsPlP;
	return 0;
}

int
SmoothIMK::revertToStart(void)
{
	EnergyP = 0;	//by SAJalali
	E0p = pf1 / (pd1 != 0 ? pd1 : 1);
	E0n = nf1 / (nd1 != 0 ? nd1 : 1);
	FydP = pf1;
	FydN = nf1;
	FcP = pf2;
	FcN = nf2;
	EshP = (FcP - FydP) / (pd2 - pd1);
	EshN = (FcN - FydN) / (nd2 - nd1);
	FrP = pf3;
	FrN = nf3;
	eP = E0p;
	epsP = sigini / E0p;
	sigP = sigini;
	sig = 0.0;
	eps = 0.0;
	e = E0p;

	isPosDirP = isPosDir = true;
	branchP = branch = precap;
	epsmaxP = pd1;
	epsminP = nd1;
	epsLimitP = pd1;
	epss0P = 0.0;
	sigs0P = 0.0;
	epssrP = 0.0;
	sigsrP = 0.0;
	R0P = r0;
	epsPl = epsPlP = 0;
	ExcurEnergy = 0;
	slopeRat = slopeRatP = 0;
	onEnvelope = onEnvelopeP = true;
	updateAsymptote();
	R0P = R0;
	epss0P = epss0;
	sigs0P = sigs0;
	slopeRatP = slopeRat;
	epsLimitP = epsLimit;
	return 0;
}

int
SmoothIMK::sendSelf(int commitTag, Channel& theChannel)
{
	static Vector data(37);
	int n = -1;
	data(n++) = this->getTag();	//0
	data(n++) = pd1;
	data(n++) = pd2;
	data(n++) = pd3;
	data(n++) = pdu;
	data(n++) = nd1;
	data(n++) = nd2;
	data(n++) = nd3;
	data(n++) = ndu;
	data(n++) = pf1;
	data(n++) = pf2;
	data(n++) = pf3;
	data(n++) = nf1;
	data(n++) = nf2;
	data(n++) = nf3;
	data(n++) = FailEnerg;
	data(n++) = c;
	data(n++) = sigini;
	data(n++) = epsP;
	data(n++) = sigP;
	data(n++) = eP;
	data(n++) = EnergyP;
	data(n++) = epsminP;
	data(n++) = epsmaxP;
	data(n++) = epsLimitP;
	data(n++) = epss0P;
	data(n++) = sigs0P;
	data(n++) = epssrP;
	data(n++) = sigsrP;
	data(n++) = isPosDirP;
	data(n++) = branchP;
	data(n++) = FydP;
	data(n++) = FydN;
	data(n++) = ExcurEnergy;
	data(n++) = slopeRatP;
	data(n++) = onEnvelopeP;
	data(n++) = R0P; // 36
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
	static Vector data(37);	//editted by SAJalali for energy

	if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "SmoothIMK::recvSelf() - failed to recvSelf\n";
		return -1;
	}
	int n = -1;
	this->setTag(data(n++));
	pd1 = data(n++);
	pd2 = data(n++);
	pd3 = data(n++);
	pdu = data(n++);
	nd1 = data(n++);
	nd2 = data(n++);
	nd3 = data(n++);
	ndu = data(n++);
	pf1 = data(n++);
	pf2 = data(n++);
	pf3 = data(n++);
	nf1 = data(n++);
	nf2 = data(n++);
	nf3 = data(n++);
	FailEnerg = data(n++);
	c = data(n++);
	sigini = data(n++);
	epsP = data(n++);
	sigP = data(n++);
	eP = data(n++);
	EnergyP = data(n++);
	epsminP = data(n++);
	epsmaxP = data(n++);
	epsLimitP = data(n++);
	epss0P = data(n++);
	sigs0P = data(n++);
	epssrP = data(n++);
	sigsrP = data(n++);
	isPosDirP = data(n++);
	branchP = (ebranch)data(n++);
	FydP = data(n++);
	FydN = data(n++);
	ExcurEnergy = data(n++);
	slopeRatP = data(n++);
	onEnvelopeP = data(n++);
	R0P = data(n++);
	e = eP;
	sig = sigP;
	eps = epsP;

	return 0;
}

void
SmoothIMK::Print(OPS_Stream& s, int flag)
{
	//    s << "SmoothIMK:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
	const char* endStr = endln;
	if (flag == OPS_PRINT_PRINTMODEL_JSON)
	{
		endStr = "";
		s << "\t\t\t{";
	}
	s << "\"name \": \"" << this->getTag() << "\", " << endStr;
	s << "\"type\": \"SmoothIMK\", " << endStr;
	s << "\"pd1 \": " << pd1 << ", " << endStr;
	s << "\"pd2 \": " << pd2 << ", " << endStr;
	s << "\"pd3 \": " << pd3 << ", " << endStr;
	s << "\"nd1 \": " << nd1 << ", " << endStr;
	s << "\"nd2 \": " << nd2 << ", " << endStr;
	s << "\"nd3 \": " << nd3 << ", " << endStr;
	s << "\"pf1 \": " << pf1 << ", " << endStr;
	s << "\"pf2 \": " << pf2 << ", " << endStr;
	s << "\"pf3 \": " << pf3 << ", " << endStr;
	s << "\"nf1 \": " << nf1 << ", " << endStr;
	s << "\"nf2 \": " << nf2 << ", " << endStr;
	s << "\"nf3 \": " << nf3 << ", " << endStr;
	s << "\"sigini \": " << sigini << endStr;
	if (flag == OPS_PRINT_PRINTMODEL_JSON)
		s << "}";

}


Response* SmoothIMK::setResponse(const char** argv, int argc, OPS_Stream& theOutput)
{
	Response* theResponse = UniaxialMaterial::setResponse(argv, argc, theOutput);
	if (theResponse != 0)
		return theResponse;
	if (strcmp(*argv, "branch") == 0)
	{
		theResponse = new MaterialResponse(this, 11, 0);
		return theResponse;
	}
	return 0;
}

int SmoothIMK::getResponse(int responseID, Information& matInformation)
{
	int res = 0;
	res = UniaxialMaterial::getResponse(responseID, matInformation);
	if (res == 0)
		return 0;
	switch (responseID)
	{
	case 11:
		matInformation.setInt(branchP);
		return 0;
	default:
		return -1;
	}
}


void SmoothIMK::updateDamage()
{
	if (sigP > 0 && sig < 0)
	{
		double zeroSigEps = epsP - sigP / (pf1 / pd1);
		double dE = 0.5 * sigP * (zeroSigEps - epsP);
		EnergyP += dE;
		if (EnergyP < 0) EnergyP = 0.;
		ExcurEnergy += dE;
		if (ExcurEnergy < 0) ExcurEnergy = 0.;
		if (branch == failing)
		{
			FydP = FrP;
			FcP = FrP;
			EshP = (FcP - FydP) / (pd2 - pd1);
		}
		else if (branch == failed)
		{
			FydP = 0;
			FcP = FrP;
			EshP = (FcP - FydP) / (pd2 - pd1);
		}
		else if (FailEnerg / pd1 / pf1 < 9999)
		{

			double beta = pow(ExcurEnergy / (FailEnerg - EnergyP), c);
			if (beta < 0)
				beta = 0;
			if (beta > 0.9999)
			{
				opserr << "\nSmoothIMK:" << this->getTag() << " WARNING! Maximum Energy Absorbance Capacity Reached\n" << endln;
				beta = 0.9999;
			}
			FydP = (1. - beta) * FydP + beta * pf3 / pf1 * FydP;
			FcP = pf2 - pf1 + FydP;
			EshP = (FcP - FydP) / (pd2 - pd1);
			//FrP = FcP - pf2 + pf3;
			//if (FrP < 0)
			//	FrP = 0;
		}
		ExcurEnergy = 0.0;
	}
	else if (sigP < 0 && sig > 0)
	{
		double zeroSigEps = epsP - sigP / (nf1 / nd1);
		double dE = 0.5 * sigP * (zeroSigEps - epsP);
		EnergyP += dE;
		if (EnergyP < 0) EnergyP = 0.;
		ExcurEnergy += dE;
		if (ExcurEnergy < 0) ExcurEnergy = 0.;
		if (branch == failing)
		{
			FydN = FrN;
			FcN = FrN;
			EshN = (FcN - FydN) / (nd2 - nd1);
		}
		else if (branch == failed)
		{
			FydN = 0;
			FcN = 0;
			EshN = (FcN - FydN) / (nd2 - nd1);
		}
		else if (FailEnerg / nd1 / nf1 < 9999)
		{
			double beta = pow(ExcurEnergy / (FailEnerg - EnergyP), c);
			if (beta > 0.999 || beta < 0)
			{
				opserr << "\nSmoothIMK:" << this->getTag() << " WARNING! Maximum Energy Absorbance Capacity Reached\n" << endln;
				beta = 0.999;
			}
			FydN = (1. - beta) * FydN + beta * nf3 / nf1 * FydN;
			FcN = nf2 - nf1 + FydN;
			EshN = (FcN - FydN) / (nd2 - nd1);
			//FrN = FcN - nf2 + nf3;
			//if (FrN > 0)
			//	FrN = 0;
		}
		ExcurEnergy = 0.0;
	}
	else
	{
		double dE = 0.5 * (sig + sigP) * (eps - epsP);
		ExcurEnergy += dE;
		EnergyP += dE;
	}
}

int SmoothIMK::setTrialStrain(double trialStrain, double strainRate)
{
	revertToLastCommit();
	eps = trialStrain + (sigini > 0 ? sigini / E0p : sigini / E0n);
	const double deps = eps - epsP;
	if (fabs(deps) < 1.e-20)
		return 0;
	if (isPosDirP && deps < 0.0)
	{
		isPosDir = false;
		epsPl = epsP - sigP / E0p;
		if (epsPl > 0.01 * pd1)
			onEnvelope = false;
		branch = cyclicRule == 1 ? precap : (cyclicRule == 2 ? peakOriented : pinching);
		if (epsP > epsmaxP)
			epsmax = epsP;
		updateAsymptote();
	}
	else if (deps > 0.0 && !isPosDirP) {
		isPosDir = true;
		epsPl = epsP - sigP / E0p;
		if (epsPl < 0.01 * nd1)
			onEnvelope = false;
		branch = cyclicRule == 1 ? precap : (cyclicRule == 2 ? peakOriented : pinching);
		if (epsP < epsminP)
			epsmin = epsP;
		updateAsymptote();
	}
	else if ((isPosDirP && eps > epsLimit) || (!isPosDirP && eps < epsLimit))
	{
		branch = nextBranch(branch);
		updateAsymptote();
	}
	double epsrat = (eps - epsr) / (epss0 - epsr);
	double dum1 = 1.0 + pow(fabs(epsrat), R0);
	double dum2 = pow(dum1, (1 / R0));

	sig = slopeRat * epsrat + (1.0 - slopeRat) * epsrat / dum2;
	sig = sig * (sigs0 - sigr) + sigr;

	e = slopeRat + (1.0 - slopeRat) / (dum1 * dum2);
	e *= (sigs0 - sigr) / (epss0 - epsr);
	return 0;
}

void SmoothIMK::updateAsymptote()
{
	switch (cyclicRule)
	{
	case 1:
		bilinAsymptote();
		return;
	case 2:
		peakOrientedAsymptote();
		return;
	case 3:
		pinchedAsymptote();
		return;
	default:
		opserr << "SmoothIMK::ERROR Unrecognized cyclic rule: " << cyclicRule << endln;
		exit(-1);

	}
}

void SmoothIMK::bilinAsymptote()
{
	//updates: epsr, sigr, epss0, sigs0, epsLimit, slopeRat, R0, stat
	const double& Esh = isPosDir ? EshP : EshN;
	const double& Fy = isPosDir ? FydP : FydN;
	const double& Fc = isPosDir ? FcP : FcN;
	const double& Fr = isPosDir ? FrP : FrN;
	const double& dy = isPosDir ? pd1 : nd1;
	const double& dc = isPosDir ? pd2 : nd2;
	const double& dr = isPosDir ? pd3 : nd3;
	const double& du = isPosDir ? pdu : ndu;
	const double& E1 = isPosDir ? E0p : E0n;
	const double& E2 = branchP == peakOriented ? eP : (isPosDir ? E0n : E0p);
	epsr = epsP;
	sigr = sigP;
	double Ec = (Fr - Fc) / (dr - dc);
	double k2 = 0;
	switch (branch)
	{
	case precap:
		epss0 = (Fy - Esh * dy - sigr * E1 / E2 + E1 * epsr) / (E1 - Esh);
		if ((!isPosDir && epss0 > epsr) || (isPosDir && epss0 < epsr))
		{
			//we should skip this branch
			branch = postcap;
			//epss0 = dc;
			//sigs0 = Fc;
			//epsLimit = (dc + dr) / 2;
			//k2 = Ec;
		}
		else
		{
			sigs0 = Fy + Esh * (epss0 - dy);
			epsLimit = (dy + dc) / 2;
			k2 = Esh;
			break;
		}
	case postcap:
		epss0 = dc;
		sigs0 = Fc;
		epsLimit = (dc + dr) / 2;
		k2 = Ec;
		break;
	case residual:
		epss0 = dr;
		sigs0 = Fr;
		epsLimit = (dr + du) / 2;
		k2 = 0;
		break;
	case failing:
		epss0 = du;
		sigs0 = Fr;
		epsLimit = (3 * du - dr) / 2;
		k2 = Ec;
		break;
	case failed:
		epss0 = (4 * du - 2 * dr) / 2;
		sigs0 = 0;
		epsLimit = 1000 * dy;
		k2 = 0;
	}
	double k1 = (sigs0 - sigr) / (epss0 - epsr);
	slopeRat = k2 / k1;
	R0 = r0;
	if (r1 != 0.0)
	{
		double xi_1 = (1 - slopeRat) * (sigs0 - sigr) / (epss0 - epsr) / E1;
		double xi_2 = (epss0 - epsr) / dy;
		R0 += r1 * xi_1 + r2 * xi_2;
		R0 = R0 < 100 ? R0 : 100;
		//if (print)
		//	opserr << epsP << " " << xi_1 << " " << xi_2 << " " << R0 << endln;
	}

}

void SmoothIMK::peakOrientedAsymptote()
{
	//updates: epsr, sigr, epss0, sigs0, epsLimit, slopeRat, R0, stat
	if (onEnvelope)
		return bilinAsymptote();
	//const double& E2 = isPosDir ? E0n : E0p;
	const double& epsPeak = isPosDir ? epsmax : epsmin;
	const double& dy = isPosDir ? pd1 : nd1;
	const double& E1 = isPosDir ? E0p : E0n;
	//peakOriented branch
	double tmp;
	double k2;
	epsr = epsP;
	sigr = sigP;
	if (branch == peakOriented)
	{
		epss0 = epsPl;
		sigs0 = 0;
		epsLimit = (epss0 + epsPeak) / 2;
		double sigmax;
		ebranch targBranch; // not used
		double tmp;
		getEnvelope(epsPeak, sigmax, targBranch, tmp, tmp);
		k2 = sigmax / (epsPeak - epss0);
	}
	else
	{
		epss0 = epsPeak;
		getEnvelope(epss0, sigs0, branch, k2, epsLimit);
	}
	double k1 = (sigs0 - sigr) / (epss0 - epsr);
	slopeRat = k2 / k1;
	R0 = r0;
	if (r1 != 0.0)
	{
		double xi_1 = (1 - slopeRat) * (sigs0 - sigr) / (epss0 - epsr) / E1;
		double xi_2 = (epss0 - epsr) / dy;
		R0 += r1 * xi_1 + r2 * xi_2;
		R0 = R0 < 100 ? R0 : 100;
		//if (print)
		//	opserr << epsP << " " << xi_1 << " " << xi_2 << " " << R0 << endln;
	}

}

void SmoothIMK::pinchedAsymptote()
{
	if (onEnvelope)
		return bilinAsymptote();
	const double& E2 = isPosDir ? E0n : E0p;
	const double& epsPeak = isPosDir ? epsmax : epsmin;
	const double& dy = isPosDir ? pd1 : nd1;
	const double& E1 = isPosDir ? E0p : E0n;
	const double pinchX = isPosDir ? pinchXPos : pinchXNeg;
	const double pinchY = isPosDir ? pinchYPos : pinchYNeg;
	//peakOriented branch
	double tmp;
	double k2;
	epsr = epsP;
	sigr = sigP;
	if (branch == pinching)
	{
		epss0 = epsPl;
		sigs0 = 0;
		double x = epsPl + pinchX * (epsPeak - epsPl);
		epsLimit = (epss0 + x) / 2;
		double sigmax;
		ebranch targBranch; // not used
		double tmp;
		getEnvelope(epsPeak, sigmax, targBranch, tmp, tmp);
		double y = pinchY * sigmax;
		k2 = y / (x - epss0);
	}
	else if (branch == peakOriented)
	{
		epss0 = epsPl + pinchX * (epsPeak - epsPl);
		epsLimit = (epss0 + epsPeak) / 2;
		double sigmax;
		ebranch targBranch; // not used
		double tmp;
		getEnvelope(epsPeak, sigmax, targBranch, tmp, tmp);
		sigs0 = pinchY * sigmax;
		k2 = (sigmax - sigs0) / (epsPeak - epss0);
	}
	else
	{
		epss0 = epsPeak;
		getEnvelope(epss0, sigs0, branch, k2, epsLimit);
	}
	double k1 = (sigs0 - sigr) / (epss0 - epsr);
	slopeRat = k2 / k1;
	R0 = r0;
	if (r1 != 0.0)
	{
		double xi_1 = fabs(k1 - k2) / E1;
		double xi_2 = (epss0 - epsr) / dy;
		R0 += r1 * xi_1 + r2 * xi_2;
		R0 = R0 < 100 ? R0 : 100;
		//if (print)
		//	opserr << epsP << " " << xi_1 << " " << xi_2 << " " << R0 << endln;
	}
}

SmoothIMK::ebranch SmoothIMK::nextBranch(ebranch branch)
{
	onEnvelope = true;
	const double& epsPeak = isPosDir ? epsmax : epsmin;
	switch (branch)
	{
	case precap:
		return postcap;
	case postcap:
		return residual;
	case residual:
		if (isPosDir)
		{
			if (FrP == 0)
				return failed;
			else
				return failed;
		}
		else
		{
			if (FrN == 0)
				return failed;
			else
				return failed;
		}
	case failing:
		return failed;
	case failed:
		return failed;
	case peakOriented:
		onEnvelope = false;
		double tmp;
		ebranch peakBranch;
		getEnvelope(epsPeak, tmp, peakBranch, tmp, tmp);
		return peakBranch;
	case pinching:
		onEnvelope = false;
		return peakOriented;
	}
}

void SmoothIMK::getEnvelope(double eps, double& targStress, SmoothIMK::ebranch& targBranch, double& k, double& limitEps)
{
	const double& Esh = eps > 0 ? EshP : EshN;
	const double& Fy = eps > 0 ? FydP : -FydN;
	const double& Fc = eps > 0 ? FcP : -FcN;
	const double& Fr = eps > 0 ? FrP : -FrN;
	const double& dy = eps > 0 ? pd1 : -nd1;
	const double& dc = eps > 0 ? pd2 : -nd2;
	const double& dr = eps > 0 ? pd3 : -nd3;
	const double& du = eps > 0 ? pdu : -ndu;
	const double& E = eps > 0 ? E0p : E0n;
	int sgn = eps > 0 ? 1 : -1;
	eps = fabs(eps);
	if (eps < dy)
	{
		targStress = sgn * E * eps;
		targBranch = precap;
		k = E;
		limitEps = sgn * (dy + dc) / 2.;
		return;
	}
	else if (eps < dc)
	{
		targStress = sgn * (Fy + Esh * (eps - dy));
		k = Esh;
		targBranch = precap;
		limitEps = sgn * (eps + dc) / 2;
		return;
	}
	else if (eps < dr)
	{
		k = (Fr - Fc) / (dr - dc);
		targStress = sgn * (Fc + k * (eps - dc));
		targBranch = postcap;
		limitEps = sgn * (eps + dr) / 2;
		return;
	}
	else if (eps < du)
	{
		targStress = sgn * Fr;
		k = 0;
		targBranch = residual;
		limitEps = sgn * (eps + du) / 2;
		return;
	}
	else // >= du
	{
		targStress = 0;
		targBranch = failed;
		k = 0;
		limitEps = sgn * 1000 * dy;
	}

}
