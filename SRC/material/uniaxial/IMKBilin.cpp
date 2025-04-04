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

#include <math.h>
#include <IMKBilin.h>
#include <elementAPI.h>
#include <Vector.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <algorithm>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
using namespace std;

static int numIMKBilinMaterials = 0;

void*
OPS_IMKBilin(void)
{
	if (numIMKBilinMaterials == 0) {
		numIMKBilinMaterials++;
		opserr << "IMK with Bilinear Response - Code by AE_KI (Nov22)\n";
	}

	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;

	int    iData[1];
	double dData[21];
	int numData = 1;

	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial IMKBilin tag" << endln;
		return 0;
	}

	numData = 21;

	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid Args want: uniaxialMaterial IMKBilin tag? Ke? ";
		opserr << "dp_pos? dpc_pos? du_pos? Fy_pos? FmaxFy_pos? FresFy_pos? ";
		opserr << "dp_neg? dpc_neg? du_neg? Fy_neg? FmaxFy_neg? FresFy_neg? ";
		opserr << "LamdaS?  LamdaC? LamdaK? Cs? Cc? Ck? D_pos? D_neg? ";
		return 0;
	}


	// Parsing was successful, allocate the material
	theMaterial = new IMKBilin(iData[0],
		dData[0],
		dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
		dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
		dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20]);

	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type IMKBilin Material\n";
		return 0;
	}

	return theMaterial;
}

IMKBilin::IMKBilin(int tag, double p_Ke,
	double p_posUp_0, double p_posUpc_0, double p_posUu_0, double p_posFy_0, double p_posFcapFy_0, double p_posFresFy_0,
	double p_negUp_0, double p_negUpc_0, double p_negUu_0, double p_negFy_0, double p_negFcapFy_0, double p_negFresFy_0,
	double p_LAMBDA_S, double p_LAMBDA_C, double p_LAMBDA_K, double p_c_S, double p_c_C, double p_c_K, double p_D_pos, double p_D_neg)
	:UniaxialMaterial(tag, MAT_TAG_IMKBilin), Ke(p_Ke),
	posUp_0(p_posUp_0), posUpc_0(p_posUpc_0), posUu_0(p_posUu_0), posFy_0(p_posFy_0), posFcapFy_0(p_posFcapFy_0), posFresFy_0(p_posFresFy_0),
	negUp_0(p_negUp_0), negUpc_0(p_negUpc_0), negUu_0(p_negUu_0), negFy_0(p_negFy_0), negFcapFy_0(p_negFcapFy_0), negFresFy_0(p_negFresFy_0),
	LAMBDA_S(p_LAMBDA_S), LAMBDA_C(p_LAMBDA_C), LAMBDA_K(p_LAMBDA_K), c_S(p_c_S), c_C(p_c_C), c_K(p_c_K), D_pos(p_D_pos), D_neg(p_D_neg)
{
	// Make sure these are all positive 
	if (negUp_0 < 0)
		negUp_0 = -negUp_0;
	if (negUpc_0 < 0)
		negUpc_0 = -negUpc_0;
	if (negUu_0 < 0)
		negUu_0 = -negUu_0;
	if (negFy_0 < 0)
		negFy_0 = -negFy_0;

	this->revertToStart();
}

IMKBilin::IMKBilin()
	:UniaxialMaterial(0, MAT_TAG_IMKBilin), Ke(0),
	posUp_0(0), posUpc_0(0), posUu_0(0), posFy_0(0), posFcapFy_0(0), posFresFy_0(0),
	negUp_0(0), negUpc_0(0), negUu_0(0), negFy_0(0), negFcapFy_0(0), negFresFy_0(0),
	LAMBDA_S(0), LAMBDA_C(0), LAMBDA_K(0), c_S(0), c_C(0), c_K(0), D_pos(0), D_neg(0)
{
	this->revertToStart();
}

IMKBilin::~IMKBilin()
{
	// does nothing
}

int IMKBilin::setTrialStrain(double strain, double strainRate)
{
	//all variables to the last commit
	this->revertToLastCommit();

	//state determination algorithm: defines the current force and tangent stiffness
	double Ui_1 = Ui;
	double Fi_1 = Fi;
	U = strain; //set trial displacement
	Ui = U;
	double dU = Ui - Ui_1;    // Incremental deformation at current step
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////  MAIN CODE //////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (Failure_State == 5) {     // When a failure has already occured
		Fi = 0;
	}
	else if (dU == 0) {   // When deformation doesn't change from the last
		Fi = Fi_1;
	}
	else {
		double  betaS = 0, betaC = 0, betaK = 0;
		bool    FailS = false, FailC = false, FailK = false;
		///////////////////////////////////////////////////////////////////////////////////////////
		/////////////////// WHEN REVERSAL /////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		if ((onBackbone && Fi_1 * dU < 0) || (onBackbone && Fi_1 == 0 && Ui_1 * dU <= 0)) {
			onBackbone = false;
			/////////////////////////// UPDATE PEAK POINTS ////////////////////////////////////////////
			/////////////////// UPDATE UNLOADING STIFFNESS ////////////////////////////////////////////
			double  EpjK = engAcml - 0.5 * (Fi_1 / Kunload) * Fi_1;
			double  EiK = engAcml - engDspt - 0.5 * (Fi_1 / Kunload) * Fi_1;
			betaK = pow((EiK / (engRefK - EpjK)), c_K);
			FailK = (betaK > 1);
			betaK = betaK < 0 ? 0 : (betaK > 1 ? 1 : betaK);
			Kunload *= (1 - betaK);
			if (Failure_State > 1) {
				Kunload = 0.5 * Ke;
			}
			KgetTangent = Kunload;
		}
		Fi = Fi_1 + Kunload * dU;
		///////////////////////////////////////////////////////////////////////////////////////////
		/////////////////// WHEN NEW EXCURSION /////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		if (!onBackbone && Fi_1 * Fi <= 0.0 && Failure_State > 0) {
			/////////////////// UPDATE BACKBONE CURVE /////////////////////////////////////////////////
			double Ei = max(0.0, engAcml - engDspt);
			betaS = pow((Ei / (engRefS - engAcml)), c_S);
			betaC = pow((Ei / (engRefC - engAcml)), c_C);
			FailS = (betaS > 1);
			FailC = (betaC > 1);
			betaS = betaS < 0 ? 0 : (betaS > 1 ? 1 : betaS);
			betaC = betaC < 0 ? 0 : (betaC > 1 ? 1 : betaC);
			engDspt = engAcml;
			// Positive
			if (dU > 0) {
				double FcapProj = posFcap - posKpc * posUcap;
				// Yield Point
				posFy *= (1 - betaS * D_pos);
				posKp *= (1 - betaS * D_pos); // Post-Yield Stiffness
				FcapProj *= (1 - betaC * D_pos);
				posUy = posFy / Ke;
				posKpc = posFy < posFres ? 0 : -posKpc_0 * (posFy - posFres) / (posFy_0 - posFres);
				// Capping Point
				double FyProj = posFy - posKp * posUy;
				posUcap = posKp <= posKpc ? 0 : (FcapProj - FyProj) / (posKp - posKpc);
				posFcap = FyProj + posKp * posUcap;
				// When a part of backbone is beneath the residual strength
				double candidateKp = (posFcap - posFres) / (posUcap - Ui_1 - (posFres - Fi_1) / Kunload);
				if (posFcap < posFres) {
					posFy = posFres;
					posFcap = posFres;
					posKp = 0;
					posKpc = 0;
					posUy = posFy / Ke;
					posUcap = 0;
				}
				else if (candidateKp > 0 && posKp > candidateKp) {
					posKp = candidateKp;
				}
				// Global Peak on the Updated Backbone
			}
			// Negative
			else {
				double FcapProj = negFcap - negKpc * negUcap;
				// Yield Point
				negFy *= (1 - betaS * D_neg);
				negKp *= (1 - betaS * D_neg); // Post-Yield Stiffness
				FcapProj *= (1 - betaC * D_neg);
				negUy = negFy / Ke;
				negKpc = negFy > negFres ? 0 : -negKpc_0 * (negFy - negFres) / (-negFy_0 - negFres);
				// Capping Point
				double FyProj = negFy - negKp * negUy;
				negUcap = negKp <= negKpc ? 0 : (FcapProj - FyProj) / (negKp - negKpc);
				negFcap = FyProj + negKp * negUcap;
				// When a part of backbone is beneath the residual strength
				double candidateKp = (negFcap - negFres) / (negUcap - Ui_1 - (negFres - Fi_1) / Kunload);
				if (negFcap > negFres) {
					negFy = negFres;
					negFcap = negFres;
					negKp = 0;
					negKpc = 0;
					negUy = negFy / Ke;
					negUcap = 0;
				}
				else if (candidateKp > 0 && negKp > candidateKp) {
					negKp = candidateKp;
				}
				// Global Peak on the Updated Backbone
			}
			////////////////////////// RELOADING TARGET DETERMINATION /////////////////////////////////
		}
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////// COMPUTE FORCE ASSUMING IT'S ON BACKBONE ///////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		// Positive
		if (dU > 0) {   // Backbone in Positive
			double Kbackbone = (Ui < posUcap) ? posKp : posKpc;
			double Fi_backbone = posFcap + (Ui - posUcap) * Kbackbone;
			if (Fi_backbone < posFres || Failure_State == 4) {
				Fi_backbone = posFres;
			}
			if (Failure_State == 3) {
				Fi_backbone = 0;
			}
			if (Fi > Fi_backbone) {
				Fi = Fi_backbone;
				onBackbone = true;
			}
		}
		// Negative
		else {          // Backbone in Negative
			double Kbackbone = (negUcap < Ui) ? negKp : negKpc;
			double Fi_backbone = negFcap + (Ui - negUcap) * Kbackbone;
			if (Fi_backbone > negFres || Failure_State == 3) {
				Fi_backbone = negFres;
			}
			if (Failure_State == 4) {
				Fi_backbone = 0;
			}
			if (Fi < Fi_backbone) {
				Fi = Fi_backbone;
				onBackbone = true;
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		// CHECK FOR FAILURE
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		// Failure_State
		//     0: Elastic
		//     1: Experienced Yield in either Positive or Negative.
		//     2: Experienced Residual Branch in either Positive or Negative.
		//     3: Strength Lost in Positive. Negative Backbone Strength will be set the Residual Strength.
		//     4: Strength Lost in Negative. Positive Backbone Strength will be set the Residual Strength.
		//     5: Strength Lost in Both Positive and Negative.
		bool    ResP, ResN;
		if (posKpc == 0) {    // Kpc can be zero when deteriorated Fy is smaller than Fres
			ResP = (Fi == posFres);
		}
		else {
			double  posUres = (posFres - posFcap + posKpc * posUcap) / posKpc;
			ResP = (Ui >= posUres);
		}
		if (negKpc == 0) {
			ResN = (Fi == negFres);
		}
		else {
			double  negUres = (negFres - negFcap + negKpc * negUcap) / negKpc;
			ResN = (Ui <= negUres);
		}

		bool	FailDp = (dU > 0 && Ui >= posUu_0);
		bool	FailDn = (dU < 0 && Ui <= -negUu_0);
		bool	FailRp = (dU > 0 && onBackbone && Fi <= 0);
		bool	FailRn = (dU < 0 && onBackbone && Fi >= 0);
		// int exFailure_State = Failure_State;
		if (FailS || FailC || FailK) {
			Fi = 0;
			Failure_State = 5;
		}
		else if (FailDp || FailRp) {
			Fi = 0;
			if (Failure_State == 4) {
				Failure_State = 5;
			}
			else {
				Failure_State = 3;
			}
		}
		else if (FailDn || FailRn) {
			Fi = 0;
			if (Failure_State == 3) {
				Failure_State = 5;
			}
			else {
				Failure_State = 4;
			}
		}
		else if (Failure_State < 2 && (ResP || ResN)) {
			Failure_State = 2;
		}
		else if (Failure_State == 0 && onBackbone) {
			Failure_State = 1;
		}
		// Failure_State Change check
				// if (Failure_State!=exFailure_State) {
				//     std::cout << exFailure_State << " -> " << Failure_State << "\n";
				// }

		engAcml += 0.5 * (Fi + Fi_1) * dU;   // Internal energy increment

		KgetTangent = (Fi - Fi_1) / dU;
	}
	if (KgetTangent == 0) {
		KgetTangent = 1e-6;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////// END OF MAIN CODE ///////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	return 0;
}

double IMKBilin::getStress(void)
{
	//cout << " getStress" << endln;
	return (Fi);
}

double IMKBilin::getTangent(void)
{
	//cout << " getTangent" << endln;
	return (KgetTangent);
}

double IMKBilin::getInitialTangent(void)
{
	//cout << " getInitialTangent" << endln;
	return (Ke);
}

double IMKBilin::getStrain(void)
{
	//cout << " getStrain" << endln;
	return (U);
}

int IMKBilin::commitState(void)
{
	//cout << " commitState" << endln;

	//commit trial  variables
// 7 Pos U and F
	cPosUy = posUy;
	cPosFy = posFy;
	cPosUcap = posUcap;
	cPosFcap = posFcap;
	cPosFres = posFres;
	cPosKp = posKp;
	cPosKpc = posKpc;
	// 7 Neg U and F
	cNegUy = negUy;
	cNegFy = negFy;
	cNegUcap = negUcap;
	cNegFcap = negFcap;
	cNegFres = negFres;
	cNegKp = negKp;
	cNegKpc = negKpc;
	// 3 State
	cU = U;
	cUi = Ui;
	cFi = Fi;
	// 1 Stiffness
	cKunload = Kunload;
	// 2 Energy
	cEngAcml = engAcml;
	cEngDspt = engDspt;
	// 2 Flag
	cFailure_State = Failure_State;
	cOnBackbone = onBackbone;
	return 0;
}

int IMKBilin::revertToLastCommit(void)
{
	//cout << " revertToLastCommit" << endln;
	//the opposite of commit trial history variables
// 7 Positive U and F
	posUy = cPosUy;
	posFy = cPosFy;
	posUcap = cPosUcap;
	posFcap = cPosFcap;
	posFres = cPosFres;
	posKp = cPosKp;
	posKpc = cPosKpc;
	// 7 Negative U and F
	negUy = cNegUy;
	negFy = cNegFy;
	negUcap = cNegUcap;
	negFcap = cNegFcap;
	negFres = cNegFres;
	negKp = cNegKp;
	negKpc = cNegKpc;
	// 3 State Variables
	U = cU;
	Ui = cUi;
	Fi = cFi;
	// 1 Stiffness
	Kunload = cKunload;
	// 2 Energy
	engAcml = cEngAcml;
	engDspt = cEngDspt;
	// 2 Flag
	Failure_State = cFailure_State;
	onBackbone = cOnBackbone;
	return 0;
}

int IMKBilin::revertToStart(void)
{
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\\
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ONE TIME CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\\
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	// 14 Initial Values
	posUy_0 = posFy_0 / Ke;
	posUcap_0 = posUy_0 + posUp_0;
	posFcap_0 = posFcapFy_0 * posFy_0;
	posKp_0 = (posFcap_0 - posFy_0) / posUp_0;
	posKpc_0 = posFcap_0 / posUpc_0;
	negUy_0 = negFy_0 / Ke;
	negUcap_0 = negUy_0 + negUp_0;
	negFcap_0 = negFcapFy_0 * negFy_0;
	negKp_0 = (negFcap_0 - negFy_0) / negUp_0;
	negKpc_0 = negFcap_0 / negUpc_0;
	engRefS = LAMBDA_S * posFy_0;
	engRefC = LAMBDA_C * posFy_0;
	engRefK = LAMBDA_K * posFy_0;
	// 7 Positive U and F
	posUy = cPosUy = posUy_0;
	posFy = cPosFy = posFy_0;
	posUcap = cPosUcap = posUcap_0;
	posFcap = cPosFcap = posFcap_0;
	posFres = cPosFres = posFy_0 * posFresFy_0;
	posKp = cPosKp = posKp_0;
	posKpc = cPosKpc = -posKpc_0;

	// 7 Negative U and F
	negUy = cNegUy = -negUy_0;
	negFy = cNegFy = -negFy_0;
	negUcap = cNegUcap = -negUcap_0;
	negFcap = cNegFcap = -negFcap_0;
	negFres = cNegFres = -negFy_0 * negFresFy_0;
	negKp = cNegKp = negKp_0;
	negKpc = cNegKpc = -negKpc_0;
	// 3 State Values
	U = cU = 0;
	Ui = cUi = 0;
	Fi = cFi = 0;
	// 2 Stiffness
	Kunload = cKunload = Ke;
	KgetTangent = Ke;
	// 2 Energy
	engAcml = cEngAcml = 0.0;
	engDspt = cEngDspt = 0.0;
	// 2 Flag
	Failure_State = cFailure_State = 0;
	onBackbone = cOnBackbone = false;
	return 0;
}

UniaxialMaterial*
IMKBilin::getCopy(void)
{
	IMKBilin* theCopy = new IMKBilin(this->getTag(), Ke,
		posUp_0, posUpc_0, posUu_0, posFy_0, posFcapFy_0, posFresFy_0,
		negUp_0, negUpc_0, negUu_0, negFy_0, negFcapFy_0, negFresFy_0,
		LAMBDA_S, LAMBDA_C, LAMBDA_K, c_S, c_C, c_K, D_pos, D_neg);

	*theCopy = *this;

	return theCopy;
}

int IMKBilin::sendSelf(int cTag, Channel& theChannel)
{
	int res = 0;

	static Vector data(57);
	data(0) = this->getTag();
	// 21 Fixed Input Material Parameters 1-25
	data(1) = Ke;
	data(2) = posUp_0;
	data(3) = posUpc_0;
	data(4) = posUu_0;
	data(5) = posFy_0;
	data(6) = posFcapFy_0;
	data(7) = posFresFy_0;
	data(8) = negUp_0;
	data(9) = negUpc_0;
	data(10) = negUu_0;
	data(11) = negFy_0;
	data(12) = negFcapFy_0;
	data(13) = negFresFy_0;
	data(14) = LAMBDA_S;
	data(15) = LAMBDA_C;
	data(16) = LAMBDA_K;
	data(17) = c_S;
	data(18) = c_C;
	data(19) = c_K;
	data(20) = D_pos;
	data(21) = D_neg;
	// 13 Initial Values 31-44
	data(22) = posUy_0;
	data(23) = posUcap_0;
	data(24) = posFcap_0;
	data(25) = posKp_0;
	data(26) = posKpc_0;
	data(27) = negUy_0;
	data(28) = negUcap_0;
	data(29) = negFcap_0;
	data(30) = negKp_0;
	data(31) = negKpc_0;
	data(32) = engRefS;
	data(33) = engRefC;
	data(34) = engRefK;
	// 7 Positive U and F 51-62
			/*
			data(51) 	= posUy;
			data(52) 	= posFy;
			data(53) 	= posUcap;
			data(54) 	= posFcap;
			data(60) 	= posFres;
			data(61) 	= posKp;
			data(62) 	= posKpc;
	// 3 State Variables 63-65
			data(63)    = U;
			data(64) 	= Ui;
			data(65) 	= Fi;
	// 1 Stiffness 66 67
			data(67) 	= Kunload;
	// 2 Energy 68 69
			data(68) 	= engAcml;
			data(69) 	= engDspt;
	// 7 Negative U and F 71-82
			data(71) 	= negUy;
			data(72) 	= negFy;
			data(73) 	= negUcap;
			data(74) 	= negFcap;
			data(80) 	= negFres;
			data(81) 	= negKp;
			data(82) 	= negKpc;
	// 2 Flag 85 86
			data(85)	= Failure_State;
			data(86) 	= onBackbone;
			*/
			// 7 Positive U and F 101-112
	data(35) = cPosUy;
	data(36) = cPosFy;
	data(37) = cPosUcap;
	data(38) = cPosFcap;
	data(39) = cPosFres;
	data(40) = cPosKp;
	data(41) = cPosKpc;
	// 3 State Variables 113-115
	data(42) = cU;
	data(43) = cUi;
	data(44) = cFi;
	// 2 Stiffness 116 117
	data(45) = cKunload;
	// 2 Energy 118 119
	data(46) = cEngAcml;
	data(47) = cEngDspt;
	// 7 Negative U and F 121-132
	data(48) = cNegUy;
	data(49) = cNegFy;
	data(50) = cNegUcap;
	data(51) = cNegFcap;
	data(52) = cNegFres;
	data(53) = cNegKp;
	data(54) = cNegKpc;
	// 2 Flag 135 136
	data(55) = cFailure_State;
	data(56) = cOnBackbone ? 1.0 : -1.0;
	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "IMKBilin::sendSelf() - failed to send data\n";

	return res;
}

int IMKBilin::recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res = 0;
	static Vector data(57);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);

	if (res < 0) {
		opserr << "IMKBilin::recvSelf() - failed to receive data\n";
		this->setTag(0);
	}
	else {
		this->setTag((int)data(0));
		// 21 Fixed Input Material Parameters
		Ke = data(1);
		posUp_0 = data(2);
		posUpc_0 = data(3);
		posUu_0 = data(4);
		posFy_0 = data(5);
		posFcapFy_0 = data(6);
		posFresFy_0 = data(7);
		negUp_0 = data(8);
		negUpc_0 = data(9);
		negUu_0 = data(10);
		negFy_0 = data(11);
		negFcapFy_0 = data(12);
		negFresFy_0 = data(13);
		LAMBDA_S = data(14);
		LAMBDA_C = data(15);
		LAMBDA_K = data(16);
		c_S = data(17);
		c_C = data(18);
		c_K = data(19);
		D_pos = data(20);
		D_neg = data(21);
		// 13 Initial Values
		posUy_0 = data(22);
		posUcap_0 = data(23);
		posFcap_0 = data(24);
		posKp_0 = data(25);
		posKpc_0 = data(26);
		negUy_0 = data(27);
		negUcap_0 = data(28);
		negFcap_0 = data(29);
		negKp_0 = data(30);
		negKpc_0 = data(31);
		engRefS = data(32);
		engRefC = data(33);
		engRefK = data(34);
		/*
			// 3 State Variables
					U               = data(63);
					Ui				= data(64);
					Fi				= data(65);
			// 1 Stiffness
					Kunload	    	= data(67);
			// 7 Positive U and F
					posUy			= data(51);
					posFy			= data(52);
					posUcap		    = data(53);
					posFcap		   	= data(54);
					posFres		   	= data(60);
					posKp			= data(61);
					posKpc		    = data(62);
			// 7 Negative U and F
					negUy			= data(71);
					negFy			= data(72);
					negUcap		    = data(73);
					negFcap		   	= data(74);
					negFres		    = data(80);
					negKp			= data(81);
					negKpc		    = data(82);
			// 2 Flag
					Failure_State	= data(85);
					onBackbone		= data(86);
			// 2 Energy
					engAcml			= data(68);
					engDspt		    = data(69);
		*/
		// 7 Positive U and F
		cPosUy = data(35);
		cPosFy = data(36);
		cPosUcap = data(37);
		cPosFcap = data(38);
		cPosFres = data(39);
		cPosKp = data(40);
		cPosKpc = data(41);
		// 3 State Variables
		cU = data(42);
		cUi = data(43);
		cFi = data(44);
		// 1 Stiffness
		cKunload = data(45);
		// 2 Energy
		cEngAcml = data(46);
		cEngDspt = data(47);
		// 7 Negative U and F
		cNegUy = data(48);
		cNegFy = data(49);
		cNegUcap = data(50);
		cNegFcap = data(51);
		cNegFres = data(52);
		cNegKp = data(53);
		cNegKpc = data(54);
		// 2 Flag
		cFailure_State = int(data(55));
		cOnBackbone = data(56) > 0.0 ? true : false;

		this->revertToLastCommit();
	}

	return res;
}

void IMKBilin::Print(OPS_Stream& s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"IMKPeakOriented\", ";
		s << "\"Ke\": " << Ke << ", ";
		s << "\"dp_pos\": " << posUp_0 << ", ";
		s << "\"dpc_pos\": " << posUpc_0 << ", ";
		s << "\"du_pos\": " << posUu_0 << ", ";
		s << "\"Fy_pos\": " << posFy_0 << ", ";
		s << "\"FmaxFy_pos\": " << posFcapFy_0 << ", ";
		s << "\"FresFy_pos\": " << posFresFy_0 << ", ";
		s << "\"dp_neg\": " << negUp_0 << ", ";
		s << "\"dpc_neg\": " << negUpc_0 << ", ";
		s << "\"du_neg\": " << negUu_0 << ", ";
		s << "\"Fy_neg\": " << negFy_0 << ", ";
		s << "\"FmaxFy_neg\": " << negFcapFy_0 << ", ";
		s << "\"FresFy_neg\": " << negFresFy_0 << ", ";
		s << "\"Lamda_S\": " << LAMBDA_S << ", ";
		s << "\"Lamda_C\": " << LAMBDA_C << ", ";
		s << "\"Lamda_K\": " << LAMBDA_K << ", ";
		s << "\"c_S\": " << c_S << ", ";
		s << "\"c_C\": " << c_C << ", ";
		s << "\"c_K\": " << c_K << ", ";
		s << "\"D_pos\": " << D_pos << ", ";
		s << "\"D_neg\": " << D_neg << "}";
		return;
	}
	cout << "IMKBilin tag: " << this->getTag() << endln;
}
