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
                                                                        
// $Revision: 1.3 $
// $Date: 2007-07-27 19:01:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Pinching4Material.cpp,v $
                                                                        
                                                                        
// Written: NM (nmitra@u.washington.edu) 
// Created: January 2002
// Updates: September 2004
//
// Description: This file contains the class implementation for 
// Pinching material which is defined by 4 points on the positive and 
// negative envelopes and a bunch of damage parameters. The material accounts for
// 3 types of damage rules : Strength degradation, Stiffness degradation, 
// unloading stiffness degradation. 
// Updates: damage calculations and several bug fixes


#include <Pinching4Material.h>
#include <OPS_Globals.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <OPS_Stream.h>

#include <string.h>
#include <sstream>
#include <elementAPI.h>

#include <Information.h>
#include <Parameter.h>

#include <Channel.h>


void* OPS_Pinching4Material()
{
    int numData = OPS_GetNumRemainingInputArgs();
    if (numData!=40 && numData!=29) {
	opserr << "WARNING: Insufficient arguments\n";
	return 0;
    }

    int tag;
    numData = 1;
    if (OPS_GetIntInput(&numData,&tag) < 0) {
	return 0;
    }

    UniaxialMaterial* mat = 0;
    int tDmg = -1;
    if (OPS_GetNumRemainingInputArgs() == 39) {
	double data[38];
	numData = 38;
	if (OPS_GetDoubleInput(&numData,data)) {
	    return 0;
	}
	const char* type = OPS_GetString();
	if (strcmp(type,"cycle")==0||strcmp(type,"Cycle")==0||strcmp(type,"DamageCycle")==0||strcmp(type,"damageCycle")==0) {
	    tDmg = 1;
	} else if (strcmp(type,"energy")==0||strcmp(type,"Energy")==0||strcmp(type,"DamageEnergy")==0||strcmp(type,"damageEnergy")==0) {
	    tDmg = 0;
	} else {
	    opserr << "WARNING invalid type of damage calculation specified\n";
	    opserr << "Pinching4 material: " << tag << endln;
	    return 0;
	}
	
	mat = new Pinching4Material(tag,data[0],data[1],data[2],data[3],data[4],data[5],
				    data[6],data[7],data[8],data[9],data[10],data[11],
				    data[12],data[13],data[14],data[15],
				    data[16],data[17],data[18],data[19],data[20],
				    data[21],data[22],data[23],
				    data[24],data[25],data[26],data[27],data[28],data[29],
				    data[30],data[31],
				    data[32],data[33],data[34],data[35],data[36],data[37],tDmg);
    } else if (OPS_GetNumRemainingInputArgs() == 28) {
	double data[27];
	numData = 27;
	if (OPS_GetDoubleInput(&numData,data)) {
	    return 0;
	}
	const char* type = OPS_GetString();
	if (strcmp(type,"cycle")==0||strcmp(type,"Cycle")==0||strcmp(type,"DamageCycle")==0||strcmp(type,"damageCycle")==0) {
	    tDmg = 1;
	} else if (strcmp(type,"energy")==0||strcmp(type,"Energy")==0||strcmp(type,"DamageEnergy")==0||strcmp(type,"damageEnergy")==0) {
	    tDmg = 0;
	} else {
	    opserr << "WARNING invalid type of damage calculation specified\n";
	    opserr << "Pinching4 material: " << tag << endln;
	    return 0;
	}
	
	mat = new Pinching4Material(tag,data[0],data[1],data[2],data[3],data[4],
				    data[5],data[6],data[7],
				    data[8],data[9],data[10],data[11],data[12],
				    data[13],data[14],data[15],
				    data[16],data[17],data[18],data[19],data[20],
				    data[21],data[22],data[23],
				    data[24],data[25],data[26],tDmg);
    }
    
    if (mat == 0) {
	opserr << "WARNING: failed to create Pinching4material material\n";
	return 0;
    }

    return mat;
}

Pinching4Material::Pinching4Material(int tag,
				     double f1p, double d1p, double f2p, double d2p,
				     double f3p, double d3p, double f4p, double d4p,
				     double f1n, double d1n, double f2n, double d2n,
				     double f3n, double d3n, double f4n, double d4n,
				     double mdp, double mfp, double msp,
				     double mdn, double mfn, double msn,
				     double gk1, double gk2, double gk3, double gk4, double gklim,
				     double gd1, double gd2, double gd3, double gd4, double gdlim,
				     double gf1, double gf2, double gf3, double gf4, double gflim, double ge, int dc):
  UniaxialMaterial(tag, MAT_TAG_Pinching4),
  stress1p(f1p), strain1p(d1p), stress2p(f2p), strain2p(d2p),
  stress3p(f3p), strain3p(d3p), stress4p(f4p), strain4p(d4p),
  stress1n(f1n), strain1n(d1n), stress2n(f2n), strain2n(d2n),
  stress3n(f3n), strain3n(d3n), stress4n(f4n), strain4n(d4n),
  envlpPosStress(6), envlpPosStrain(6), envlpNegStress(6), envlpNegStrain(6),
  gammaK1(gk1), gammaK2(gk2), gammaK3(gk3), gammaK4(gk4), gammaKLimit(gklim),
  gammaD1(gd1), gammaD2(gd2), gammaD3(gd3), gammaD4(gd4), gammaDLimit(gdlim),
  gammaF1(gf1), gammaF2(gf2), gammaF3(gf3), gammaF4(gf4), gammaFLimit(gflim),
  gammaE(ge), TnCycle(0.0), CnCycle(0.0), DmgCyc(dc),
  rDispP(mdp), rForceP(mfp), uForceP(msp), rDispN(mdn), rForceN(mfn), uForceN(msn),
  state3Stress(4), state3Strain(4), state4Stress(4), state4Strain(4), 
  envlpPosDamgdStress(6), envlpNegDamgdStress(6), Cenergy(0)
#ifdef _CSS
	 ,energy(0)
#endif // _CSS

{
	bool error = false;

	// Positive backbone parameters
	if (strain1p <= 0.0)
		error = true;

	if (strain2p <= 0.0)
		error = true;

	if (strain3p <= 0.0)
		error = true;

	if (strain4p <= 0.0)
		error = true;

	// Negative backbone parameters
	if (strain1n >= 0.0)
		error = true;

	if (strain2n >= 0.0)
		error = true;

	if (strain3n >= 0.0)
		error = true;

	if (strain4n >= 0.0)
		error = true;

	if (error){
		opserr << "ERROR: -- input backbone is not unique (one-to-one) , Pinching4Material::Pinching4Material" << "\a";
	}

	envlpPosStress.Zero(); envlpPosStrain.Zero(); envlpNegStress.Zero(); envlpNegStrain.Zero(); 
	energyCapacity = 0.0; kunload = 0.0; elasticStrainEnergy = 0.0;


	// set envelope slopes
	this->SetEnvelope();

	envlpPosDamgdStress = envlpPosStress; envlpNegDamgdStress = envlpNegStress;
	state3Stress.Zero(); state3Strain.Zero(); state4Stress.Zero(); state4Strain.Zero();
	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();
}

Pinching4Material::Pinching4Material(int tag,
				     double f1p, double d1p, double f2p, double d2p,
				     double f3p, double d3p, double f4p, double d4p,
				     double mdp, double mfp, double msp,
				     double gk1, double gk2, double gk3, double gk4, double gklim,
				     double gd1, double gd2, double gd3, double gd4, double gdlim,
				     double gf1, double gf2, double gf3, double gf4, double gflim, double ge, int dc):
  UniaxialMaterial(tag, MAT_TAG_Pinching4),
  stress1p(f1p), strain1p(d1p), stress2p(f2p), strain2p(d2p),
  stress3p(f3p), strain3p(d3p), stress4p(f4p), strain4p(d4p),
  envlpPosStress(6), envlpPosStrain(6), envlpNegStress(6), envlpNegStrain(6),
  gammaK1(gk1), gammaK2(gk2), gammaK3(gk3), gammaK4(gk4), gammaKLimit(gklim),
  gammaD1(gd1), gammaD2(gd2), gammaD3(gd3), gammaD4(gd4), gammaDLimit(gdlim),
  gammaF1(gf1), gammaF2(gf2), gammaF3(gf3), gammaF4(gf4), gammaFLimit(gflim),
  gammaE(ge), TnCycle(0.0), CnCycle(0.0), DmgCyc(dc),
  rDispP(mdp), rForceP(mfp), uForceP(msp),
  state3Stress(4), state3Strain(4), state4Stress(4), state4Strain(4), 
  envlpPosDamgdStress(6), envlpNegDamgdStress(6)
#ifdef _CSS
	 , energy(0)
#endif // _CSS

{
	bool error = false;

	// Positive backbone parameters
	if (strain1p <= 0.0)
		error = true;

	if (strain2p <= 0.0)
		error = true;

	if (strain3p <= 0.0)
		error = true;

	if (strain4p <= 0.0)
		error = true;

	if (error){
		opserr << "ERROR: -- input backbone is not unique (one-to-one) , Pinching4Material::Pinching4Material" << "\a";
	}

	strain1n = -strain1p; stress1n = -stress1p; strain2n = -strain2p; stress2n = -stress2p;
	strain3n = -strain3p; stress3n = -stress3p; strain4n = -strain4p; stress4n = -stress4p;
	rDispN = rDispP; rForceN = rForceP; uForceN = uForceP;

	envlpPosStress.Zero(); envlpPosStrain.Zero(); envlpNegStress.Zero(); envlpNegStrain.Zero(); 
	energyCapacity = 0.0; kunload = 0.0; elasticStrainEnergy = 0.0;
	state3Stress.Zero(); state3Strain.Zero(); state4Stress.Zero(); state4Strain.Zero();
	// set envelope slopes
	this->SetEnvelope();

	envlpPosDamgdStress = envlpPosStress; envlpNegDamgdStress = envlpNegStress;

	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();
}


Pinching4Material::Pinching4Material():
  UniaxialMaterial(0, MAT_TAG_Pinching4),
  stress1p(0.0), strain1p(0.0), stress2p(0.0), strain2p(0.0),
  stress3p(0.0), strain3p(0.0), stress4p(0.0), strain4p(0.0),
  stress1n(0.0), strain1n(0.0), stress2n(0.0), strain2n(0.0),
  stress3n(0.0), strain3n(0.0), stress4n(0.0), strain4n(0.0),
  gammaK1(0.0), gammaK2(0.0), gammaK3(0.0), gammaKLimit(0.0),
  gammaD1(0.0), gammaD2(0.0), gammaD3(0.0), gammaDLimit(0.0),
  gammaF1(0.0), gammaF2(0.0), gammaF3(0.0), gammaFLimit(0.0), gammaE(0.0),
  rDispP(0.0), rForceP(0.0), uForceP(0.0), rDispN(0.0), rForceN(0.0), uForceN(0.0)
#ifdef _CSS
	 , energy(0)
#endif // _CSS
{

}

Pinching4Material::~Pinching4Material()
{

}

int Pinching4Material::setTrialStrain(double strain, double CstrainRate)
{

	Tstate = Cstate;
	Tenergy = Cenergy;
	Tstrain = strain;
	lowTstateStrain = lowCstateStrain;
	hghTstateStrain = hghCstateStrain;
	lowTstateStress = lowCstateStress;
	hghTstateStress = hghCstateStress;
	TminStrainDmnd = CminStrainDmnd;
	TmaxStrainDmnd = CmaxStrainDmnd;
	TgammaF = CgammaF;
	TgammaK = CgammaK; 
	TgammaD = CgammaD;

	dstrain = Tstrain - Cstrain;
	if (dstrain<1e-12 && dstrain>-1e-12){
		dstrain = 0.0;
	}

	// determine new state if there is a change in state
	getstate(Tstrain,dstrain);  

	switch (Tstate)
	{ 

	case 0:
		Ttangent = envlpPosStress(0)/envlpPosStrain(0);
		Tstress = Ttangent*Tstrain;
		break;
	case 1:
		Tstress = posEnvlpStress(strain);
		Ttangent = posEnvlpTangent(strain);
		break;
	case 2:
		Ttangent = negEnvlpTangent(strain);
		Tstress = negEnvlpStress(strain);
		break;
	case 3:
		kunload = (hghTstateStrain<0.0) ? kElasticNegDamgd:kElasticPosDamgd; 	
			state3Strain(0) = lowTstateStrain;
			state3Strain(3) = hghTstateStrain;
			state3Stress(0) = lowTstateStress;
			state3Stress(3) = hghTstateStress;

		getState3(state3Strain,state3Stress,kunload);
		Ttangent = Envlp3Tangent(state3Strain,state3Stress,strain);
		Tstress = Envlp3Stress(state3Strain,state3Stress,strain);
		
		//Print(opserr,1);
		break;
	case 4:
		kunload = (lowTstateStrain<0.0) ? kElasticNegDamgd:kElasticPosDamgd;
			state4Strain(0) = lowTstateStrain;
			state4Strain(3) = hghTstateStrain;
			state4Stress(0) = lowTstateStress;
			state4Stress(3) = hghTstateStress;

		getState4(state4Strain,state4Stress,kunload);
		Ttangent = Envlp4Tangent(state4Strain,state4Stress,strain);
		Tstress = Envlp4Stress(state4Strain,state4Stress,strain);
		break;
	}

	double denergy = 0.5*(Tstress+Cstress)*dstrain;
	elasticStrainEnergy = (Tstrain>0.0) ? 0.5*Tstress/kElasticPosDamgd*Tstress:0.5*Tstress/kElasticNegDamgd*Tstress;

	Tenergy = Cenergy + denergy;

	updateDmg(Tstrain,dstrain);
	return 0;
}

double Pinching4Material::getStrain(void)
{
	return Tstrain;
}

double Pinching4Material::getStress(void)
{
	return Tstress;
}

double Pinching4Material::getTangent(void)
{
	return Ttangent;
}

double Pinching4Material::getInitialTangent(void)
{
	return envlpPosStress(0)/envlpPosStrain(0);
}

int Pinching4Material::commitState(void)  {
	Cstate = Tstate;

	if (dstrain>1e-12||dstrain<-(1e-12)) {
		CstrainRate = dstrain;}
	else {
		CstrainRate = TstrainRate;}

	lowCstateStrain = lowTstateStrain;
	lowCstateStress = lowTstateStress;
	hghCstateStrain = hghTstateStrain;
	hghCstateStress = hghTstateStress;
	CminStrainDmnd = TminStrainDmnd;
	CmaxStrainDmnd = TmaxStrainDmnd;
	Cenergy = Tenergy;

#ifdef _CSS
	energy += 0.5 * (Cstress + Tstress) * (Tstrain - Cstrain);
#endif // _CSS
	Cstress = Tstress;
	Cstrain = Tstrain;

	CgammaK = TgammaK;
	CgammaD = TgammaD;
	CgammaF = TgammaF;
	
	// define adjusted strength and stiffness parameters
	kElasticPosDamgd = kElasticPos*(1 - gammaKUsed);
	kElasticNegDamgd = kElasticNeg*(1 - gammaKUsed);

	uMaxDamgd = TmaxStrainDmnd*(1 + CgammaD);   
	uMinDamgd = TminStrainDmnd*(1 + CgammaD);

	envlpPosDamgdStress = envlpPosStress*(1-gammaFUsed);
	envlpNegDamgdStress = envlpNegStress*(1-gammaFUsed);

	CnCycle = TnCycle; // number of cycles of loading

	return 0;
}

int Pinching4Material::revertToLastCommit(void)
{
	
	Tstate = Cstate;

	TstrainRate = CstrainRate;

	lowTstateStrain = lowCstateStrain;
	lowTstateStress = lowCstateStress;
	hghTstateStrain = hghCstateStrain;
	hghTstateStress = hghCstateStress;
	TminStrainDmnd = CminStrainDmnd;
	TmaxStrainDmnd = CmaxStrainDmnd;
	Tenergy = Cenergy;

	Tstrain = Cstrain; Tstress = Cstress;

	TgammaD = CgammaD;
	TgammaK = CgammaK;
	TgammaF = CgammaF;

	TnCycle = CnCycle;

	return 0;
}

int Pinching4Material::revertToStart(void)
{
#ifdef _CSS
	 energy = 0;
#endif // _CSS
	 Cstate = 0;
	Cstrain = 0.0;
	Cstress = 0.0;
	CstrainRate = 0.0;

	this->SetEnvelope();
	
	lowCstateStrain = envlpNegStrain(0);
	lowCstateStress = envlpNegStress(0);
	hghCstateStrain = envlpPosStrain(0);
	hghCstateStress = envlpPosStress(0);
	CminStrainDmnd = envlpNegStrain(1);
	CmaxStrainDmnd = envlpPosStrain(1);
	Cenergy = 0.0;
	CgammaK = 0.0;
	CgammaD = 0.0;
	CgammaF = 0.0;
	CnCycle = 0.0;
	
	TnCycle = CnCycle;
	Tstrain = Cstrain;
	Tstress = Cstress;
	Ttangent = envlpPosStress(0)/envlpPosStrain(0);
	dstrain = 0.0;       
	gammaKUsed = 0.0;
	gammaFUsed = 0.0;

	state3Stress.Zero();
	state3Strain.Zero();
	state4Stress.Zero();
	state4Strain.Zero();	

	envlpPosDamgdStress = envlpPosStress;
	envlpNegDamgdStress = envlpNegStress;	

	kElasticPosDamgd = kElasticPos;
	kElasticNegDamgd = kElasticNeg;
	uMaxDamgd = CmaxStrainDmnd;
	uMinDamgd = CminStrainDmnd;

	return 0;
}

UniaxialMaterial* Pinching4Material::getCopy(void)
{
	Pinching4Material *theCopy = new Pinching4Material (this->getTag(),
		stress1p, strain1p, stress2p, strain2p, stress3p, strain3p, stress4p, strain4p,
        stress1n, strain1n, stress2n, strain2n, stress3n, strain3n, stress4n, strain4n,
        rDispP, rForceP, uForceP, rDispN, rForceN, uForceN,gammaK1,gammaK2,gammaK3,gammaK4,
		gammaKLimit,gammaD1,gammaD2,gammaD3,gammaD4,gammaDLimit,gammaF1,gammaF2,gammaF3,gammaF4,
		gammaFLimit,gammaE,DmgCyc);
	
#ifdef _CSS
	theCopy->energy = 0;
#endif // _CSS
	theCopy->rDispN = rDispN;
	theCopy->rDispP = rDispP;
	theCopy->rForceN = rForceN;
	theCopy->rForceP = rForceP;
	theCopy->uForceN = uForceN;
	theCopy->uForceP = uForceP;

	// Trial state variables
	theCopy->Tstress = Tstress;
	theCopy->Tstrain = Tstrain;
	theCopy->Ttangent = Ttangent;

	// Coverged material history parameters
	theCopy->Cstate = Cstate;
	theCopy->Cstrain = Cstrain;
	theCopy->Cstress = Cstress;
	theCopy->CstrainRate = CstrainRate;

	theCopy->lowCstateStrain = lowCstateStrain;
	theCopy->lowCstateStress = lowCstateStress;
	theCopy->hghCstateStrain = hghCstateStrain;
	theCopy->hghCstateStress = hghCstateStress;
	theCopy->CminStrainDmnd = CminStrainDmnd;
	theCopy->CmaxStrainDmnd = CmaxStrainDmnd;
	theCopy->Cenergy = Cenergy;
	theCopy->CgammaK = CgammaK;
	theCopy->CgammaD = CgammaD;
	theCopy->CgammaF = CgammaF;
	theCopy->CnCycle = CnCycle;
	theCopy->gammaKUsed = gammaKUsed;
	theCopy->gammaFUsed = gammaFUsed;
	theCopy->DmgCyc = DmgCyc;

	// trial material history parameters
	theCopy->Tstate = Tstate;
	theCopy->dstrain = dstrain;
	theCopy->lowTstateStrain = lowTstateStrain;
	theCopy->lowTstateStress = lowTstateStress;
	theCopy->hghTstateStrain = hghTstateStrain;
	theCopy->hghTstateStress = hghTstateStress;
	theCopy->TminStrainDmnd = TminStrainDmnd;
	theCopy->TmaxStrainDmnd = TmaxStrainDmnd;
	theCopy->Tenergy = Tenergy;
	theCopy->TgammaK = TgammaK;
	theCopy->TgammaD = TgammaD;
	theCopy->TgammaF = TgammaF;
	theCopy->TnCycle = TnCycle;

	// Strength and stiffness parameters
	theCopy->kElasticPos = kElasticPos;
	theCopy->kElasticNeg = kElasticNeg;
	theCopy->kElasticPosDamgd = kElasticPosDamgd;
	theCopy->kElasticNegDamgd = kElasticNegDamgd;
	theCopy->uMaxDamgd = uMaxDamgd;
	theCopy->uMinDamgd = uMinDamgd;

	for (int i = 0; i<6; i++)
	{
		theCopy->envlpPosStrain(i) = envlpPosStrain(i);
		theCopy->envlpPosStress(i) = envlpPosStress(i);
		theCopy->envlpNegStrain(i) = envlpNegStrain(i);
		theCopy->envlpNegStress(i) = envlpNegStress(i);
		theCopy->envlpNegDamgdStress(i) = envlpNegDamgdStress(i);
		theCopy->envlpPosDamgdStress(i) = envlpPosDamgdStress(i);
	}

	for (int j = 0; j<4; j++)
	{
		theCopy->state3Strain(j) = state3Strain(j);
		theCopy->state3Stress(j) = state3Stress(j);
		theCopy->state4Strain(j) = state4Strain(j);
		theCopy->state4Stress(j) = state4Stress(j);
	}

	theCopy->energyCapacity = energyCapacity;
	theCopy->kunload = kunload;
	theCopy->elasticStrainEnergy = elasticStrainEnergy;

	return theCopy;
}

int Pinching4Material::sendSelf(int commitTag, Channel &theChannel)
{

    // Instantiate a Vector to store the relevant class attributes
    static Vector data(119);

    // Fill the Vector with class attributes.
	int indx = 0;

	data(indx++) = this->getTag();

	// Backbone parameters
	data(indx++) = stress1p;
	data(indx++) = strain1p;
	data(indx++) = stress2p;
	data(indx++) = strain2p;
	data(indx++) = stress3p;
	data(indx++) = strain3p;
	data(indx++) = stress4p;
	data(indx++) = strain4p;
	data(indx++) = stress1n;
	data(indx++) = strain1n;
	data(indx++) = stress2n;
	data(indx++) = strain2n;
	data(indx++) = stress3n;
	data(indx++) = strain3n;
	data(indx++) = stress4n;
	data(indx++) = strain4n;

	// Damage parameters
	data(indx++) = gammaK1;
	data(indx++) = gammaK2;
	data(indx++) = gammaK3;
	data(indx++) = gammaK4;
	data(indx++) = gammaKLimit;
	data(indx++) = gammaD1;
	data(indx++) = gammaD2;
	data(indx++) = gammaD3;
	data(indx++) = gammaD4;
	data(indx++) = gammaDLimit;
	data(indx++) = gammaF1;
	data(indx++) = gammaF2;
	data(indx++) = gammaF3;
	data(indx++) = gammaF4;
	data(indx++) = gammaFLimit;
	data(indx++) = gammaE;
	data(indx++) = TnCycle;
	data(indx++) = CnCycle;
	data(indx++) = DmgCyc;

	// unloading-reloading parameters
	data(indx++) = rDispP;
	data(indx++) = rForceP;
	data(indx++) = uForceP;
	data(indx++) = rDispN;
	data(indx++) = rForceN;
	data(indx++) = uForceN;

	// Converged Material History parameters
	data(indx++) = Cstate;
	data(indx++) = Cstrain;
	data(indx++) = Cstress;
	data(indx++) = CstrainRate;
	data(indx++) = lowCstateStrain;
	data(indx++) = lowCstateStress;
	data(indx++) = hghCstateStrain;
	data(indx++) = hghCstateStress;
	data(indx++) = CminStrainDmnd;
	data(indx++) = CmaxStrainDmnd;
	data(indx++) = Cenergy;
	data(indx++) = CgammaK;
	data(indx++) = CgammaD;
	data(indx++) = CgammaF;
	data(indx++) = gammaKUsed;
	data(indx++) = gammaFUsed;

	// strength and stiffness parameters
	data(indx++) = kElasticPos;
	data(indx++) = kElasticNeg;
	data(indx++) = kElasticPosDamgd;
	data(indx++) = kElasticNegDamgd;
	data(indx++) = uMaxDamgd;
	data(indx++) = uMinDamgd;

	// energy parameters
	data(indx++) = energyCapacity;
	data(indx++) = kunload;
	data(indx++) = elasticStrainEnergy;

	// Vector-type variable contents
	for (int k = 0; k<6; k++){
	  data(indx++) = envlpPosStress[k];
	  data(indx++) = envlpPosStrain[k];
	  data(indx++) = envlpNegStress[k];
	  data(indx++) = envlpNegStrain[k];          
	  data(indx++) = envlpPosDamgdStress[k];
	  data(indx++) = envlpNegDamgdStress[k];
	}
	for (int k = 0; k<4; k++){
	  data(indx++) = state3Stress[k];
	  data(indx++) = state3Strain[k];
	  data(indx++) = state4Stress[k];
	  data(indx++) = state4Strain[k];
	}

	// Send the data vector
    int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0) 
    {
        opserr << "Pinching4Material::sendSelf() - failed to send data\n";
        return res;
    }

    return res;
}

int Pinching4Material::recvSelf(
    int commitTag, Channel &theChannel,
	FEM_ObjectBroker & theBroker)
{

	// Instantiate a Vector to store the relevant class attributes
    static Vector data(119);

    int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) 
    {
        opserr << "Pinching4Material::recvSelf() - failed to receive data\n";
        return res;
    }

    // Assign the received values to the class attributes

	int indx = 0;

	this->setTag((int)data(indx++));
	
	// Backbone parameters
	stress1p = data(indx++);
	strain1p = data(indx++);
	stress2p = data(indx++);
	strain2p = data(indx++);
	stress3p = data(indx++);
	strain3p = data(indx++);
	stress4p = data(indx++);
	strain4p = data(indx++);
	stress1n = data(indx++);
	strain1n = data(indx++);
	stress2n = data(indx++);
	strain2n = data(indx++);
	stress3n = data(indx++);
	strain3n = data(indx++);
	stress4n = data(indx++);
	strain4n = data(indx++);

	// Damage parameters
	gammaK1 = data(indx++);
	gammaK2 = data(indx++);
	gammaK3 = data(indx++);
	gammaK4 = data(indx++);
	gammaKLimit = data(indx++);
	gammaD1 = data(indx++);
	gammaD2 = data(indx++);
	gammaD3 = data(indx++);
	gammaD4 = data(indx++);
	gammaDLimit = data(indx++);
	gammaF1 = data(indx++);
	gammaF2 = data(indx++);
	gammaF3 = data(indx++);
	gammaF4 = data(indx++);
	gammaFLimit = data(indx++);
	gammaE = data(indx++);
	TnCycle = data(indx++);
	CnCycle = data(indx++);
	DmgCyc = (int)data(indx++);

	// unloading-reloading parameters
	rDispP = data(indx++);
	rForceP = data(indx++);
	uForceP = data(indx++);
	rDispN = data(indx++);
	rForceN = data(indx++);
	uForceN = data(indx++);

	// Converged Material History parameters
	Cstate = (int)data(indx++);
	Cstrain = data(indx++);
	Cstress = data(indx++);
	CstrainRate = data(indx++);
	lowCstateStrain = data(indx++);
	lowCstateStress = data(indx++);
	hghCstateStrain = data(indx++);
	hghCstateStress = data(indx++);
	CminStrainDmnd = data(indx++);
	CmaxStrainDmnd = data(indx++);
	Cenergy = data(indx++);
	CgammaK = data(indx++);
	CgammaD = data(indx++);
	CgammaF = data(indx++);
	gammaKUsed = data(indx++);
	gammaFUsed = data(indx++);

	// strength and stiffness parameters
	kElasticPos = data(indx++);
	kElasticNeg = data(indx++);
	kElasticPosDamgd = data(indx++);
	kElasticNegDamgd = data(indx++);
	uMaxDamgd = data(indx++);
	uMinDamgd = data(indx++);

	// energy parameters
	energyCapacity = data(indx++);
	kunload = data(indx++);
	elasticStrainEnergy = data(indx++);

	// Vector-type variable contents
	for (int k = 0; k<6; k++){
	  envlpPosStress[k] = data(indx++);
	  envlpPosStrain[k] = data(indx++);
	  envlpNegStress[k] = data(indx++);
	  envlpNegStrain[k] = data(indx++);
	  envlpPosDamgdStress[k] = data(indx++);
	  envlpNegDamgdStress[k] = data(indx++);
	}
	for (int k = 0; k<4; k++){
	  state3Stress[k] = data(indx++);
	  state3Strain[k] = data(indx++);
	  state4Stress[k] = data(indx++);
	  state4Strain[k] = data(indx++);
	}

	// Set trial variables to the last converged values
	Tstress = Cstress;
	Tstrain = Cstrain;
	Ttangent = 0.0;
	Tstate = Cstate;
	dstrain = 0.0;
	TstrainRate = CstrainRate;
	lowTstateStrain = lowCstateStrain;
	lowTstateStress = lowCstateStress;
	hghTstateStrain = hghCstateStrain;
	hghTstateStress = hghCstateStress;
	TminStrainDmnd = CminStrainDmnd;
	TmaxStrainDmnd = CmaxStrainDmnd;
	Tenergy = Cenergy;
	TgammaK = CgammaK;
	TgammaD = CgammaD;
	TgammaF = CgammaF;
	
    return res;
}

void Pinching4Material::Print(OPS_Stream &s, int flag)
{
	// Create a string stream in order to flush at the end everything together (useful for parallel running)
	std::stringstream ss;
	ss << "Pinching4Material, tag: " << this->getTag() << endln;
	ss << "strain: " << Tstrain << endln;
	ss << "stress: " << Tstress << endln;
	ss << "state: " << Tstate << endln;
	if (flag == 2) {
		ss << "Envelope positive: (" << this->strain1p << ", " << this->stress1p << "), ";
		ss << "(" << this->strain2p << ", " << this->stress2p << "), ";
		ss << "(" << this->strain3p << ", " << this->stress3p << "), ";
		ss << "(" << this->strain4p << ", " << this->stress4p << ")\n";
		ss << "Envelope negative: (" << this->strain1n << ", " << this->stress1n << "), ";
		ss << "(" << this->strain2n << ", " << this->stress2n << "), ";
		ss << "(" << this->strain3n << ", " << this->stress3n << "), ";
		ss << "(" << this->strain4n << ", " << this->stress4n << ")\n";
		ss << "Hysteretic parameters: " << endln;
		ss << "rDispP = " << rDispP << endln;
		ss << "rForceP = " << rForceP << endln;
		ss << "uForceP = " << uForceP << endln;
		ss << "rDispN = " << rDispN << endln;
		ss << "rForceN = " << rForceN << endln;
		ss << "uForceN = " << uForceN << endln;
	}
	else if (flag == 3) {
		ss << "Cstrain: " << Cstrain << endln;
		ss << "Cstress:" << Cstress << endln;
		ss << "Cstate: " << Cstate << endln;
		ss << "Tstrain: " << Tstrain << endln;
		ss << "Cstrainrate: " << CstrainRate << endln;
		ss << "dstrain: " << dstrain << endln;
		int cid = 0;
		if (dstrain * CstrainRate <= 0.0) {
			cid = 1;
		}
		ss << "change in direction? " << dstrain * CstrainRate << "( " << cid << ")\n";
		ss << "Tstate: " << Tstate << endln;
		ss << "plt.plot([" << this->strain4n << ", " << this->strain3n << ", " << this->strain2n << ", " << this->strain1n << ", ";
		ss << this->strain1p << ", " << this->strain2p << ", " << this->strain3p << ", " << this->strain4p << "], [";
		ss << this->stress4n << ", " << this->stress3n << ", " << this->stress2n << ", " << this->stress1n << ", ";
		ss << this->stress1p << ", " << this->stress2p << ", " << this->stress3p << ", " << this->stress4p << "],'-o')\n";
		ss << "points st4:\n";
		ss << "plt.plot([" << this->state4Strain(0) << ", " << this->state4Strain(1) << ", " << this->state4Strain(2) << ", " << this->state4Strain(3) << "], ["; 
		ss << this->state4Stress(0) << ", " << this->state4Stress(1) << ", " << this->state4Stress(2) << ", " << this->state4Stress(3) << "],'-o')\n";
		ss << "Commited point:\n";
		ss << "plt.plot(" << this->Cstrain << ", " << this->Cstress << ", 's')\n";
		ss << "Trial point:\n";
		ss << "plt.plot(" << this->Tstrain << ", " << this->Tstress << ", '^')\n";
	}

	std::string str = ss.str();
	s << str.c_str();

}

void Pinching4Material::SetEnvelope(void)
{ 
	double kPos = stress1p/strain1p;
	double kNeg = stress1n/strain1n;
	double k = (kPos>kNeg) ? kPos:kNeg;
	double u = (strain1p>-strain1n) ? 1e-4*strain1p:-1e-4*strain1n;
    
	envlpPosStrain(0) = u;
	envlpPosStress(0) = u*k;
	envlpNegStrain(0) = -u;
	envlpNegStress(0) = -u*k;

	envlpPosStrain(1) = strain1p;
	envlpPosStrain(2) = strain2p;
	envlpPosStrain(3) = strain3p;
	envlpPosStrain(4) = strain4p;

	envlpNegStrain(1) = strain1n;
	envlpNegStrain(2) = strain2n;
	envlpNegStrain(3) = strain3n;
	envlpNegStrain(4) = strain4n;

	envlpPosStress(1) = stress1p;
	envlpPosStress(2) = stress2p;
	envlpPosStress(3) = stress3p;
	envlpPosStress(4) = stress4p;

	envlpNegStress(1) = stress1n;
	envlpNegStress(2) = stress2n;
	envlpNegStress(3) = stress3n;
	envlpNegStress(4) = stress4n;

	double k1 = (stress4p - stress3p)/(strain4p - strain3p);
	double k2 = (stress4n - stress3n)/(strain4n - strain3n);


	envlpPosStrain(5) = 1e+6*strain4p;
	envlpPosStress(5) = (k1>0.0)? stress4p+k1*(envlpPosStrain(5) - strain4p):stress4p*1.1;
	envlpNegStrain(5) = 1e+6*strain4n;
	envlpNegStress(5) = (k2>0.0)? stress4n+k2*(envlpNegStrain(5) - strain4n):stress4n*1.1;
	
	// define critical material properties
	kElasticPos = envlpPosStress(1)/envlpPosStrain(1);	
	kElasticNeg = envlpNegStress(1)/envlpNegStrain(1);

	double energypos = 0.5*envlpPosStrain(0)*envlpPosStress(0);

	for (int jt = 0; jt<4; jt++){
		energypos += 0.5*(envlpPosStress(jt) + envlpPosStress(jt+1))*(envlpPosStrain(jt+1)-envlpPosStrain(jt));
	}

	double energyneg = 0.5*envlpNegStrain(0)*envlpNegStress(0);

	for (int jy = 0; jy<4; jy++){
		energyneg += 0.5*(envlpNegStress(jy) + envlpNegStress(jy+1))*(envlpNegStrain(jy+1)-envlpNegStrain(jy));
	}

	double max_energy = (energypos>energyneg) ? energypos:energyneg;

	energyCapacity = gammaE*max_energy;

	
}

void Pinching4Material::getstate(double u,double du)
{
	int cid = 0;
	int cis = 0;
	int newState = 0; 
	if (du*CstrainRate<=0.0){   
		cid = 1;
	}
	if (u<lowTstateStrain || u>hghTstateStrain || cid) {                
		if (Tstate == 0) {                                              
			if (u>hghTstateStrain) {
				cis = 1;
				newState = 1;
				lowTstateStrain = envlpPosStrain(0);
				lowTstateStress = envlpPosStress(0);
				hghTstateStrain = envlpPosStrain(5);
				hghTstateStress = envlpPosStress(5);
			}
			else if (u<lowTstateStrain){
				cis = 1;
				newState = 2;
				lowTstateStrain = envlpNegStrain(5);
				lowTstateStress = envlpNegStress(5);
				hghTstateStrain = envlpNegStrain(0);
				hghTstateStress = envlpNegStress(0);
			}
		}
		else if (Tstate==1 && du<0.0) {
			cis = 1;
			if (Cstrain>TmaxStrainDmnd) {
				TmaxStrainDmnd = u - du;
			}
			if (TmaxStrainDmnd<uMaxDamgd) {
				TmaxStrainDmnd = uMaxDamgd;
			}
			if (u<uMinDamgd) {
				newState = 2;
				gammaFUsed = CgammaF;     
				for (int i=0; i<=5; i++) {
					envlpNegDamgdStress(i) = envlpNegStress(i)*(1.0-gammaFUsed);
				}
				lowTstateStrain = envlpNegStrain(5);
				lowTstateStress = envlpNegStress(5);
				hghTstateStrain = envlpNegStrain(0);
				hghTstateStress = envlpNegStress(0);
			}
			else {
				newState = 3;
				lowTstateStrain = uMinDamgd;
				gammaFUsed = CgammaF;        
				for (int i=0; i<=5; i++) {
					envlpNegDamgdStress(i) = envlpNegStress(i)*(1.0-gammaFUsed);
				}
				lowTstateStress = negEnvlpStress(uMinDamgd); 
				hghTstateStrain = Cstrain;
				hghTstateStress = Cstress;
			}
			gammaKUsed = CgammaK;
			kElasticPosDamgd = kElasticPos*(1.0-gammaKUsed);
		}
		else if (Tstate ==2 && du>0.0){
			cis = 1;
			if (Cstrain<TminStrainDmnd) {
				TminStrainDmnd = Cstrain;
			}
			if (TminStrainDmnd>uMinDamgd) {
				TminStrainDmnd = uMinDamgd;
			}
			if (u>uMaxDamgd) {
				newState = 1;
				gammaFUsed = CgammaF;      
				for (int i=0; i<=5; i++) {
					envlpPosDamgdStress(i) = envlpPosStress(i)*(1.0-gammaFUsed);
				}
				lowTstateStrain = envlpPosStrain(0);
				lowTstateStress = envlpPosStress(0);
				hghTstateStrain = envlpPosStrain(5);
				hghTstateStress = envlpPosStress(5);
			}
			else {
				newState = 4;
				lowTstateStrain = Cstrain;
				lowTstateStress = Cstress;
				hghTstateStrain = uMaxDamgd;
				gammaFUsed = CgammaF;         
				for (int i=0; i<=5; i++) {
					envlpPosDamgdStress(i) = envlpPosStress(i)*(1.0-gammaFUsed);
				}
				hghTstateStress = posEnvlpStress(uMaxDamgd);
			}
			gammaKUsed = CgammaK;
			kElasticNegDamgd = kElasticNeg*(1.0-gammaKUsed);
		}
			else if (Tstate ==3) {
				if (u<lowTstateStrain){
					cis = 1;
					newState = 2;
					lowTstateStrain = envlpNegStrain(5);
					hghTstateStrain = envlpNegStrain(0);
					lowTstateStress = envlpNegDamgdStress(5);
					hghTstateStress = envlpNegDamgdStress(0);
				}
				else if (u>uMaxDamgd && du>0.0) {
					cis = 1;
					newState = 1;
					lowTstateStrain = envlpPosStrain(0);
					lowTstateStress = envlpPosStress(0);
					hghTstateStrain = envlpPosStrain(5);
					hghTstateStress = envlpPosStress(5);
				}
				else if (du>0.0) {
					cis = 1;
					newState = 4;
					lowTstateStrain = Cstrain;
					lowTstateStress = Cstress;
					hghTstateStrain = uMaxDamgd;
					gammaFUsed = CgammaF;
					for (int i=0; i<=5; i++) {
					envlpPosDamgdStress(i) = envlpPosStress(i)*(1.0-gammaFUsed);
				    }
					hghTstateStress = posEnvlpStress(uMaxDamgd);
					gammaKUsed = CgammaK;
					kElasticNegDamgd = kElasticNeg*(1.0-gammaKUsed); 
				}
			}
			else if (Tstate == 4){
				if (u>hghTstateStrain){
					cis = 1;
					newState = 1;
					lowTstateStrain = envlpPosStrain(0);
					lowTstateStress = envlpPosDamgdStress(0);
					hghTstateStrain = envlpPosStrain(5);
					hghTstateStress = envlpPosDamgdStress(5);
				}
				else if (u<uMinDamgd && du <0.0) {
					cis = 1;
					newState = 2;
					lowTstateStrain = envlpNegStrain(5);
					lowTstateStress = envlpNegDamgdStress(5);
					hghTstateStrain = envlpNegStrain(0);
					hghTstateStress = envlpNegDamgdStress(0);
				}
				else if (du<0.0) { 
					cis = 1;
					newState = 3;
					lowTstateStrain = uMinDamgd;
					gammaFUsed = CgammaF;         
					for (int i=0; i<=5; i++) {
					envlpNegDamgdStress(i) = envlpNegStress(i)*(1.0-gammaFUsed);
				     }
					lowTstateStress = negEnvlpStress(uMinDamgd);
					hghTstateStrain = Cstrain;
					hghTstateStress = Cstress;
					gammaKUsed = CgammaK;
					kElasticPosDamgd = kElasticPos*(1.0-gammaKUsed);
				}
			}
			}
			if (cis) {
				Tstate = newState;
				
			}
		}

 double Pinching4Material::posEnvlpStress(double u)
		{
			double k = 0.0;
			int i = 0;
			double f = 0.0;
			while (k==0.0 && i<=4){   
				 
				 if (u<=envlpPosStrain(i+1)){
					 k = (envlpPosDamgdStress(i+1)-envlpPosDamgdStress(i))/(envlpPosStrain(i+1)-envlpPosStrain(i));
					 f = envlpPosDamgdStress(i) + (u-envlpPosStrain(i))*k;
				 }
				 i++;
			}
            

			if (k==0.0){
				k = (envlpPosDamgdStress(5) - envlpPosDamgdStress(4))/(envlpPosStrain(5) - envlpPosStrain(4));
				f = envlpPosDamgdStress(5) + k*(u-envlpPosStrain(5));
			}

			return f;

		}

double Pinching4Material::posEnvlpTangent(double u)
	 {
			double k = 0.0;
			int i = 0;
			while (k==0.0 && i<=4){        
				 
				 if (u<=envlpPosStrain(i+1)){
					 k = (envlpPosDamgdStress(i+1)-envlpPosDamgdStress(i))/(envlpPosStrain(i+1)-envlpPosStrain(i));
			}
				 i++;
			}

			if (k==0.0){
				k = (envlpPosDamgdStress(5) - envlpPosDamgdStress(4))/(envlpPosStrain(5) - envlpPosStrain(4));
		   }

			return k;

		}


 double Pinching4Material::negEnvlpStress(double u)
		{
			double k = 0.0;
			int i = 0;
			double f = 0.0;
			while (k==0.0 && i<=4){      				 
				 if (u>=envlpNegStrain(i+1)){
					 k = (envlpNegDamgdStress(i)-envlpNegDamgdStress(i+1))/(envlpNegStrain(i)-envlpNegStrain(i+1));
					 f = envlpNegDamgdStress(i+1)+(u-envlpNegStrain(i+1))*k;
				 }
				 i++;
			}

			if (k==0.0){
				k = (envlpNegDamgdStress(4) - envlpNegDamgdStress(5))/(envlpNegStrain(4)-envlpNegStrain(5));
				f = envlpNegDamgdStress(5) + k*(u-envlpNegStrain(5));
			}
			return f;

		}

double Pinching4Material::negEnvlpTangent(double u)
		{
			double k = 0.0;
			int i = 0;
			while (k==0.0 && i<=4){              
				 
				 if (u>=envlpNegStrain(i+1)){
					 k = (envlpNegDamgdStress(i)-envlpNegDamgdStress(i+1))/(envlpNegStrain(i)-envlpNegStrain(i+1));
					}
				 i++;
			}

			if (k==0.0){
				k = (envlpNegDamgdStress(4) - envlpNegDamgdStress(5))/(envlpNegStrain(4)-envlpNegStrain(5));
				}
			return k;

		}


void Pinching4Material::getState3(Vector& state3Strain, Vector& state3Stress, double kunload)
		{

			double kmax = (kunload>kElasticNegDamgd) ? kunload:kElasticNegDamgd;

			if (state3Strain(0)*state3Strain(3) <0.0){
				// trilinear unload reload path expected, first define point for reloading
				state3Strain(1) = lowTstateStrain*rDispN;
				if (rForceN-uForceN > 1e-8) {
					state3Stress(1) = lowTstateStress*rForceN;
				}
				else {
					if (TminStrainDmnd < envlpNegStrain(3)) {
						double st1 = lowTstateStress*uForceN*(1.0+1e-6);
						double st2 = envlpNegDamgdStress(4)*(1.0+1e-6);
						state3Stress(1) = (st1<st2) ? st1:st2;
					}
					else {
						double st1 = envlpNegDamgdStress(3)*uForceN*(1.0+1e-6);
						double st2 = envlpNegDamgdStress(4)*(1.0+1e-6);
						state3Stress(1) = (st1<st2) ? st1:st2;
					}
				}
				// if reload stiffness exceeds unload stiffness, reduce reload stiffness to make it equal to unload stiffness
				if ((state3Stress(1)-state3Stress(0))/(state3Strain(1)-state3Strain(0)) > kElasticNegDamgd) {
					state3Strain(1) = lowTstateStrain + (state3Stress(1)-state3Stress(0))/kElasticNegDamgd;
				}
				// check that reloading point is not behind point 4 
				if (state3Strain(1)>state3Strain(3)) {
					// path taken to be a straight line between points 1 and 4
					double du = state3Strain(3) - state3Strain(0);
					double df = state3Stress(3) - state3Stress(0);
					state3Strain(1) = state3Strain(0) + 0.33*du;
					state3Strain(2) = state3Strain(0) + 0.67*du;
					state3Stress(1) = state3Stress(0) + 0.33*df;
					state3Stress(2) = state3Stress(0) + 0.67*df;
				}
				else {
					if (TminStrainDmnd < envlpNegStrain(3)) {
						state3Stress(2) = uForceN*envlpNegDamgdStress(4);
					}
					else {
						state3Stress(2) = uForceN*envlpNegDamgdStress(3);
					}
					state3Strain(2) = hghTstateStrain - (hghTstateStress-state3Stress(2))/kunload;

					if (state3Strain(2) > state3Strain(3)) {
						// point 3 should be along a line between 2 and 4
						double du = state3Strain(3) - state3Strain(1);
						double df = state3Stress(3) - state3Stress(1);
						state3Strain(2) = state3Strain(1) + 0.5*du;
						state3Stress(2) = state3Stress(1) + 0.5*df;
					}
					else if ((state3Stress(2) - state3Stress(1))/(state3Strain(2) - state3Strain(1)) > kmax) {
						// linear unload-reload path expected
						double du = state3Strain(3) - state3Strain(0);
						double df = state3Stress(3) - state3Stress(0);
						state3Strain(1) = state3Strain(0) + 0.33*du;
						state3Strain(2) = state3Strain(0) + 0.67*du;
						state3Stress(1) = state3Stress(0) + 0.33*df;
						state3Stress(2) = state3Stress(0) + 0.67*df;
					}
					else if ((state3Strain(2) < state3Strain(1))||((state3Stress(2)-state3Stress(1))/(state3Strain(2)-state3Strain(1))<0)) {
						if (state3Strain(2)<0.0) {
							// pt 3 should be along a line between 2 and 4
							double du = state3Strain(3)-state3Strain(1);
							double df = state3Stress(3)-state3Stress(1);
							state3Strain(2) = state3Strain(1) + 0.5*du;
							state3Stress(2) = state3Stress(1) + 0.5*df;
						}
						else if (state3Strain(1) > 0.0) {
							// pt 2 should be along a line between 1 and 3
							double du = state3Strain(2)-state3Strain(0);
							double df = state3Stress(2)-state3Stress(0);
							state3Strain(1) = state3Strain(0) + 0.5*du;
							state3Stress(1) = state3Stress(0) + 0.5*df;
						}
						else {
							double avgforce = 0.5*(state3Stress(2) + state3Stress(1));
							double dfr = 0.0;
							if (avgforce < 0.0){
								dfr = -avgforce/100;
							}
							else {
								dfr = avgforce/100;
							}
							double slope12 = (state3Stress(1) - state3Stress(0))/(state3Strain(1) - state3Strain(0));
							double slope34 = (state3Stress(3) - state3Stress(2))/(state3Strain(3) - state3Strain(2));
							state3Stress(1) = avgforce - dfr;
							state3Stress(2) = avgforce + dfr;
							state3Strain(1) = state3Strain(0) + (state3Stress(1) - state3Stress(0))/slope12;
							state3Strain(2) = state3Strain(3) - (state3Stress(3) - state3Stress(2))/slope34;
						}
					}
				}
			}
				else {
					// linear unload reload path is expected		 
					double du = state3Strain(3)-state3Strain(0);
					double df = state3Stress(3)-state3Stress(0);
					state3Strain(1) = state3Strain(0) + 0.33*du;
					state3Strain(2) = state3Strain(0) + 0.67*du;
					state3Stress(1) = state3Stress(0) + 0.33*df;
					state3Stress(2) = state3Stress(0) + 0.67*df;
				}
			
				
				double checkSlope = state3Stress(0)/state3Strain(0);
				double slope = 0.0;

				// final check
				int i = 0;
				while (i<3) {
					double du = state3Strain(i+1)-state3Strain(i);
					double df = state3Stress(i+1)-state3Stress(i);
					if (du<0.0 || df<0.0) {
						double du = state3Strain(3)-state3Strain(0);
						double df = state3Stress(3)-state3Stress(0);
						state3Strain(1) = state3Strain(0) + 0.33*du;
						state3Strain(2) = state3Strain(0) + 0.67*du;
						state3Stress(1) = state3Stress(0) + 0.33*df;
						state3Stress(2) = state3Stress(0) + 0.67*df;
						slope = df/du;
						i = 3;
					}
					if (slope > 1e-8 && slope < checkSlope) {
						state3Strain(1) = 0.0; state3Stress(1) = 0.0;
						state3Strain(2) = state3Strain(3)/2; state3Stress(2) = state3Stress(3)/2;
					} 
					i++;
				}


			}
				
void Pinching4Material::getState4(Vector& state4Strain,Vector& state4Stress, double kunload)
		{

			double kmax = (kunload>kElasticPosDamgd) ? kunload:kElasticPosDamgd;

			if (state4Strain(0)*state4Strain(3) <0.0){
				// trilinear unload reload path expected
				state4Strain(2) = hghTstateStrain*rDispP;
				if (uForceP==0.0){
					state4Stress(2) = hghTstateStress*rForceP;
				}
				else if (rForceP-uForceP > 1e-8) {
					state4Stress(2) = hghTstateStress*rForceP;
				}
				else {
					if (TmaxStrainDmnd > envlpPosStrain(3)) {
						double st1 = hghTstateStress*uForceP*(1.0+1e-6);
						double st2 = envlpPosDamgdStress(4)*(1.0+1e-6);
						state4Stress(2) = (st1>st2) ? st1:st2;
					}
					else {
						double st1 = envlpPosDamgdStress(3)*uForceP*(1.0+1e-6);
						double st2 = envlpPosDamgdStress(4)*(1.0+1e-6);
						state4Stress(2) = (st1>st2) ? st1:st2;
					}
				}
				// if reload stiffness exceeds unload stiffness, reduce reload stiffness to make it equal to unload stiffness
				if ((state4Stress(3)-state4Stress(2))/(state4Strain(3)-state4Strain(2)) > kElasticPosDamgd) {
					state4Strain(2) = hghTstateStrain - (state4Stress(3)-state4Stress(2))/kElasticPosDamgd;
				}
				// check that reloading point is not behind point 1 
				if (state4Strain(2)<state4Strain(0)) {
					// path taken to be a straight line between points 1 and 4
					double du = state4Strain(3) - state4Strain(0);
					double df = state4Stress(3) - state4Stress(0);
					state4Strain(1) = state4Strain(0) + 0.33*du;
					state4Strain(2) = state4Strain(0) + 0.67*du;
					state4Stress(1) = state4Stress(0) + 0.33*df;
					state4Stress(2) = state4Stress(0) + 0.67*df;
				}
				else {
					if (TmaxStrainDmnd > envlpPosStrain(3)) {
						state4Stress(1) = uForceP*envlpPosDamgdStress(4);
					}
					else {
						state4Stress(1) = uForceP*envlpPosDamgdStress(3);
					}
					state4Strain(1) = lowTstateStrain + (-lowTstateStress+state4Stress(1))/kunload;

					if (state4Strain(1) < state4Strain(0)) {
						// point 2 should be along a line between 1 and 3
						double du = state4Strain(2) - state4Strain(0);
						double df = state4Stress(2) - state4Stress(0);
						state4Strain(1) = state4Strain(0) + 0.5*du;
						state4Stress(1) = state4Stress(0) + 0.5*df;
					}
					else if ((state4Stress(2) - state4Stress(1))/(state4Strain(2) - state4Strain(1)) > kmax) {
						// linear unload-reload path expected
						double du = state4Strain(3) - state4Strain(0);
						double df = state4Stress(3) - state4Stress(0);
						state4Strain(1) = state4Strain(0) + 0.33*du;
						state4Strain(2) = state4Strain(0) + 0.67*du;
						state4Stress(1) = state4Stress(0) + 0.33*df;
						state4Stress(2) = state4Stress(0) + 0.67*df;
					}
					else if ((state4Strain(2) < state4Strain(1))||((state4Stress(2)-state4Stress(1))/(state4Strain(2)-state4Strain(1))<0)) {
						if (state4Strain(1)>0.0) {
							// pt 2 should be along a line between 1 and 3
							double du = state4Strain(2)-state4Strain(0);
							double df = state4Stress(2)-state4Stress(0);
							state4Strain(1) = state4Strain(0) + 0.5*du;
							state4Stress(1) = state4Stress(0) + 0.5*df;
						}
						else if (state4Strain(2) < 0.0) {
							// pt 2 should be along a line between 2 and 4
							double du = state4Strain(3)-state4Strain(1);
							double df = state4Stress(3)-state4Stress(1);
							state4Strain(2) = state4Strain(1) + 0.5*du;
							state4Stress(2) = state4Stress(1) + 0.5*df;
						}
						else {
							double avgforce = 0.5*(state4Stress(2) + state4Stress(1));
							double dfr = 0.0;
							if (avgforce < 0.0){
								dfr = -avgforce/100;
							}
							else {
								dfr = avgforce/100;
							}
							double slope12 = (state4Stress(1) - state4Stress(0))/(state4Strain(1) - state4Strain(0));
							double slope34 = (state4Stress(3) - state4Stress(2))/(state4Strain(3) - state4Strain(2));
							state4Stress(1) = avgforce - dfr;
							state4Stress(2) = avgforce + dfr;
							state4Strain(1) = state4Strain(0) + (state4Stress(1) - state4Stress(0))/slope12;
							state4Strain(2) = state4Strain(3) - (state4Stress(3) - state4Stress(2))/slope34;
						}
					}
				}
			}
				else {
					// linear unload reload path is expected
					double du = state4Strain(3)-state4Strain(0);
					double df = state4Stress(3)-state4Stress(0);
					state4Strain(1) = state4Strain(0) + 0.33*du;
					state4Strain(2) = state4Strain(0) + 0.67*du;
					state4Stress(1) = state4Stress(0) + 0.33*df;
					state4Stress(2) = state4Stress(0) + 0.67*df;
				}
			

				
				double checkSlope = state4Stress(0)/state4Strain(0);
				double slope = 0.0;

				// final check
				int i = 0;
				while (i<3) {
					double du = state4Strain(i+1)-state4Strain(i);
					double df = state4Stress(i+1)-state4Stress(i);
					if (du<0.0 || df<0.0) {
						double du = state4Strain(3)-state4Strain(0);
						double df = state4Stress(3)-state4Stress(0);
						state4Strain(1) = state4Strain(0) + 0.33*du;
						state4Strain(2) = state4Strain(0) + 0.67*du;
						state4Stress(1) = state4Stress(0) + 0.33*df;
						state4Stress(2) = state4Stress(0) + 0.67*df;
						slope = df/du;
						i = 3;
					}
					if (slope > 1e-8 && slope < checkSlope) {
						state4Strain(1) = 0.0; state4Stress(1) = 0.0;
						state4Strain(2) = state4Strain(3)/2; state4Stress(2) = state4Stress(3)/2;
					} 

					i++;
				}
			}

double Pinching4Material::Envlp3Tangent(Vector s3Strain, Vector s3Stress, double u)
			{
				double k = 0.0;
				int i = 0;
				while ((k==0.0||i<=2) && (i<=2)) 
				{
					if (u>= s3Strain(i)) {
						k = (s3Stress(i+1)-s3Stress(i))/(s3Strain(i+1)-s3Strain(i));
					}
					i++;
				}
				if (k==0.0) {
					if (u<s3Strain(0)) {
						i = 0;
					}
					else {
						i = 2;
					}
					k = (s3Stress(i+1)-s3Stress(i))/(s3Strain(i+1)-s3Strain(i));
					
				}
				return k;
			}

double Pinching4Material::Envlp4Tangent(Vector s4Strain, Vector s4Stress, double u)
			{
				double k = 0.0;
				int i = 0;
				while ((k==0.0||i<=2) && (i<=2)) 
				{
					if (u>= s4Strain(i)) {
						k = (s4Stress(i+1)-s4Stress(i))/(s4Strain(i+1)-s4Strain(i));
					}
					i++;
				}
				if (k==0.0) {
					if (u<s4Strain(0)) {
						i = 0;
					}
					else {
						i = 2;
					}
					k = (s4Stress(i+1)-s4Stress(i))/(s4Strain(i+1)-s4Strain(i));
					
				}
				return k;
			}


  double Pinching4Material::Envlp3Stress(Vector s3Strain, Vector s3Stress, double u)
			{
				double k = 0.0;
				int i = 0;
				double f = 0.0;
				while ((k==0.0||i<=2) && (i<=2)) 
				{
					if (u>= s3Strain(i)) {
						k = (s3Stress(i+1)-s3Stress(i))/(s3Strain(i+1)-s3Strain(i));
						f = s3Stress(i)+(u-s3Strain(i))*k;
					}
					i++;
				}
				if (k==0.0) {
					if (u<s3Strain(0)) {
						i = 0;
					}
					else {
						i = 2;
					}
					k = (s3Stress(i+1)-s3Stress(i))/(s3Strain(i+1)-s3Strain(i));
					f = s3Stress(i)+(u-s3Strain(i))*k;
				}
				return f;
			}

double Pinching4Material::Envlp4Stress(Vector s4Strain, Vector s4Stress, double u)
			{
				double k = 0.0;
				int i = 0;
				double f = 0.0;
				while ((k==0.0||i<=2) && (i<=2)) 
				{
					if (u>= s4Strain(i)) {
						k = (s4Stress(i+1)-s4Stress(i))/(s4Strain(i+1)-s4Strain(i));
						f = s4Stress(i)+(u-s4Strain(i))*k;
					}
					i++;
				}
				if (k==0.0) {
					if (u<s4Strain(0)) {
						i = 0;
					}
					else {
						i = 2;
					}
					k = (s4Stress(i+1)-s4Stress(i))/(s4Strain(i+1)-s4Strain(i));
					f = s4Stress(i)+(u-s4Strain(i))*k;
				}
				return f;
			}

void Pinching4Material::updateDmg(double strain, double dstrain)
	{
		double tes = 0.0;
		double umaxAbs = (TmaxStrainDmnd>-TminStrainDmnd) ? TmaxStrainDmnd:-TminStrainDmnd;
		double uultAbs = (envlpPosStrain(4)>-envlpNegStrain(4)) ? envlpPosStrain(4):-envlpNegStrain(4);
		TnCycle = CnCycle + fabs(dstrain)/(4*umaxAbs);
		if ((strain<uultAbs && strain>-uultAbs)&& Tenergy< energyCapacity)
		{
			TgammaK = gammaK1*pow((umaxAbs/uultAbs),gammaK3);
			TgammaD = gammaD1*pow((umaxAbs/uultAbs),gammaD3);
			TgammaF = gammaF1*pow((umaxAbs/uultAbs),gammaF3);

			if (Tenergy>elasticStrainEnergy && DmgCyc == 0) {
				tes = ((Tenergy-elasticStrainEnergy)/energyCapacity);
				TgammaK = TgammaK + gammaK2*pow(tes,gammaK4);
				TgammaD = TgammaD + gammaD2*pow(tes,gammaD4);
				TgammaF = TgammaF + gammaF2*pow(tes,gammaF4);
			} else if (DmgCyc == 1) {
				TgammaK = TgammaK + gammaK2*pow(TnCycle,gammaK4);
				TgammaD = TgammaD + gammaD2*pow(TnCycle,gammaD4);
				TgammaF = TgammaF + gammaF2*pow(TnCycle,gammaF4);
			}
			double kminP = (posEnvlpStress(TmaxStrainDmnd)/TmaxStrainDmnd);
			double kminN = (negEnvlpStress(TminStrainDmnd)/TminStrainDmnd);
			double kmin = ((kminP/kElasticPos)>(kminN/kElasticNeg)) ? (kminP/kElasticPos):(kminN/kElasticNeg);
			double gammaKLimEnv = (0.0>(1.0-kmin)) ? 0.0:(1.0-kmin);
			
			double k1 = (TgammaK<gammaKLimit) ? TgammaK:gammaKLimit;
			TgammaK = (k1<gammaKLimEnv) ? k1:gammaKLimEnv;
			TgammaD = (TgammaD<gammaDLimit) ? TgammaD:gammaDLimit;
			TgammaF = (TgammaF<gammaFLimit) ? TgammaF:gammaFLimit;
		}
		else if (strain<uultAbs && strain>-uultAbs) {
			double kminP = (posEnvlpStress(TmaxStrainDmnd)/TmaxStrainDmnd);
			double kminN = (negEnvlpStress(TminStrainDmnd)/TminStrainDmnd);
			double kmin = ((kminP/kElasticPos)>=(kminN/kElasticNeg)) ? (kminP/kElasticPos):(kminN/kElasticNeg);
			double gammaKLimEnv = (0.0>(1.0-kmin)) ? 0.0:(1.0-kmin);
			
			TgammaK = (gammaKLimit<gammaKLimEnv) ? gammaKLimit:gammaKLimEnv;    
			TgammaD = gammaDLimit;
			TgammaF = gammaFLimit;
		}
		
	}

int
Pinching4Material::setParameter(const char** argv, int argc, Parameter& param)
{
	// Parameters for backbone control points
	if (strcmp(argv[0], "f1p") == 0 || strcmp(argv[0], "stress1p") == 0) {
		param.setValue(stress1p);
		return param.addObject(1, this);
	}
	if (strcmp(argv[0], "d1p") == 0 || strcmp(argv[0], "strain1p") == 0) {
		param.setValue(strain1p);
		return param.addObject(2, this);
	}
	if (strcmp(argv[0], "f2p") == 0 || strcmp(argv[0], "stress2p") == 0) {
		param.setValue(stress2p);
		return param.addObject(3, this);
	}
	if (strcmp(argv[0], "d2p") == 0 || strcmp(argv[0], "strain2p") == 0) {
		param.setValue(strain2p);
		return param.addObject(4, this);
	}
	if (strcmp(argv[0], "f3p") == 0 || strcmp(argv[0], "stress3p") == 0) {
		param.setValue(stress3p);
		return param.addObject(5, this);
	}
	if (strcmp(argv[0], "d3p") == 0 || strcmp(argv[0], "strain3p") == 0) {
		param.setValue(strain3p);
		return param.addObject(6, this);
	}
	if (strcmp(argv[0], "f4p") == 0 || strcmp(argv[0], "stress4p") == 0) {
		param.setValue(stress4p);
		return param.addObject(7, this);
	}
	if (strcmp(argv[0], "d4p") == 0 || strcmp(argv[0], "strain4p") == 0) {
		param.setValue(strain4p);
		return param.addObject(8, this);
	}
	if (strcmp(argv[0], "f1n") == 0 || strcmp(argv[0], "stress1n") == 0) {
		param.setValue(stress1n);
		return param.addObject(9, this);
	}
	if (strcmp(argv[0], "d1n") == 0 || strcmp(argv[0], "strain1n") == 0) {
		param.setValue(strain1n);
		return param.addObject(10, this);
	}
	if (strcmp(argv[0], "f2n") == 0 || strcmp(argv[0], "stress2n") == 0) {
		param.setValue(stress2n);
		return param.addObject(11, this);
	}
	if (strcmp(argv[0], "d2n") == 0 || strcmp(argv[0], "strain2n") == 0) {
		param.setValue(strain2n);
		return param.addObject(12, this);
	}
	if (strcmp(argv[0], "f3n") == 0 || strcmp(argv[0], "stress3n") == 0) {
		param.setValue(stress3n);
		return param.addObject(13, this);
	}
	if (strcmp(argv[0], "d3n") == 0 || strcmp(argv[0], "strain3n") == 0) {
		param.setValue(strain3n);
		return param.addObject(14, this);
	}
	if (strcmp(argv[0], "f4n") == 0 || strcmp(argv[0], "stress4n") == 0) {
		param.setValue(stress4n);
		return param.addObject(15, this);
	}
	if (strcmp(argv[0], "d4n") == 0 || strcmp(argv[0], "strain4n") == 0) {
		param.setValue(strain4n);
		return param.addObject(16, this);
	}

	// Parameters for hysteretic rules
	if (strcmp(argv[0], "rDispP") == 0) {
		param.setValue(rDispP);
		return param.addObject(17, this);
	}
	if (strcmp(argv[0], "rForceP") == 0) {
		param.setValue(rForceP);
		return param.addObject(18, this);
	}
	if (strcmp(argv[0], "uForceP") == 0) {
		param.setValue(uForceP);
		return param.addObject(19, this);
	}
	if (strcmp(argv[0], "rDispN") == 0) {
		param.setValue(rDispN);
		return param.addObject(20, this);
	}
	if (strcmp(argv[0], "rForceN") == 0) {
		param.setValue(rForceN);
		return param.addObject(21, this);
	}
	if (strcmp(argv[0], "uForceN") == 0) {
		param.setValue(uForceN);
		return param.addObject(22, this);
	}

	return -1;
}



int
Pinching4Material::updateParameter(int parameterID, Information& info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->stress1p = info.theDouble;
		break;
	case 2:
		this->strain1p = info.theDouble;
		break;
	case 3:
		this->stress2p = info.theDouble;
		break;
	case 4:
		this->strain2p = info.theDouble;
		break;
	case 5:
		this->stress3p = info.theDouble;
		break;
	case 6:
		this->strain3p = info.theDouble;
		break;
	case 7:
		this->stress4p = info.theDouble;
		break;
	case 8:
		this->strain4p = info.theDouble;
		break;
	case 9:
		this->stress1n = info.theDouble;
		break;
	case 10:
		this->strain1n = info.theDouble;
		break;
	case 11:
		this->stress2n = info.theDouble;
		break;
	case 12:
		this->strain2n = info.theDouble;
		break;
	case 13:
		this->stress3n = info.theDouble;
		break;
	case 14:
		this->strain3n = info.theDouble;
		break;
	case 15:
		this->stress4n = info.theDouble;
		break;
	case 16:
		this->strain4n = info.theDouble;
		break;
	case 17:
		this->rDispP = info.theDouble;
		break;
	case 18:
		this->rForceP = info.theDouble;
		break;
	case 19:
		this->uForceP = info.theDouble;
		break;
	case 20:
		this->rDispN = info.theDouble;
		break;
	case 21:
		this->rForceN = info.theDouble;
		break;
	case 22:
		this->uForceN = info.theDouble;
		break;
	default:
		return -1;
	}

	// Changed a parameter: we need to update the envelope
	this->SetEnvelope();

	// Then we need to force the update of damaged envelope
	envlpPosDamgdStress = envlpPosStress * (1 - gammaFUsed);
	envlpNegDamgdStress = envlpNegStress * (1 - gammaFUsed);

	// setTrialStrain and the commit to store history variables
	this->setTrialStrain(this->getStrain());
	this->commitState();


	return 0;
}
