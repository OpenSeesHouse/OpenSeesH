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
                                                                        
// $Revision: 1.19 $
// $Date: 2008-12-18 23:40:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/BucklingMaterial.cpp,v $

// Written: MHS
// Created: July 2000
//
// Description: This file contains the implementation of 
// BucklingMaterial.  BucklingMaterial is
// a one-dimensional hysteretic model with pinching of both
// force and deformation, damage due to deformation and energy, and
// degraded unloading stiffness based on maximum ductility.  This
// is a modified implementation of Hyster2.f90 by Filippou.
#include <OPS_Globals.h>
#ifdef _CSS
#include <stdlib.h>
#include <BucklingMaterial.h>
#include <math.h>
#include <float.h>
#include <Channel.h>

#include <elementAPI.h>

void *
OPS_BucklingMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs != 20 && numArgs != 19 && numArgs != 16 && numArgs != 15) {
    opserr << "Want: uniaxialMaterial BucklingMaterial tag? mom1p? rot1p? mom2p? rot2p? <mom3p? rot3p?> "
	   << "\nmom1n? rot1n? mom2n? rot2n? <mom3n? rot3n?> pinchX? pinchY? damfc1pp? damfc2pp? damfc1pn? damfc2pn? <beta?>";
    return 0;
  }
  
  int iData[1];
  double dData[19];
  for (int i=0; i<19; i++) 
    dData[i] = 0.0;

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial BucklingMaterial" << endln;
    return 0;
  }

  numData = numArgs-1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxial BucklingMaterial " << iData[0] << endln;
    return 0;	
  }

  // Parsing was successful, allocate the material
  if (numData > 15) 
    theMaterial = new BucklingMaterial(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
					 dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
  					 dData[13], dData[14], dData[15], dData[16], dData[17], dData[18]);
  else
    theMaterial = new BucklingMaterial(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
					 dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type BucklingMaterial\n";
    return 0;
  }

  return theMaterial;
}



BucklingMaterial::BucklingMaterial(int tag,
				       double m1p, double r1p, double m2p, double r2p, double m3p, double r3p,
				       double m1n, double r1n, double m2n, double r2n, double m3n, double r3n,
				       double px, double py, double d1p, double d2p, double d1n, double d2n, double b):
UniaxialMaterial(tag, MAT_TAG_Hysteretic),
pinchX(px), pinchY(py), damfc1p(d1p), damfc2p(d2p), damfc1n(d1n), damfc2n(d2n), beta(b),
mom1p(m1p), rot1p(r1p), mom2p(m2p), rot2p(r2p), mom3p(m3p), rot3p(r3p),
mom1n(m1n), rot1n(r1n), mom2n(m2n), rot2n(r2n), mom3n(m3n), rot3n(r3n)
{
	 Tdamfcn = 0;
  bool error = false;
  // Positive backbone parameters
  if (rot1p <= 0.0)
    error = true;
  
  if (rot2p <= rot1p)
    error = true;
  
  if (rot3p <= rot2p)
    error = true;
  
  // Negative backbone parameters
  if (rot1n >= 0.0)
    error = true;
  
  if (rot2n >= rot1n)
    error = true;
  
  if (rot3n >= rot2n)
    error = true;
  
  if (error) {
    opserr << "BucklingMaterial::BucklingMaterial -- input backbone is not unique (one-to-one)\n";
    exit(-1);
  }		
  

  energyA = 0.5 * (rot1p*mom1p + (rot2p-rot1p)*(mom2p+mom1p) + (rot3p-rot2p)*(mom3p+mom2p) +
		   rot1n*mom1n + (rot2n-rot1n)*(mom2n+mom1n) + (rot3n-rot2n)*(mom3n+mom2n));
  
  // Set envelope slopes
  this->setEnvelope();
  
  // Initialize history variables
  this->revertToStart();
  this->revertToLastCommit();
  
}

BucklingMaterial::BucklingMaterial(int tag,
			double m1p, double r1p, double m2p, double r2p,
			double m1n, double r1n, double m2n, double r2n,
			double px, double py, double d1, double d2, double d1n, double d2n, double b):
UniaxialMaterial(tag, MAT_TAG_Hysteretic),
pinchX(px), pinchY(py), damfc1p(d1), damfc2p(d2), damfc1n(d1n), damfc2n(d2n), beta(b),
mom1p(m1p), rot1p(r1p), mom3p(m2p), rot3p(r2p),
mom1n(m1n), rot1n(r1n), mom3n(m2n), rot3n(r2n)
{
	 Tdamfcn = 0;
	 bool error = false;
	
	// Positive backbone parameters
	if (rot1p <= 0.0)
		error = true;

	if (rot3p <= rot1p)
		error = true;

	// Negative backbone parameters
	if (rot1n >= 0.0)
		error = true;

	if (rot3n >= rot1n)
		error = true;

	if (error) {
	  opserr << "BucklingMaterial::BucklingMaterial -- input backbone is not unique (one-to-one)\n";
	  exit(-1);
	}

				      

	energyA = 0.5 * (rot1p*mom1p + (rot3p-rot1p)*(mom3p+mom1p) +
		rot1n*mom1n + (rot3n-rot1n)*(mom3n+mom1n));

	mom2p = 0.5*(mom1p+mom3p);
	mom2n = 0.5*(mom1n+mom3n);

	rot2p = 0.5*(rot1p+rot3p);
	rot2n = 0.5*(rot1n+rot3n);

	// Set envelope slopes
	this->setEnvelope();

	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();
}

BucklingMaterial::BucklingMaterial():
UniaxialMaterial(0, MAT_TAG_Hysteretic),
pinchX(0.0), pinchY(0.0), damfc1p(0.0), damfc2p(0.0), beta(0.0),
mom1p(0.0), rot1p(0.0), mom2p(0.0), rot2p(0.0), mom3p(0.0), rot3p(0.0),
mom1n(0.0), rot1n(0.0), mom2n(0.0), rot2n(0.0), mom3n(0.0), rot3n(0.0)
{
	 Tdamfcn = 0;

}

BucklingMaterial::~BucklingMaterial()
{
	// Nothing to do
}

int
BucklingMaterial::setTrialStrain(double strain, double strainRate)
{
	 if (TloadIndicator == 0 && strain == 0.0)
		  return 0;

	 TrotMax = CrotMax;
	 TrotMin = CrotMin;
	 TenergyD = CenergyD;
	 TrotPu = CrotPu;
	 TrotNu = CrotNu;
	 Tdamfcn = Cdamfcn;

	 Tstrain = strain;
	 double dStrain = Tstrain - Cstrain;

	 if (fabs(dStrain) < DBL_EPSILON)
		  return 0;

	 TloadIndicator = CloadIndicator;

	 if (TloadIndicator == 0)
		  TloadIndicator = (dStrain < 0.0) ? 2 : 1;

	 if (Tstrain >= CrotMax) {
		  TrotMax = Tstrain;
		  Ttangent = posEnvlpTangent(Tstrain);
		  Tstress = posEnvlpStress(Tstrain);
		  TloadIndicator = 1;
	 }
	 else if (Tstrain <= CrotMin) {
		  TrotMin = Tstrain;
		  Ttangent = negEnvlpTangent(Tstrain-TrotPu);
		  Tstress = negEnvlpStress(Tstrain - TrotPu);
		  TloadIndicator = 2;
	 }
	 else {
		  if (dStrain < 0.0)
				negativeIncrement(dStrain);
		  else if (dStrain > 0.0)
				positiveIncrement(dStrain);
	 }

	 TenergyD = CenergyD + 0.5 * (Cstress + Tstress) * dStrain;


	 //  if (this->getTag() == 40)
	 //    opserr << "setTrial: " << Tstrain << " " << Ttangent << " " << Tstress << endln;

	 return 0;
}


double
BucklingMaterial::getStrain(void)
{
	return Tstrain;
}

double
BucklingMaterial::getStress(void)
{
	return Tstress;
}

double
BucklingMaterial::getTangent(void)
{
  return Ttangent;
}

void
BucklingMaterial::positiveIncrement(double dStrain)
{
	double kn = pow((CrotMin- TrotPu )/rot1n,beta);
	kn = (kn < 1.0) ? 1.0 : 1.0/kn;
	double kp = pow(CrotMax/rot1p,beta);
	kp = (kp < 1.0) ? 1.0 : 1.0/kp;

	if (TloadIndicator == 2) {
		TloadIndicator = 1;
		if (Cstress <= 0.0) {
			TrotNu = Cstrain - Cstress/(Eun*kn);
			double energy = CenergyD - 0.5*Cstress/(Eun*kn)*Cstress;
			double damfc = 0.0;
			if (CrotMin- TrotPu < rot1n) {
				damfc = damfc2p*energy/energyA;
				damfc += damfc1p*(CrotMin- TrotPu -rot1n)/rot1n;
			}

			TrotMax *= (1.0+damfc);
		}
	}

  TloadIndicator = 1;

  if (TrotMax > POS_INF_STRAIN)
    TrotMax = POS_INF_STRAIN;

	TrotMax = (TrotMax > rot1p) ? TrotMax : rot1p;

	double maxmom = posEnvlpStress(TrotMax);
	double rotlim = negEnvlpRotlim(CrotMin- TrotPu);
	double rotrel = (rotlim > TrotNu) ? rotlim : TrotNu;

	// rotrel = TrotNu;
	// if (negEnvlpStress(CrotMin) >= 0.0)
	//    rotrel = rotlim;
	
	//	double rotmp1 = rotrel + pinchY*(TrotMax-rotrel);

	double rotmp2 = TrotMax - (1.0-pinchY)*maxmom/(Eup*kp);
	//double rotmp2 = TrotMax-(1-pinchY)*maxmom/Eup;
	//	double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;
	double rotch = rotrel + (rotmp2-rotrel)*pinchX;                   // changed on 7/11/2006

	double tmpmo1;
	double tmpmo2;

	if (Tstrain < TrotNu) {
		Ttangent = Eun*kn;
		Tstress = Cstress + Ttangent*dStrain;
		if (Tstress >= 0.0) {
			Tstress = 0.0;
			Ttangent = Eun*1.0e-9;
		}
	}

	else if (Tstrain < rotch) {
		if (Tstrain <= rotrel) {
			Tstress = 0.0;
			Ttangent = Eup*1.0e-9;
		}
		else {
			Ttangent = maxmom*pinchY/(rotch-rotrel);
			tmpmo1 = Cstress + Eup*kp*dStrain;
			tmpmo2 = (Tstrain-rotrel)*Ttangent;
			if (tmpmo1 < tmpmo2) {
				Tstress = tmpmo1;
				Ttangent = Eup*kp;
			}
			else
				Tstress = tmpmo2;
		}
	}

	else {
		Ttangent = (1.0-pinchY)*maxmom/(TrotMax-rotch);
		tmpmo1 = Cstress + Eup*kp*dStrain;
		tmpmo2 = pinchY*maxmom + (Tstrain-rotch)*Ttangent;
		if (tmpmo1 < tmpmo2) {
			Tstress = tmpmo1;
			Ttangent = Eup*kp;
		}
		else
			Tstress = tmpmo2;
	}
}

void
BucklingMaterial::negativeIncrement(double dStrain)
{
	 double kn = pow((CrotMin- TrotPu) / rot1n, beta);
	 kn = (kn < 1.0) ? 1.0 : 1.0 / kn;
	 double kp = pow(CrotMax / rot1p, beta);
	 kp = (kp < 1.0) ? 1.0 : 1.0 / kp;
	 if (TloadIndicator == 1) {
		  TloadIndicator = 2;
		  if (Cstress >= 0.0) {
				TrotPu = Cstrain - Cstress / (Eup * kp);
				double energy = CenergyD - 0.5 * Cstress / (Eup * kp) * Cstress;

				if (CrotMax > rot1p) {
					 Tdamfcn = damfc2n * energy / energyA;
					 Tdamfcn += damfc1n * (CrotMax - rot1p) / rot1p;
				}


				TrotMin = rot1n + TrotPu;
		  }
	 }
	 else
		  TrotMin = rot1n * (1.0 + Tdamfcn) + TrotPu;
	 double py = 0.999, px = 0.5;
	 TloadIndicator = 2;

	 if (TrotMin < NEG_INF_STRAIN)
		  TrotMin = NEG_INF_STRAIN;

	 TrotMin = (TrotMin < rot1n + TrotPu) ? TrotMin : rot1n + TrotPu;
	 if (Tstrain < TrotMin)
		  TrotMin = Tstrain;
	 double minmom = negEnvlpStress(TrotMin - TrotPu);
	 double rotlim = posEnvlpRotlim(CrotMax);
	 double rotrel = (rotlim < TrotPu) ? rotlim : TrotPu;

	 //rotrel = TrotPu;
	 //if (posEnvlpStress(CrotMax) <= 0.0)
	 //  rotrel = rotlim;

	 //double rotmp1 = rotrel + py*(TrotMin-rotrel);
	 double rotmp2 = TrotMin - (1.0 - py) * minmom / (Eun * kn);
	 //double rotmp2 = TrotMin-(1-py)*minmom/Eun;	
	 //double rotch = rotmp1 + (rotmp2-rotmp1)*px;
	 double rotch = rotrel + (rotmp2 - rotrel) * px;                   // changed on 7/11/2006
	 double tmpmo1;
	 double tmpmo2;

	 if (Tstrain > TrotPu) {
		  Ttangent = Eup * kp;
		  Tstress = Cstress + Ttangent * dStrain;
		  if (Tstress <= 0.0) {
				Tstress = 0.0;
				Ttangent = Eup * 1.0e-9;
		  }
	 }

	 else if (Tstrain <= TrotPu && Tstrain > rotch) {
		  if (Tstrain >= rotrel) {
				Tstress = 0.0;
				Ttangent = Eun * 1.0e-9;
		  }
		  else {
				Ttangent = minmom * py / (rotch - rotrel);
				tmpmo1 = Cstress + Eun * kn * dStrain;
				tmpmo2 = (Tstrain - rotrel) * Ttangent;
				if (tmpmo1 > tmpmo2) {
					 Tstress = tmpmo1;
					 Ttangent = Eun * kn;
				}
				else
					 Tstress = tmpmo2;
		  }
		  //Tstress *= (1. - Tdamfcn);
		  //Ttangent *= (1. - Tdamfcn);
	 }

	 else {
		  Ttangent = (1.0 - py) * minmom / (TrotMin - rotch);
		  tmpmo1 = Cstress + Eun * kn * dStrain;
		  tmpmo2 = py * minmom + (Tstrain - rotch) * Ttangent;
		  if (tmpmo1 > tmpmo2) {
				Tstress = tmpmo1;
				Ttangent = Eun * kn;
		  }
		  else
				Tstress = tmpmo2;
		  //Tstress *= (1. - Tdamfcn);
		  //Ttangent *= (1. - Tdamfcn);
	 }
}

int
BucklingMaterial::commitState(void)
{
	CrotMax = TrotMax;
	CrotMin = TrotMin;
	CrotPu = TrotPu;
	CrotNu = TrotNu;
	CenergyD = TenergyD;
	CloadIndicator = TloadIndicator;

	Cstress = Tstress;
	Cstrain = Tstrain;
	Cdamfcn = Tdamfcn;

	return 0;
}

int
BucklingMaterial::revertToLastCommit(void)
{
	TrotMax = CrotMax;
	TrotMin = CrotMin;
	TrotPu = CrotPu;
	TrotNu = CrotNu;
	TenergyD = CenergyD;
	TloadIndicator = CloadIndicator;

	Tstress = Cstress;
	Tstrain = Cstrain;
	Tdamfcn = Cdamfcn;

	return 0;
}

int
BucklingMaterial::revertToStart(void)
{
	CrotMax = 0.0;
	CrotMin = 0.0;
	CrotPu = 0.0;
	CrotNu = 0.0;
	CenergyD = 0.0;
	CloadIndicator = 0;

	Cstress = 0.0;
	Cstrain = 0.0;

	Tstrain = 0;
	Tstress = 0;
	Ttangent = E1p;

	return 0;
}

UniaxialMaterial*
BucklingMaterial::getCopy(void)
{
	BucklingMaterial *theCopy = new BucklingMaterial (this->getTag(),
		mom1p, rot1p, mom2p, rot2p, mom3p, rot3p,
		mom1n, rot1n, mom2n, rot2n, mom3n, rot3n,
		pinchX, pinchY, damfc1p, damfc2p, damfc1n, damfc2n, beta);

	theCopy->CrotMax = CrotMax;
	theCopy->CrotMin = CrotMin;
	theCopy->CrotPu = CrotPu;
	theCopy->CrotNu = CrotNu;
	theCopy->CenergyD = CenergyD;
	theCopy->CloadIndicator = CloadIndicator;
	theCopy->Cstress = Cstress;
	theCopy->Cstrain = Cstrain;
	theCopy->Ttangent = Ttangent;

	return theCopy;
}

int
BucklingMaterial::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(27);
  
  data(0) = this->getTag();
  data(1) = mom1p;
  data(2) = rot1p;
  data(3) = mom2p;
  data(4) = rot2p;
  data(5) = mom3p;
  data(6) = rot3p;
  data(7) = mom1n;
  data(8) = rot1n;
  data(9) = mom2n;
  data(10) = rot2n;
  data(11) = mom3n;
  data(12) = rot3n;
  data(13) = pinchX;
  data(14) = pinchY;
  data(15) = damfc1p;
  data(16) = damfc2p;
  data(17) = beta;
  data(18) = CrotMax;
  data(19) = CrotMin;
  data(20) = CrotPu;
  data(21) = CrotNu;
  data(22) = CenergyD;
  data(23) = CloadIndicator;
  data(24) = Cstress;
  data(25) = Cstrain;
  data(26) = Ttangent;

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) 
    opserr << "BucklingMaterial::sendSelf() - failed to send data\n";


  return res;
}

int
BucklingMaterial::recvSelf(int commitTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(27);
  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
      opserr << "BucklingMaterial::recvSelf() - failed to receive data\n";
      return res;
  }
  else {
    this->setTag((int)data(0));
    mom1p = data(1);
    rot1p = data(2);
    mom2p = data(3);
    rot2p = data(4);
    mom3p = data(5);
    rot3p = data(6);
    mom1n = data(7);
    rot1n = data(8);
    mom2n = data(9);
    rot2n = data(10);
    mom3n = data(11);
    rot3n = data(12);
    pinchX = data(13);
    pinchY = data(14);
    damfc1p = data(15);
    damfc2p = data(16);
    beta = data(17);

    CrotMax = data(18);
    CrotMin = data(19);
    CrotPu = data(20);
    CrotNu = data(21);
    CenergyD = data(22);
    CloadIndicator = int(data(23));
    Cstress = data(24);
    Cstrain = data(25);
    Ttangent = data(26);

    // set the trial values
    TrotMax = CrotMax;
    TrotMin = CrotMin;
    TrotPu = CrotPu;
    TrotNu = CrotNu;
    TenergyD = CenergyD;
    TloadIndicator = CloadIndicator;
    Tstress = Cstress;
    Tstrain = Cstrain;

  }

  // Set envelope slopes
  this->setEnvelope();
  
  return 0;
}
    
void
BucklingMaterial::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "HBucklingMaterial, tag: " << this->getTag() << endln;
        s << "s1p: " << mom1p << endln;
        s << "e1p: " << rot1p << endln;
        s << "E1p: " << E1p << endln;
        s << "s2p: " << mom2p << endln;
        s << "e2p: " << rot2p << endln;
        s << "E2p: " << E2p << endln;
        s << "s3p: " << mom3p << endln;
        s << "e3p: " << rot3p << endln;
        s << "E3p: " << E3p << endln;
        
        s << "s1n: " << mom1n << endln;
        s << "e1n: " << rot1n << endln;
        s << "E1n: " << E1n << endln;
        s << "s2n: " << mom2n << endln;
        s << "e2n: " << rot2n << endln;
        s << "E2n: " << E2n << endln;
        s << "s3n: " << mom3n << endln;
        s << "e3n: " << rot3n << endln;
        s << "E3n: " << E3n << endln;
        
        s << "pinchX: " << pinchX << endln;
        s << "pinchY: " << pinchY << endln;
        s << "damfc1p: " << damfc1p << endln;
        s << "damfc2p: " << damfc2p << endln;
        s << "energyA: " << energyA << endln;
        s << "beta: " << beta << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"BucklingMaterial\", ";
        s << "\"s1p\": " << mom1p << ", ";
        s << "\"e1p\": " << rot1p << ", ";
        s << "\"E1p\": " << E1p << ", ";
        s << "\"s2p\": " << mom2p << ", ";
        s << "\"e2p\": " << rot2p << ", ";
        s << "\"E2p\": " << E2p << ", ";
        s << "\"s3p\": " << mom3p << ", ";
        s << "\"e3p\": " << rot3p << ", ";
        s << "\"E3p\": " << E3p << ", ";
        
        s << "\"s1n\": " << mom1n << ", ";
        s << "\"e1n\": " << rot1n << ", ";
        s << "\"E1n\": " << E1n << ", ";
        s << "\"s2n\": " << mom2n << ", ";
        s << "\"e2n\": " << rot2n << ", ";
        s << "\"E2n\": " << E2n << ", ";
        s << "\"s3n\": " << mom3n << ", ";
        s << "\"e3n\": " << rot3n << ", ";
        s << "\"E3n\": " << E3n << ", ";
        
        s << "\"pinchX\": " << pinchX << ", ";
        s << "\"pinchY\": " << pinchY << ", ";
        s << "\"damfc1p\": " << damfc1p << ", ";
        s << "\"damfc2p\": " << damfc2p << ", ";
        s << "\"energyA\": " << energyA << ", ";
        s << "\"beta\": " << beta << "}";
    }
}

void
BucklingMaterial::setEnvelope(void)
{
	E1p = mom1p/rot1p;
	E2p = (mom2p-mom1p)/(rot2p-rot1p);
	E3p = (mom3p-mom2p)/(rot3p-rot2p);

	E1n = mom1n/rot1n;
	E2n = (mom2n-mom1n)/(rot2n-rot1n);
	E3n = (mom3n-mom2n)/(rot3n-rot2n);

	Eup = E1p;
	if (E2p > Eup) Eup = E2p;
	if (E3p > Eup) Eup = E3p;

	Eun = E1n;
	if (E2n > Eun) Eun = E2n;
	if (E3n > Eun) Eun = E3n;
}

double
BucklingMaterial::posEnvlpStress(double strain)
{
	if (strain <= 0.0)
		return 0.0;
	else if (strain <= rot1p)
		return E1p*strain;
	else if (strain <= rot2p)
		return mom1p + E2p*(strain-rot1p);
	else if (strain <= rot3p || E3p > 0.0)
		return mom2p + E3p*(strain-rot2p);
	else
		return mom3p;
}

double
BucklingMaterial::negEnvlpStress(double strain)
{
	if (strain >= 0)
		return 0.0;
	else if (strain >= rot1n)
		return E1n*strain;
	else if (strain >= rot2n)
		return mom1n + E2n*(strain-rot1n);
	else if (strain >= rot3n || E3n > 0.0)
		return mom2n + E3n*(strain-rot2n);
	else
		return mom3n;
}

double
BucklingMaterial::posEnvlpTangent(double strain)
{
  if (strain < 0.0)
    return E1p*1.0e-9;
  else if (strain <= rot1p)
    return E1p;
  else if (strain <= rot2p)
    return E2p;
  else if (strain <= rot3p || E3p > 0.0)
    return E3p;
  else
    return E1p*1.0e-9;
}

double
BucklingMaterial::negEnvlpTangent(double strain)
{
  if (strain > 0)
    return E1n*1.0e-9;
  else if (strain >= rot1n)
    return E1n;
  else if (strain >= rot2n)
    return E2n;
  else if (strain >= rot3n || E3n > 0.0)
    return E3n;
  else
    return E1n*1.0e-9;
}

double
BucklingMaterial::posEnvlpRotlim(double strain)
{
  double strainLimit = POS_INF_STRAIN;

  if (strain <= rot1p)
    return POS_INF_STRAIN;
  if (strain > rot1p && strain <= rot2p && E2p < 0.0)
    strainLimit = rot1p - mom1p/E2p;
  if (strain > rot2p && E3p < 0.0)
    strainLimit = rot2p - mom2p/E3p;

  if (strainLimit == POS_INF_STRAIN)
    return POS_INF_STRAIN;
  else if (posEnvlpStress(strainLimit) > 0)
    return POS_INF_STRAIN;
  else
    return strainLimit;
}

double
BucklingMaterial::negEnvlpRotlim(double strain)
{
  double strainLimit = NEG_INF_STRAIN;

  if (strain >= rot1n)
    return NEG_INF_STRAIN;
  if (strain < rot1n && strain >= rot2n && E2n < 0.0)
    strainLimit = rot1n - mom1n/E2n;
  if (strain < rot2n && E3n < 0.0)
    strainLimit = rot2n - mom2n/E3n;

  if (strainLimit == NEG_INF_STRAIN)
    return NEG_INF_STRAIN;
  else if (negEnvlpStress(strainLimit) < 0)
    return NEG_INF_STRAIN;
  else
    return strainLimit;
}

#endif //_CSS
