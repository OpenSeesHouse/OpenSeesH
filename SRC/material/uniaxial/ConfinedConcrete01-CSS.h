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
// $Date: 2010-02-10 23:31:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ConfinedConcrete01.h,v $
                                                                        
#ifndef ConfinedConcrete01_h
#define ConfinedConcrete01_h

// Description: This file contains the class definition for ConfinedConcrete01.
// Compressive envelope curve:     concrete model for confined concrete proposed by F. Braga, R. Gigliotti, M. Laterza in  
//                                 "Analytical Stress-Strain Relationship for Concrete Confined by Steel Stirrups
//                                 and/or FRP Jackets", ASCE, Journal of Structural Engineering, Vol. 132, No. 9, 
//                                 September 1, 2006 
// Unloading and reloanding curve: stress-strain model proposed by D. Karsan, and J. Jirsa in 
//                                 "Behavior of concrete under compressive loadings"
//                                 Journal of the Structural Division, ASCE, Vol. 95, No. ST12, December, 1969
// Tensile envelope strength:      model has no tensile strength
//                                 


// Written: Michele D'Amato, University of Basilicata, Potenza, Italy -email: damato.mic@gmail.com
//          Newton Le, University of California, Davis, CA, USA -email: newle@ucdavis.edu
// Created: 07/2008


#include <vector>
#include "UniaxialMaterial.h"


class ConfinedConcrete01 : public UniaxialMaterial
{
 public:
	  ConfinedConcrete01(int tag, Vector e, Vector s);

  ConfinedConcrete01(int tag, int secType, int dim, Vector semiLength, 
		     Vector phis, Vector S, Vector fyh, Vector Es0, Vector haRatio, Vector mueps, 
		     Vector As, Vector Is, double rhos, double fpc, double stRatio, double Ec, 
		     int epscuOption, double epscu, double epscuLimit, int nuOption, double nuc, double phiLon, int concrType, 
		     int aggrType, double tol, int maxNumIter, double facToMpa
  );
  ConfinedConcrete01 ();
  ~ConfinedConcrete01();
  
  const char *getClassType(void) const {return "ConfinedConcrete01";};
  
  int setTrialStrain(double strain, double strainRate = 0.0); 
  int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
  double getStrain(void);      
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void) {return 2.0*fpc/epsc0;}
  
  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);        
  
  UniaxialMaterial *getCopy(void);

  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);
  
  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int    setParameter             (const char **argv, int argc, Parameter &param);
  int    updateParameter          (int parameterID, Information &info);
  int    activateParameter        (int parameterID);
  double getStressSensitivity     (int gradNumber, bool conditional);
  int    commitSensitivity        (double strainGradient, int gradNumber, int numGrads);
  // AddingSensitivity:END ///////////////////////////////////////////
  
  int getVariable(const char *variable, Information &);
  
 protected:
  
 private:

  /*** Envelope curve ***/
 Vector eps, sigmac;
  
  /*** Material Properties ***/
  double fpc;    // Compressive strength
  double epsc0;  // Strain at compressive strength
  double fpcu;   // Crushing strength
  double epscu;  // Strain at crushing strength
  double Ec0; // Initial tangent of envelope curve
  
  /*** CONVERGED History Variables ***/
  double CminStrain;   // Smallest previous concrete strain (compression)
  double CunloadSlope; // Unloading (reloading) slope from CminStrain
  double CendStrain;   // Strain at the end of unloading from CminStrain
  
  /*** CONVERGED State Variables ***/
  double Cstrain;
  double Cstress;   
  double Ctangent;	// Don't need Ctangent other than for revert and sendSelf/recvSelf
  // Storing it is better than recomputing it!!!
  
  /*** TRIAL History Variables ***/
  double TminStrain;
  double TunloadSlope;
  double TendStrain;
  
  /*** TRIAL State Variables ***/
  double Tstrain;
  double Tstress;
  double Ttangent; // Not really a state variable, but declared here
  // for convenience
    
  void reload();
  void unload();
  void envelope();
  
  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int parameterID;
  Matrix *SHVs;
  // AddingSensitivity:END ///////////////////////////////////////////

  double convFacToMpa;
};


#endif


