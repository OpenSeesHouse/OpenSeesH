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

// $Revision: 1.8 $
// $Date: 2006-01-10 18:13:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/NormDispAndUnbalance.cpp,v $


#include <NormDispAndUnbalance.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>
#include <elementAPI.h>

void* OPS_NormDispAndUnbalance()
{
    double data[2] = {0,0};
    int idata[4] = {0,0,2,-1};

    if (OPS_GetNumRemainingInputArgs() >= 2) {
	int numData = 2;
	if (OPS_GetDoubleInput(&numData, data) < 0) {
	    opserr << "WARNING NormDispAndUnbalance failed to get tolerance\n";
	    return 0;
	}
	numData = 2;
	if (OPS_GetIntInput(&numData, idata) < 0) {
	    opserr << "WARNING NormDispAndUnbalance failed to get max iter and print flag\n";
	    return 0;
	}	
    }
    else {
      opserr << "ERROR - NormDispAndUnbalance insufficient arguments" << endln;
      return 0;
    }

    if (OPS_GetNumRemainingInputArgs() >= 2) {
	int numData = 2;
	if (OPS_GetIntInput(&numData, idata) < 0) {
	    opserr << "WARNING NormDispAndUnbalance failed to get int values\n";
	    return 0;
	}
    }

    return new NormDispAndUnbalance(data[0],data[1],idata[0],idata[1],idata[2],idata[3]);
}

NormDispAndUnbalance::NormDispAndUnbalance()	    	
  : ConvergenceTest(CONVERGENCE_TEST_NormDispAndUnbalance),
    theSOE(0), tolDisp(0), tolUnbalance(0),
    maxNumIter(0), currentIter(0), printFlag(0), 
    norms(25), nType(2), maxIncr(0), numIncr(0)
{
    
}


NormDispAndUnbalance::NormDispAndUnbalance(double theTolDisp, double theTolUnbalance, int maxIter, int printIt, int normType, int maxincr)
  : ConvergenceTest(CONVERGENCE_TEST_NormDispAndUnbalance),
    theSOE(0), tolDisp(theTolDisp), tolUnbalance(theTolUnbalance),
    maxNumIter(maxIter), currentIter(0), printFlag(printIt),
    norms(2*maxIter), nType(normType), maxIncr(maxincr), numIncr(0)
{
    if(maxIncr < 0) maxIncr = maxNumIter;
}


NormDispAndUnbalance::~NormDispAndUnbalance()
{
  
}


ConvergenceTest* NormDispAndUnbalance::getCopy(int iterations)
{
  NormDispAndUnbalance *theCopy ;
  theCopy = new NormDispAndUnbalance(this->tolDisp, 
				     this->tolUnbalance,
				     iterations, 
				     this->printFlag, 
				     this->nType,
                                     this->maxIncr) ;
  
  theCopy->theSOE = this->theSOE ;
  
  return theCopy ;
}


void NormDispAndUnbalance::setTolerance(double newTolDisp)
{
    tolDisp = newTolDisp;
}


int NormDispAndUnbalance::setEquiSolnAlgo(EquiSolnAlgo &theAlgo)
{
    theSOE = theAlgo.getLinearSOEptr();
    return 0;
}


int NormDispAndUnbalance::test(void)
{
    // check to ensure the SOE has been set - this should not happen if the 
    // return from start() is checked
    if (theSOE == 0) {
		 opserr << "WARNING: NormDispAndUnbalance::test() - no soe set.\n";	
        return -2;
	}
    
    // check to ensure the algo does invoke start() - this is needed otherwise
    // may never get convergence later on in analysis!
    if (currentIter == 0) {
        opserr << "WARNING: NormDispAndUnbalance::test() - start() was never invoked.\n";	
        return -2;
    }
    
    // get the X vector & determine it's norm & save the value in norms vector
    const Vector &x = theSOE->getX();
    double normX = x.pNorm(nType);

    const Vector &b = theSOE->getB();
    double normB = b.pNorm(nType);

    if((currentIter>1 && norms(currentIter-2)<normX) || (currentIter>1 && norms(maxNumIter+currentIter-2)<normB)) {
        numIncr++;
    }

    if (currentIter <= maxNumIter) {
        norms(currentIter-1) = normX;
        norms(maxNumIter+currentIter-1) = normB;
    }
    
    // print the data if required
    if (printFlag == 1) {
        opserr << "NormDispAndUnbalance::test() - iteration: " << currentIter;
        opserr << " current NormX: " << normX;
        opserr << ", NormB: " << normB  << ", NormIncr: " << numIncr << "\n";
    } 
    if (printFlag == 4) {
        opserr << "NormDispAndUnbalance::test() - iteration: " << currentIter;
        opserr << " current NormX: " << normX;
        opserr << ", NormB: " << normB  << ", NormIncr: " << numIncr << "\n";
        opserr << "\tdeltaX: " << x << "\tdeltaR: " << theSOE->getB();
    } 
    
    //
    // check if the algorithm converged
    //
    
    // if converged - print & return ok
    if (normX <= tolDisp && normB <= tolUnbalance) { 
        
        // do some printing first
        if (printFlag != 0) {
            if (printFlag == 1 || printFlag == 4) 
                opserr << endln;
            else if (printFlag == 2 || printFlag == 6) {
                opserr << "NormDispAndUnbalance::test() - iteration: " << currentIter;
		opserr << " current NormX: " << normX;
		opserr << ", NormB: " << normB  << ", NormIncr: " << numIncr << "\n";
            }
        }
        
        // return the number of times test has been called
        return currentIter;
    }
    
    // algo failed to converged after specified number of iterations - but RETURN OK
    else if ((printFlag == 5 || printFlag == 6) && (currentIter >= maxNumIter || numIncr > maxIncr)) {
        opserr << "WARNING: NormDispAndUnbalance::test() - failed to converge but going on - ";
	opserr << " current NormX: " << normX;
	opserr << ", NormB: " << normB  << ", NormIncr: " << numIncr << "\n";
        return currentIter;
    }
    
    // algo failed to converged after specified number of iterations - return FAILURE -2
    else if (currentIter >= maxNumIter || numIncr > maxIncr) { // fails to converge
        opserr << "WARNING: NormDispAndUnbalance::test() - failed to converge \n";
        opserr << "after: " << currentIter << " iterations\n";	
        currentIter++;    
        return -2;
    } 
    
    // algorithm not yet converged - increment counter and return -1
    else {
        currentIter++;    
        return -1;
    }
}


int NormDispAndUnbalance::start(void)
{
    if (theSOE == 0) {
        opserr << "WARNING: NormDispAndUnbalance::test() - no SOE returning true\n";
        return -1;
    }
    
    // set iteration count = 1
    norms.Zero();
    currentIter = 1;
    numIncr = 0;
    return 0;
}


int NormDispAndUnbalance::getNumTests()
{
    return currentIter;
}


int NormDispAndUnbalance::getMaxNumTests(void)
{
    return maxNumIter;
}


double NormDispAndUnbalance::getRatioNumToMax(void)
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector& NormDispAndUnbalance::getNorms() 
{
    return norms;
}


int NormDispAndUnbalance::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  Vector x(6);
  x(0) = tolDisp;
  x(4) = tolUnbalance;
  x(1) = maxNumIter;
  x(2) = printFlag;
  x(3) = nType;
  x(5) = maxIncr;
  res = theChannel.sendVector(this->getDbTag(), cTag, x);
  if (res < 0) 
    opserr << "NormDispAndUnbalance::sendSelf() - failed to send data\n";
    
  return res;
}

int 
NormDispAndUnbalance::recvSelf(int cTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  Vector x(6);
  res = theChannel.recvVector(this->getDbTag(), cTag, x);    


  if (res < 0) {
    opserr << "NormDispAndUnbalance::sendSelf() - failed to send data\n";
    tolDisp = 1.0e-8;
    maxNumIter = 25;
    printFlag = 0;
    nType = 2;
    maxIncr = 3;
    norms.resize(maxNumIter);
  } else {
    tolDisp = x(0);
    tolUnbalance = x(4);
    maxNumIter = (int)x(1);
    printFlag = (int)x(2);
    nType = (int)x(3);
    maxIncr = (int)x(5);
    norms.resize(maxNumIter);
  } 
  return res;
}


