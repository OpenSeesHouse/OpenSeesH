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
                                                                        
// $Revision: 1.4 $
// $Date: 2010-04-23 22:47:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/Recorder.cpp,v $
                                                                        

// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class implementation for Recorder.
//
// What: "@(#) Recorder.cpp, revA"

#include <Recorder.h>
#include <OPS_Globals.h>
#include <ID.h>
#include <math.h>
#include <stdlib.h>

int Recorder::lastRecorderTag(0);
const double Recorder::relDeltaTTol = 1e-5;

#ifdef _CSS
int
Recorder::getNumProcOuts(int nVals, int procGrpNum)
{
	if (procGrpNum == -1)
		return 1;
	int nProcOuts = nVals / procGrpNum;
	if (nProcOuts * procGrpNum < nVals)
		nProcOuts++;
	return nProcOuts;
}

int
Recorder::getFinalProcOuts(int nVals, const ID& methods, const ID& grpNums, bool warn)
{
	if (methods.Size() == 0)
		return nVals;
	int n = nVals;
	for (int s = 0; s < methods.Size(); s++) {
		n = getNumProcOuts(n, grpNums(s));
		if (n == 1 && s + 1 < methods.Size()) {
			if (warn) {
				opserr << "WARNING recorder: process stage " << (s + 1)
					<< " (grpNum=" << grpNums(s)
					<< ") reduced output to 1 column; subsequent process operations are skipped.\n";
			}
			return 1;
		}
	}
	return n;
}

int
Recorder::applyProcDataStage(int method, int grpNum, const double* inVals, int nIn, double* outVals)
{
	int nOut = getNumProcOuts(nIn, grpNum);
	for (int i = 0; i < nOut; i++)
		outVals[i] = 0.0;

	int iGrpN = 0;
	int nextGrpN = grpNum;
	double* val = &outVals[iGrpN];
	for (int i = 0; i < nIn; i++) {
		double val1 = inVals[i];
		bool newGrp = false;
		if (grpNum != -1 && i == nextGrpN) {
			iGrpN++;
			nextGrpN += grpNum;
			val = &outVals[iGrpN];
			newGrp = true;
		}
		if ((i == 0 || newGrp) && method != 1 && method != 6)
			*val = fabs(val1);
		if (method == 1)
			*val += val1;
		else if (method == 2 && val1 > *val)
			*val = val1;
		else if (method == 3 && val1 < *val)
			*val = val1;
		else if (method == 4 && fabs(val1) > *val)
			*val = fabs(val1);
		else if (method == 5 && fabs(val1) < *val)
			*val = fabs(val1);
		else if (method == 6)
			*val += val1 * val1;
	}
	if (method == 6) {
		for (int i = 0; i < nOut; i++)
			outVals[i] = sqrt(outVals[i]);
	}
	return nOut;
}

int
Recorder::applyProcDataChain(const ID& methods, const ID& grpNums, double* vals, int nVals, bool warn)
{
	if (methods.Size() == 0)
		return nVals;

	double* tmp = new double[nVals];
	if (tmp == 0)
		return nVals;

	int n = nVals;
	for (int s = 0; s < methods.Size(); s++) {
		int nOut = applyProcDataStage(methods(s), grpNums(s), vals, n, tmp);
		for (int i = 0; i < nOut; i++)
			vals[i] = tmp[i];
		n = nOut;
		if (n == 1 && s + 1 < methods.Size()) {
			if (warn) {
				opserr << "WARNING recorder: process stage " << (s + 1)
					<< " (grpNum=" << grpNums(s)
					<< ") reduced output to 1 column; subsequent process operations are skipped.\n";
			}
			break;
		}
	}
	delete[] tmp;
	return n;
}
#endif // _CSS

Recorder::Recorder(int classTag)
  :MovableObject(classTag), TaggedObject(lastRecorderTag)
{
  lastRecorderTag++;
}

Recorder::~Recorder() 
{

}

int 
Recorder::restart(void)
{
  return 0;
}

int 
Recorder::flush(void)
{
  return 0;
}

int 
Recorder::domainChanged(void)
{
  return 0;
}

int 
Recorder::setDomain(Domain &theDomain)
{
  return 0;
}

int 
Recorder::sendSelf(int commitTag, Channel &theChannel)
{
  opserr << "Recorder::sendSelf() - not yet implemented\n";
  return 0;
}

int 
Recorder::recvSelf(int commitTag, Channel &theChannel, 
		   FEM_ObjectBroker &theBroker)
{
  opserr << "Recorder::recvSelf() - not yet implemented\n";
  return 0;
}

void Recorder::resetLastTag()
{
  lastRecorderTag = 0;
}

void
Recorder::Print(OPS_Stream &s, int flag)
{
  return;
}
