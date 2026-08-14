/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
    10|** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.1 $
// $Date: 1395-12-08 22:25:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/ResidDriftRecorder.cpp,v $

// Written: SAJalali
// Created: 08/06
//
// Description: This file contains the class definition for ResidDriftRecorder.
#include <OPS_Globals.h>
#ifdef _CSS


#include <math.h>

#include <ResidDriftRecorder.h>
#include <Domain.h>
#include <Node.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

ResidDriftRecorder::ResidDriftRecorder()
	 :Recorder(RECORDER_TAGS_ResidDriftRecorder),
	 ndI(0), ndJ(0), theNodes(0), theDofs(0), numDOF(0), perpDirn(0), oneOverL(0),
	 theDomain(0), theOutputHandler(0), data(0),
	 initializationDone(false), numNodes(0), echoTimeFlag(false), dofsFirstFlag(false)
#ifdef _CSS
	 , procDataMethods(0), procGrpNums(0)
#endif // _CSS
{

}


ResidDriftRecorder::ResidDriftRecorder(int ni,
	 int nj,
	 const ID& dofs,
	 int dirn,
	 Domain& theDom,
	 OPS_Stream* theCurrentDataOutputHandler,
	 const ID& procMethods, const ID& procGrpN,
	 bool timeFlag,
	 bool dofsFirst)
	 :Recorder(RECORDER_TAGS_ResidDriftRecorder),
	 ndI(0), ndJ(0), theNodes(0), theDofs(0), numDOF(0), perpDirn(dirn), oneOverL(0),
	 theDomain(&theDom), theOutputHandler(theCurrentDataOutputHandler), data(0),
	 initializationDone(false), numNodes(0), echoTimeFlag(timeFlag), dofsFirstFlag(dofsFirst)
#ifdef _CSS
	 , procDataMethods(procMethods), procGrpNums(procGrpN)
#endif // _CSS
{
	 ndI = new ID(1);
	 ndJ = new ID(1);

	 if (ndI != 0 && ndJ != 0) {
		  (*ndI)(0) = ni;
		  (*ndJ)(0) = nj;
	 }
	 theDofs = new ID(dofs);
	 numDOF = theDofs->Size();
}


ResidDriftRecorder::ResidDriftRecorder(const ID& nI,
	 const ID& nJ,
	 const ID& dofs,
	 int dirn,
	 Domain& theDom,
	 OPS_Stream* theDataOutputHandler,
	 const ID& procMethods, const ID& procGrpN,
	 bool timeFlag,
	 bool dofsFirst)
	 :Recorder(RECORDER_TAGS_ResidDriftRecorder),
	 ndI(0), ndJ(0), theNodes(0), theDofs(0), numDOF(0), perpDirn(dirn), oneOverL(0),
	 theDomain(&theDom), theOutputHandler(theDataOutputHandler), data(0),
	 initializationDone(false), numNodes(0), echoTimeFlag(timeFlag), dofsFirstFlag(dofsFirst)
#ifdef _CSS
	 , procDataMethods(procMethods), procGrpNums(procGrpN)
#endif // _CSS
{
	 ndI = new ID(nI);
	 ndJ = new ID(nJ);
	 theDofs = new ID(dofs);
	 numDOF = theDofs->Size();
}

ResidDriftRecorder::~ResidDriftRecorder()
{
	 //
	 // write the data
	 //
	 if (theOutputHandler != 0 && data != 0) {
		  int size = data->noCols();
		  Vector currentData(size);
		  theOutputHandler->tag("Data"); // Data

		  for (int j = 0; j < size; j++)
				currentData(j) = (*data)(0, j);
		  theOutputHandler->write(currentData);
		  theOutputHandler->endTag(); // Data
		  theOutputHandler->endTag(); // OpenSeesOutput
	 }



	 if (ndI != 0)
		  delete ndI;

	 if (ndJ != 0)
		  delete ndJ;

	 if (theDofs != 0)
		  delete theDofs;

	 if (oneOverL != 0)
		  delete oneOverL;

	 if (theNodes != 0)
		  delete[] theNodes;

	 if (theOutputHandler != 0)
		  delete theOutputHandler;
}

double
ResidDriftRecorder::computeDrift(int pairIndex, int dofIndex) const
{
	 if ((*oneOverL)(pairIndex) == 0.0)
		  return 0.0;
	 Node* nodeI = theNodes[2 * pairIndex];
	 Node* nodeJ = theNodes[2 * pairIndex + 1];
	 const Vector& dispI = nodeI->getTrialDisp();
	 const Vector& dispJ = nodeJ->getTrialDisp();
	 int d = (*theDofs)(dofIndex);
	 if (d < 0 || d >= dispI.Size() || d >= dispJ.Size())
		  return 0.0;
	 return (dispJ(d) - dispI(d)) * (*oneOverL)(pairIndex);
}

int
ResidDriftRecorder::record(int commitTag, double timeStamp)
{

	 if (theDomain == 0 || ndI == 0 || ndJ == 0) {
		  return 0;
	 }

	 if (initializationDone != true)
		  if (this->initialize() != 0) {
				opserr << "ResidDriftRecorder::record() - failed in initialize()\n";
				return -1;
		  }

	 if (numNodes == 0 || data == 0)
		  return 0;
	 int timeOffset = 0;
	 if (echoTimeFlag)
	 {
		  timeOffset = 1;
		  (*data)(0, 0) = timeStamp;
	 }
	 if (procDataMethods.Size() != 0)
	 {
		  // Materialize full raw layout, then reduce once
		  int nDof = (numDOF > 0) ? numDOF : 1;
		  int nRaw = numNodes * nDof;
		  double* raw = new double[nRaw];
		  if (dofsFirstFlag) {
				for (int i = 0; i < numNodes; i++)
					 for (int j = 0; j < nDof; j++)
						  raw[i * nDof + j] = this->computeDrift(i, j);
		  }
		  else {
				for (int j = 0; j < nDof; j++)
					 for (int i = 0; i < numNodes; i++)
						  raw[j * numNodes + i] = this->computeDrift(i, j);
		  }
		  int nOut = Recorder::applyProcDataChain(procDataMethods, procGrpNums, raw, nRaw, false);
		  for (int i = 0; i < nOut; i++)
				(*data)(0, timeOffset + i) = raw[i];
		  delete[] raw;
	 }
	 else {
		  int nDof = (numDOF > 0) ? numDOF : 1;
		  if (dofsFirstFlag) {
				for (int i = 0; i < numNodes; i++)
					 for (int j = 0; j < nDof; j++)
						  (*data)(0, timeOffset + i * nDof + j) = this->computeDrift(i, j);
		  }
		  else {
				for (int j = 0; j < nDof; j++)
					 for (int i = 0; i < numNodes; i++)
						  (*data)(0, timeOffset + j * numNodes + i) = this->computeDrift(i, j);
		  }
	 }

	 return 0;
}

int
ResidDriftRecorder::restart(void)
{
	 data->Zero();
	 return 0;
}

int
ResidDriftRecorder::setDomain(Domain& theDom)
{
	 theDomain = &theDom;
	 initializationDone = false;
	 return 0;
}

int
ResidDriftRecorder::sendSelf(int commitTag, Channel& theChannel)
{
	 static ID idData(7);
	 idData.Zero();
	 if (ndI != 0 && ndI->Size() != 0)
		  idData(0) = ndI->Size();
	 if (ndJ != 0 && ndJ->Size() != 0)
		  idData(1) = ndJ->Size();
	 idData(2) = numDOF;
	 idData(3) = perpDirn;
	 if (theOutputHandler != 0) {
		  idData(4) = theOutputHandler->getClassTag();
	 }
	 else
		 idData(4) = 0;
	 if (echoTimeFlag == true)
		  idData(5) = 0;
	 else
		  idData(5) = 1;
	 idData(6) = dofsFirstFlag ? 1 : 0;

	 if (theChannel.sendID(0, commitTag, idData) < 0) {
		  opserr << "ResidDriftRecorder::sendSelf() - failed to send idData\n";
		  return -1;
	 }

	 if (ndI != 0)
		  if (theChannel.sendID(0, commitTag, *ndI) < 0) {
				opserr << "ResidDriftRecorder::sendSelf() - failed to send dof id's\n";
				return -1;
		  }

	 if (ndJ != 0)
		  if (theChannel.sendID(0, commitTag, *ndJ) < 0) {
				opserr << "ResidDriftRecorder::sendSelf() - failed to send dof id's\n";
				return -1;
		  }

	 if (theDofs != 0)
		  if (theChannel.sendID(0, commitTag, *theDofs) < 0) {
				opserr << "ResidDriftRecorder::sendSelf() - failed to send theDofs\n";
				return -1;
		  }

	 if (theOutputHandler != 0)
		  if (theOutputHandler->sendSelf(commitTag, theChannel) < 0) {
				opserr << "ResidDriftRecorder::sendSelf() - failed to send the DataOutputHandler\n";
				return -1;
		  }

	 return 0;
}

int
ResidDriftRecorder::recvSelf(int commitTag, Channel& theChannel,
	 FEM_ObjectBroker& theBroker)
{
	 static ID idData(7);
	 if (theChannel.recvID(0, commitTag, idData) < 0) {
		  opserr << "ResidDriftRecorder::sendSelf() - failed to send idData\n";
		  return -1;
	 }

	 if (idData(0) != 0) {
		  ndI = new ID(idData(0));
		  if (ndI == 0) {
				opserr << "ResidDriftRecorder::sendSelf() - out of memory\n";
				return -1;
		  }
		  if (theChannel.recvID(0, commitTag, *ndI) < 0) {
				opserr << "ResidDriftRecorder::sendSelf() - failed to recv dof id's\n";
				return -1;
		  }
	 }

	 if (idData(1) != 0) {

		  ndJ = new ID(idData(1));
		  if (ndJ == 0) {
				opserr << "ResidDriftRecorder::sendSelf() - out of memory\n";
				return -1;
		  }
		  if (theChannel.recvID(0, commitTag, *ndJ) < 0) {
				opserr << "ResidDriftRecorder::sendSelf() - failed to recv dof id's\n";
				return -1;
		  }
	 }

	 numDOF = idData(2);
	 perpDirn = idData(3);
	 dofsFirstFlag = (idData(6) == 1);

	 if (numDOF != 0) {
		  theDofs = new ID(numDOF);
		  if (theChannel.recvID(0, commitTag, *theDofs) < 0) {
				opserr << "ResidDriftRecorder::recvSelf() - failed to recv theDofs\n";
				return -1;
		  }
	 }

	 if (idData(5) == 0)
		  echoTimeFlag = true;
	 else
		  echoTimeFlag = false;

	 if (theOutputHandler != 0)
		  delete theOutputHandler;
	 if (idData(4) != 0)
	 {
		 theOutputHandler = theBroker.getPtrNewStream(idData(4));
		 if (theOutputHandler == 0) {
			 opserr << "ResidDriftRecorder::sendSelf() - failed to get a data output handler\n";
			 return -1;
		 }

		 if (theOutputHandler->recvSelf(commitTag, theChannel, theBroker) < 0) {
			 delete theOutputHandler;
			 theOutputHandler = 0;
		 }
	 }

	 return 0;
}


int
ResidDriftRecorder::initialize(void)
{
	 if (theOutputHandler != 0)
		  theOutputHandler->tag("OpenSeesOutput");

	 initializationDone = true; // still might fail but don't want back in again

	 //
	 // clean up old memory
	 //

	 if (theNodes != 0) {
		  delete[] theNodes;
		  theNodes = 0;
	 }
	 if (oneOverL != 0) {
		  delete oneOverL;
		  oneOverL = 0;
	 }

	 //
	 // check valid node ID's
	 //

	 if (ndI == 0 || ndJ == 0) {
		  opserr << "ResidDriftRecorder::initialize() - no nodal id's set\n";
		  return -1;
	 }

	 int ndIsize = ndI->Size();
	 int ndJsize = ndJ->Size();

	 if (ndIsize == 0) {
		  opserr << "ResidDriftRecorder::initialize() - no nodal id's set\n";
		  return -1;
	 }

	 if (ndIsize != ndJsize) {
		  opserr << "ResidDriftRecorder::initialize() - error node arrays differ in size\n";
		  return -2;
	 }

	 //
	 // lets loop through & determine number of valid nodes
	 //


	 numNodes = 0;

	 for (int i = 0; i < ndIsize; i++) {
		  int ni = (*ndI)(i);
		  int nj = (*ndJ)(i);

		  Node* nodeI = theDomain->getNode(ni);
		  Node* nodeJ = theDomain->getNode(nj);

		  if (nodeI != 0 && nodeJ != 0) {
				const Vector& crdI = nodeI->getCrds();
				const Vector& crdJ = nodeJ->getCrds();

				if (crdI.Size() > perpDirn && crdJ.Size() > perpDirn)
					 if (crdI(perpDirn) != crdJ(perpDirn))
						  numNodes++;
		  }
	 }

	 if (numNodes == 0) {
		  opserr << "ResidDriftRecorder::initialize() - no valid nodes or perpendicular direction\n";
		  return 0;
	 }

	 //
	 // allocate memory
	 //
	 int timeOffset = 0;
	 if (echoTimeFlag == true)
		  timeOffset = 1;

	 int nDof = (numDOF > 0) ? numDOF : 1;
	 int nRaw = numNodes * nDof;
	 int nProcOuts = (procDataMethods.Size() == 0) ? nRaw
		  : Recorder::getFinalProcOuts(nRaw, procDataMethods, procGrpNums);

	 data = new Matrix(1, nProcOuts + timeOffset);
	 data->Zero();
	 theNodes = new Node * [2 * numNodes];
	 oneOverL = new Vector(numNodes);
	 if (theNodes == 0 || oneOverL == 0) {
		  opserr << "ResidDriftRecorder::initialize() - out of memory\n";
		  return -3;
	 }

	 //
	 // set node pointers and determine one over L
	 //

	 int counter = 0;
	 int counterI = 0;
	 int counterJ = 1;
	 for (int j = 0; j < ndIsize; j++) {
		  int ni = (*ndI)(j);
		  int nj = (*ndJ)(j);

		  Node* nodeI = theDomain->getNode(ni);
		  Node* nodeJ = theDomain->getNode(nj);

		  if (nodeI != 0 && nodeJ != 0) {
				const Vector& crdI = nodeI->getCrds();
				const Vector& crdJ = nodeJ->getCrds();

				if (crdI.Size() > perpDirn && crdJ.Size() > perpDirn)
					 if (crdI(perpDirn) != crdJ(perpDirn)) {

						  if (theOutputHandler != 0)
						  {
								theOutputHandler->tag("DriftOutput");
								theOutputHandler->attr("node1", ni);
								theOutputHandler->attr("node2", ni);
								theOutputHandler->attr("perpDirn", perpDirn);
								theOutputHandler->attr("lengthPerpDirn", fabs(crdJ(perpDirn) - crdI(perpDirn)));

								if (echoTimeFlag == true) {
									 theOutputHandler->tag("TimeOutput");
									 theOutputHandler->tag("ResponseType", "time");
									 theOutputHandler->endTag(); // TimeOutput
								}

								theOutputHandler->tag("ResponseType", "drift");
								theOutputHandler->endTag(); // DriftOutput
						  }

						  (*oneOverL)(counter) = 1.0 / fabs(crdJ(perpDirn) - crdI(perpDirn));
						  theNodes[counterI] = nodeI;
						  theNodes[counterJ] = nodeJ;
						  counterI += 2;
						  counterJ += 2;
						  counter++;
					 }
		  }
	 }


	 //
	 // mark as having been done & return
	 //

	 return 0;
}

#endif // _CSS
