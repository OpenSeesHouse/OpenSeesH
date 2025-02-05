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

// $Revision: 1.1 $
// $Date: 2006-08-17 22:25:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/EnvelopeDriftRecorder.cpp,v $

// Written: fmk
// Created: 08/06
//
// Description: This file contains the class definition for EnvelopeDriftRecorder.

#include <BinaryFileStream.h>
#include <Channel.h>
#include <DataFileStream.h>
#include <Domain.h>
#include <EnvelopeDriftRecorder.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>
#include <Matrix.h>
#include <Node.h>
#include <StandardStream.h>
#include <Vector.h>
#include <XmlFileStream.h>
#include <elementAPI.h>
#include <math.h>
#include <string.h>

enum outputMode {
  STANDARD_STREAM,
  DATA_STREAM,
  XML_STREAM,
  DATABASE_STREAM,
  BINARY_STREAM,
  DATA_STREAM_CSV,
  TCP_STREAM,
  DATA_STREAM_ADD
};

EnvelopeDriftRecorder::EnvelopeDriftRecorder()
	 :Recorder(RECORDER_TAGS_EnvelopeDriftRecorder),
	 ndI(0), ndJ(0), dof(0), perpDirn(0), oneOverL(0), currentData(0),
	 theDomain(0), theOutputHandler(0),
	 initializationDone(false), numNodes(0), echoTimeFlag(false)
{

}


EnvelopeDriftRecorder::EnvelopeDriftRecorder(int ni,
	 int nj,
	 int df,
	 int dirn,
	 Domain& theDom,
	 OPS_Stream* theCurrentDataOutputHandler,
	 int procMethod, int procGrpN,
	 bool timeFlag)
	 :Recorder(RECORDER_TAGS_EnvelopeDriftRecorder),
	 ndI(0), ndJ(0), theNodes(0), dof(df), perpDirn(dirn), oneOverL(0), currentData(0),
	 theDomain(&theDom), theOutputHandler(theCurrentDataOutputHandler),
	 initializationDone(false), numNodes(0), echoTimeFlag(timeFlag)
#ifdef _CSS
	 , procDataMethod(procMethod), procGrpNum(procGrpN), Modified(0)
#endif // _CSS
{
	 ndI = new ID(1);
	 ndJ = new ID(1);

	 if (ndI != 0 && ndJ != 0) {
		  (*ndI)(0) = ni;
		  (*ndJ)(0) = nj;
	 }
}


EnvelopeDriftRecorder::EnvelopeDriftRecorder(const ID& nI,
	 const ID& nJ,
	 int df,
	 int dirn,
	 Domain& theDom,
	 OPS_Stream* theDataOutputHandler,
	 int procMethod, int procGrpN,
	 bool timeFlag)
	 :Recorder(RECORDER_TAGS_EnvelopeDriftRecorder),
	 ndI(0), ndJ(0), theNodes(0), dof(df), perpDirn(dirn), oneOverL(0), currentData(0),
	 theDomain(&theDom), theOutputHandler(theDataOutputHandler),
	 initializationDone(false), numNodes(0), echoTimeFlag(timeFlag)
#ifdef _CSS
	 , procDataMethod(procMethod), procGrpNum(procGrpN), Modified(0)
#endif // _CSS
{
	 ndI = new ID(nI);
	 ndJ = new ID(nJ);
}

EnvelopeDriftRecorder::~EnvelopeDriftRecorder()
{
	 //
	 // write the data
	 //
	 if (theOutputHandler != 0 && currentData != 0) {

		  theOutputHandler->tag("Data"); // Data
		  for (int i = 0; i < 3; i++) {
				int size = currentData->Size();
				for (int j = 0; j < size; j++)
					 (*currentData)(j) = (*data)(i, j);
				theOutputHandler->write(*currentData);
		  }
		  theOutputHandler->endTag(); // Data
		  theOutputHandler->endTag(); // OpenSeesOutput
	 }



	 if (ndI != 0)
		  delete ndI;

	 if (ndJ != 0)
		  delete ndJ;

	 if (oneOverL != 0)
		  delete oneOverL;

	 if (currentData != 0)
		  delete currentData;

	 if (theNodes != 0)
		  delete[] theNodes;

	 if (theOutputHandler != 0)
		  delete theOutputHandler;
}

int
EnvelopeDriftRecorder::record(int commitTag, double timeStamp)
{

	 if (theDomain == 0 || ndI == 0 || ndJ == 0) {
		  return 0;
	 }
	 if (initializationDone != true)
		  if (this->initialize() != 0) {
				opserr << "EnvelopeDriftRecorder::record() - failed in initialize()\n";
				return -1;
		  }

	 if (numNodes == 0 || currentData == 0)
		  return 0;

#ifdef _CSS
	 Modified = 0;
	 if (procDataMethod)
	 {
		  int nProcOuts;
		  int nVals = numNodes;
		  if (procGrpNum == -1)
				if (procDataMethod != 0)
					 nProcOuts = 1;
				else
					 nProcOuts = nVals;
		  else {
				nProcOuts = nVals / procGrpNum;
				if (nProcOuts * procGrpNum < nVals)
					 nProcOuts++;
		  }
		  double* vals = 0, * val, val1 = 0;
		  vals = new double[nProcOuts];
		  for (int i = 0; i < nProcOuts; i++)
				vals[i] = 0;
		  int iGrpN = 0;
		  int nextGrpN = procGrpNum;
		  val = &vals[iGrpN];
		  int loc = 0;
		  for (int i = 0; i < numNodes; i++) {
				Node* nodeI = theNodes[2 * i];
				Node* nodeJ = theNodes[2 * i + 1];

				if ((*oneOverL)(i) != 0.0) {
					 const Vector& dispI = nodeI->getTrialDisp();
					 const Vector& dispJ = nodeJ->getTrialDisp();

					 double dx = dispJ(dof) - dispI(dof);

					 val1 = dx * (*oneOverL)(i);

				}
				else
					 val1 = 0.0;
				if (procGrpNum != -1 && i == nextGrpN)
				{
					 iGrpN++;
					 nextGrpN += procGrpNum;
					 val = &vals[iGrpN];
				}

				if (i == 0 && procDataMethod != 1)
					 *val = fabs(val1);
				if (procDataMethod == 1)
					 *val += val1;
				else if (procDataMethod == 2 && val1 > *val)
					 *val = val1;
				else if (procDataMethod == 3 && val1 < *val)
					 *val = val1;
				else if (procDataMethod == 4 && fabs(val1) > *val)
					 *val = fabs(val1);
				else if (procDataMethod == 5 && fabs(val1) < *val)
					 *val = fabs(val1);
		  }
		  for (int i = 0; i < nProcOuts; i++)
		  {
				val = &vals[i];
				(*currentData)(loc++) = *val;
		  }
		  delete[] vals;
	 }
	 else
#endif // _CSS
		  for (int i = 0; i < numNodes; i++) {
				Node* nodeI = theNodes[2 * i];
				Node* nodeJ = theNodes[2 * i + 1];

				if ((*oneOverL)(i) != 0.0) {
					 const Vector& dispI = nodeI->getTrialDisp();
					 const Vector& dispJ = nodeJ->getTrialDisp();

					 double dx = dispJ(dof) - dispI(dof);

					 (*currentData)(i) = dx * (*oneOverL)(i);

				}
				else
					 (*currentData)(i) = 0.0;
		  }

	 // check if currentData modifies the saved data
	 int sizeData = currentData->Size();
	 if (echoTimeFlag == false) {

		  bool writeIt = false;
		  if (first == true) {
				for (int i = 0; i < sizeData; i++) {
					 (*data)(0, i) = (*currentData)(i);
					 (*data)(1, i) = (*currentData)(i);
					 (*data)(2, i) = fabs((*currentData)(i));
					 first = false;
					 writeIt = true;
#ifdef _CSS
					 Modified = 1;
#endif // _CSS

				}
		  }
		  else {
				for (int i = 0; i < sizeData; i++) {
					 double value = (*currentData)(i);
					 if ((*data)(0, i) > value) {
						  (*data)(0, i) = value;
						  double absValue = fabs(value);
						  if ((*data)(2, i) < absValue)
#ifdef _CSS
						  {
								Modified = 1;
								(*data)(2, i) = absValue;
						  }
#else
								(*data)(2, i) = absValue;
#endif // _CSS
						  writeIt = true;
					 }
					 else if ((*data)(1, i) < value) {
						  (*data)(1, i) = value;
						  double absValue = fabs(value);
						  if ((*data)(2, i) < absValue)
#ifdef _CSS
						  {
								Modified = 1;
								(*data)(2, i) = absValue;
						  }
#else
								(*data)(2, i) = absValue;
#endif // _CSS
						  writeIt = true;
					 }
				}
		  }
	 }
	 else {
		  sizeData /= 2;
		  bool writeIt = false;
		  if (first == true) {
				for (int i = 0; i < sizeData; i++) {

					 (*data)(0, i * 2) = timeStamp;
					 (*data)(1, i * 2) = timeStamp;
					 (*data)(2, i * 2) = timeStamp;
					 (*data)(0, i * 2 + 1) = (*currentData)(i);
					 (*data)(1, i * 2 + 1) = (*currentData)(i);
					 (*data)(2, i * 2 + 1) = fabs((*currentData)(i));
					 first = false;
					 writeIt = true;
#ifdef _CSS
					 Modified = 1;
#endif // _CSS
				}
		  }
		  else {
				for (int i = 0; i < sizeData; i++) {
					 double value = (*currentData)(i);
					 if ((*data)(0, 2 * i + 1) > value) {
						  (*data)(0, i * 2) = timeStamp;
						  (*data)(0, i * 2 + 1) = value;
						  double absValue = fabs(value);
						  if ((*data)(2, i * 2 + 1) < absValue) {
								(*data)(2, i * 2 + 1) = absValue;
								(*data)(2, i * 2) = timeStamp;
#ifdef _CSS
								Modified = 1;
#endif // _CSS
						  }
						  writeIt = true;
					 }
					 else if ((*data)(1, i * 2 + 1) < value) {
						  (*data)(1, i * 2) = timeStamp;
						  (*data)(1, i * 2 + 1) = value;
						  double absValue = fabs(value);
						  if ((*data)(2, i * 2 + 1) < absValue) {
								(*data)(2, i * 2) = timeStamp;
								(*data)(2, i * 2 + 1) = absValue;
#ifdef _CSS
								Modified = 1;
#endif // _CSS
						  }
						  writeIt = true;
					 }
				}
		  }
	 }
	 return 0;
}

int
EnvelopeDriftRecorder::restart(void)
{
	 data->Zero();
	 first = true;
	 return 0;
}

int
EnvelopeDriftRecorder::setDomain(Domain& theDom)
{
	 theDomain = &theDom;
	 initializationDone = false;
	 return 0;
}

int
EnvelopeDriftRecorder::sendSelf(int commitTag, Channel& theChannel)
{
	 static ID idData(6);
	 idData.Zero();
	 if (ndI != 0 && ndI->Size() != 0)
		  idData(0) = ndI->Size();
	 if (ndJ != 0 && ndJ->Size() != 0)
		  idData(1) = ndJ->Size();
	 idData(2) = dof;
	 idData(3) = perpDirn;
	 if (theOutputHandler != 0) {
		  idData(4) = theOutputHandler->getClassTag();
	 }
	 if (echoTimeFlag == true)
		  idData(5) = 0;
	 else
		  idData(5) = 1;

	 if (theChannel.sendID(0, commitTag, idData) < 0) {
		  opserr << "EnvelopeDriftRecorder::sendSelf() - failed to send idData\n";
		  return -1;
	 }

	 if (ndI != 0)
		  if (theChannel.sendID(0, commitTag, *ndI) < 0) {
				opserr << "EnvelopeDriftRecorder::sendSelf() - failed to send dof id's\n";
				return -1;
		  }

	 if (ndJ != 0)
		  if (theChannel.sendID(0, commitTag, *ndJ) < 0) {
				opserr << "EnvelopeDriftRecorder::sendSelf() - failed to send dof id's\n";
				return -1;
		  }


	 if (theOutputHandler != 0)
		  if (theOutputHandler->sendSelf(commitTag, theChannel) < 0) {
				opserr << "EnvelopeDriftRecorder::sendSelf() - failed to send the DataOutputHandler\n";
				return -1;
		  }

	 return 0;
}

int
EnvelopeDriftRecorder::recvSelf(int commitTag, Channel& theChannel,
	 FEM_ObjectBroker& theBroker)
{
  static ID idData(6); 
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "EnvelopeDriftRecorder::sendSelf() - failed to send idData\n";
    return -1;
  }
  
  if (idData(0) != 0) {
    ndI = new ID(idData(0));
    if (ndI == 0) {
      opserr << "EnvelopeDriftRecorder::sendSelf() - out of memory\n";
      return -1;
    }
    if (theChannel.recvID(0, commitTag, *ndI) < 0) {
      opserr << "EnvelopeDriftRecorder::sendSelf() - failed to recv dof id's\n";
      return -1;
    } 
  }

	 if (idData(1) != 0) {

		  ndJ = new ID(idData(1));
		  if (ndJ == 0) {
				opserr << "EnvelopeDriftRecorder::sendSelf() - out of memory\n";
				return -1;
		  }
		  if (theChannel.recvID(0, commitTag, *ndJ) < 0) {
				opserr << "EnvelopeDriftRecorder::sendSelf() - failed to recv dof id's\n";
				return -1;
		  }
	 }

	 dof = idData(2);
	 perpDirn = idData(3);

	 if (idData(5) == 0)
		  echoTimeFlag = true;
	 else
		  echoTimeFlag = false;

	 if (theOutputHandler != 0)
		  delete theOutputHandler;

	 theOutputHandler = theBroker.getPtrNewStream(idData(4));
	 if (theOutputHandler == 0) {
		  opserr << "EnvelopeDriftRecorder::sendSelf() - failed to get a data output handler\n";
		  return -1;
	 }

	 if (theOutputHandler->recvSelf(commitTag, theChannel, theBroker) < 0) {
		  delete theOutputHandler;
		  theOutputHandler = 0;
	 }

	 return 0;
}


int
EnvelopeDriftRecorder::initialize(void)
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
	 if (currentData != 0) {
		  delete currentData;
		  currentData = 0;
	 }
	 if (oneOverL != 0) {
		  delete oneOverL;
		  oneOverL = 0;
	 }

	 //
	 // check valid node ID's
	 //

	 if (ndI == 0 || ndJ == 0) {
		  opserr << "EnvelopeDriftRecorder::initialize() - no nodal id's set\n";
		  return -1;
	 }

	 int ndIsize = ndI->Size();
	 int ndJsize = ndJ->Size();

	 if (ndIsize == 0) {
		  opserr << "EnvelopeDriftRecorder::initialize() - no nodal id's set\n";
		  return -1;
	 }

	 if (ndIsize != ndJsize) {
		  opserr << "EnvelopeDriftRecorder::initialize() - error node arrays differ in size\n";
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
		  opserr << "EnvelopeDriftRecorder::initialize() - no valid nodes or perpendicular direction\n";
		  return 0;
	 }

	 //
	 // allocate memory
	 //
#ifdef _CSS
	 int nProcOuts;
	 int nVals = numNodes;
	 if (procGrpNum == -1)
		  if (procDataMethod != 0)
				nProcOuts = 1;
		  else
				nProcOuts = nVals;
	 else {
		  nProcOuts = nVals / procGrpNum;
		  if (nProcOuts * procGrpNum < nVals)
				nProcOuts++;
	 }
	 if (echoTimeFlag == true) {
		  currentData = new Vector(nProcOuts * 2); // additional data allocated for time  
		  data = new Matrix(3, nProcOuts * 2);
	 }
	 else
	 {
		  currentData = new Vector(nProcOuts);
		  data = new Matrix(3, nProcOuts);
	 }
#else
	 if (echoTimeFlag == true) {
		  currentData = new Vector(numNodes * 2); // additional data allocated for time  
		  data = new Matrix(3, numNodes * 2);
	 }
	 else {
		  currentData = new Vector(numNodes); // data(0) allocated for time  
		  data = new Matrix(3, numNodes);
	 }
#endif // _CSS

	 data->Zero();
	 theNodes = new Node * [2 * numNodes];
	 oneOverL = new Vector(numNodes);
	 if (theNodes == 0 || oneOverL == 0 || currentData == 0) {
		  opserr << "EnvelopeDriftRecorder::initialize() - out of memory\n";
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

	 first = true;

	 //
	 // mark as having been done & return
	 //

	 return 0;
}
//added by SAJalali
double EnvelopeDriftRecorder::getRecordedValue(int clmnId, int rowOffset, bool reset)
{
	double res = 0;
	if (!initializationDone)
		return res;
	if (clmnId >= data->noCols())
		return res;
	res = (*data)(2 - rowOffset, clmnId);
	if (reset)
		first = true;
	return res;
}

int EnvelopeDriftRecorder::flush(void) {
  if (theOutputHandler != 0) {
    return theOutputHandler->flush();
  }
  return 0;
}