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

// Written: fmk 
//
// Description: This file contains the class definition for EnvelopeNodeRecorder.
// A EnvelopeNodeRecorder is used to record the envelop of specified dof responses 
// at a collection of nodes over an analysis. (between commitTag of 0 and
// last commitTag).
//
// What: "@(#) EnvelopeNodeRecorder.C, revA"

#include <EnvelopeNodeRecorder.h>
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <FE_Datastore.h>
#include <FEM_ObjectBroker.h>
#include <MeshRegion.h>
#include <TimeSeries.h>

#include <StandardStream.h>
#include <DataFileStream.h>
#include <DataFileStreamAdd.h>
#include <XmlFileStream.h>
#include <BinaryFileStream.h>
#include <DatabaseStream.h>
#include <TCP_Stream.h>

#include <elementAPI.h>

#include <string.h>
#include <stdlib.h>
#include <math.h>


EnvelopeNodeRecorder::EnvelopeNodeRecorder()
	 :Recorder(RECORDER_TAGS_EnvelopeNodeRecorder),
	 theDofs(0), theNodalTags(0), theNodes(0),
	 currentData(0), data(0),
	 theDomain(0), theHandler(0),
	 first(true), initializationDone(false),
	 numValidNodes(0), addColumnInfo(0), theTimeSeries(0), timeSeriesValues(0)
{

}

EnvelopeNodeRecorder::EnvelopeNodeRecorder(const ID& dofs,
	 const ID* nodes,
	 const char* dataToStore,
	 Domain& theDom,
	 OPS_Stream* theOutputHandler,
	 int procMethod, int procGrpN,
	 bool echoTime,
	 TimeSeries** theSeries)
	 :Recorder(RECORDER_TAGS_EnvelopeNodeRecorder),
	 theDofs(0), theNodalTags(0), theNodes(0),
	 currentData(0), data(0),
	 theDomain(&theDom), theHandler(theOutputHandler),
	 first(true), initializationDone(false), numValidNodes(0), echoTimeFlag(echoTime),
	 addColumnInfo(0), theTimeSeries(theSeries), timeSeriesValues(0)
	 , procDataMethod(procMethod), procGrpNum(procGrpN), Modified(0)
{
	 // verify dof are valid 
#ifndef _CSS
	 int numDOF = dofs.Size();
#else
	 numDOF = dofs.Size();
#endif // !_CSS
	 theDofs = new ID(0, numDOF);

	 int count = 0;
	 int i;
	 for (i = 0; i < numDOF; i++) {
		  int dof = dofs(i);
		  if (dof >= 0) {
				(*theDofs)[count] = dof;
				count++;
		  }
		  else {
				opserr << "EnvelopeNodeRecorder::EnvelopeNodeRecorder - invalid dof  " << dof;
				opserr << " will be ignored\n";
		  }
	 }


	 // 
	 // create memory to hold nodal ID's (neeed parallel)
	 //

	 if (nodes != 0) {
		  int numNode = nodes->Size();
		  if (numNode != 0) {
				theNodalTags = new ID(*nodes);
				if (theNodalTags == 0 || theNodalTags->Size() != nodes->Size()) {
					 opserr << "NodeRecorder::NodeRecorder - out of memory\n";
				}
		  }
	 }

#ifndef _CSS
	 if (theTimeSeries != 0) {
		  timeSeriesValues = new double[numDOF];
		  for (int i = 0; i < numDOF; i++)
				timeSeriesValues[i] = 0.0;
	 }
#endif // !_CSS

	 //
	 // set the data flag used as a switch to get the response in a record
	 //

	 if (dataToStore == 0 || (strcmp(dataToStore, "disp") == 0)) {
		  dataFlag = 0;
	 }
	 else if ((strcmp(dataToStore, "vel") == 0)) {
		  dataFlag = 1;
	 }
	 else if ((strcmp(dataToStore, "accel") == 0)) {
		  dataFlag = 2;
	 }
	 else if ((strcmp(dataToStore, "incrDisp") == 0)) {
		  dataFlag = 3;
	 }
	 else if ((strcmp(dataToStore, "incrDeltaDisp") == 0)) {
		  dataFlag = 4;
	 }
	 else if ((strcmp(dataToStore, "unbalance") == 0)) {
		  dataFlag = 5;
	 }
	 else if ((strcmp(dataToStore, "unbalanceInclInertia") == 0) ||
		  (strcmp(dataToStore, "unbalanceIncInertia") == 0) ||
		  (strcmp(dataToStore, "unbalanceIncludingInertia") == 0)) {
		  dataFlag = 6;
	 }
	 else if ((strcmp(dataToStore, "reaction") == 0)) {
		  dataFlag = 7;
	 }
	 else if (((strcmp(dataToStore, "reactionIncInertia") == 0))
		  || ((strcmp(dataToStore, "reactionInclInertia") == 0))
		  || ((strcmp(dataToStore, "reactionIncludingInertia") == 0))) {
		  dataFlag = 8;
	 }
	 else if (((strcmp(dataToStore, "rayleighForces") == 0))
		  || ((strcmp(dataToStore, "rayleighDampingForces") == 0))) {
		  dataFlag = 9;

	 }
	 else if ((strcmp(dataToStore, "dispNorm") == 0)) {
		  dataFlag = 10000;

	 }
#ifdef _CSS
	 else if ((strcmp(dataToStore, "motionEnergy") == 0) || (strcmp(dataToStore, "MotionEnergy") == 0)) {
		  dataFlag = 999997;
		  numDOF = 1;
	 }
	 else if ((strcmp(dataToStore, "kineticEnergy") == 0) || (strcmp(dataToStore, "KineticEnergy") == 0)) {
		  dataFlag = 999998;
		  numDOF = 1;
	 }
	 else if ((strcmp(dataToStore, "dampingEnergy") == 0) || (strcmp(dataToStore, "DampingEnergy") == 0)) {
		  dataFlag = 999999;
		  numDOF = 1;
	 }
#endif _CSS
	 else {
		  dataFlag = 10;
		  opserr << "EnvelopeNodeRecorder::EnvelopeNodeRecorder - dataToStore " << dataToStore;
		  opserr << "not recognized (disp, vel, accel, incrDisp, incrDeltaDisp)\n";
	 }

	 if (theHandler != 0)
	 {
		  if (dataFlag == 7 || dataFlag == 8 || dataFlag == 9) {
				if (echoTime == true)
					 theHandler->setAddCommon(2);
				else
					 theHandler->setAddCommon(1);
		  }
	 }
#ifdef _CSS
	 if (theTimeSeries != 0) {
		  timeSeriesValues = new double[numDOF];
		  for (int i = 0; i < numDOF; i++)
				timeSeriesValues[i] = 0.0;
	 }
#endif // !_CSS
}


EnvelopeNodeRecorder::~EnvelopeNodeRecorder()
{
	 //
	 // write the data
	 //

	 if (theHandler != 0 && data != 0) {

		  theHandler->tag("Data"); // Data

		  int numResponse = data->noCols();

		  for (int i = 0; i < 3; i++) {
				for (int j = 0; j < numResponse; j++)
					 (*currentData)(j) = (*data)(i, j);
				theHandler->write(*currentData);
		  }

		  theHandler->endTag(); // Data
	 }

	 //
	 // clean up the memory
	 //

#ifdef _CSS
	 if (theDofs != 0) {
		  delete theDofs;
	 }
	 if (theTimeSeries != 0 && theNodes != 0)
	 {
		  for (int i = 0; i < numValidNodes; i++)
				theNodes[i]->resetTimeSeries();
	 }
	 if (dataFlag == 999997 || dataFlag == 999999)
	 {
		  for (int i = 0; i < numValidNodes; i++)
				theNodes[i]->resetEnergies();
	 }
#else
	 int numDOF;
	 if (theDofs != 0) {
		  numDOF = theDofs->Size();
		  delete theDofs;
	 }
#endif // _CSS

	 if (theNodalTags != 0)
		  delete theNodalTags;

	 if (theHandler != 0)
		  delete theHandler;

	 if (currentData != 0)
		  delete currentData;

	 if (data != 0)
		  delete data;

	 if (theNodes != 0)
		  delete[] theNodes;

	 if (theTimeSeries != 0) {
		  for (int i = 0; i < numDOF; i++)
				delete theTimeSeries[i];
		  delete[] theTimeSeries;
	 }
}

Vector EnvelopeNodeRecorder::getResponse(Node* theNode)
{
	//
	// if need nodal reactions get the domain to calculate them
	// before we iterate over the nodes
	//

	if (dataFlag == 7)
		theDomain->calculateNodalReactions(0);
	else if (dataFlag == 8)
		theDomain->calculateNodalReactions(1);
	if (dataFlag == 9)
		theDomain->calculateNodalReactions(2);

	double timeSeriesTerm = 0.0;
	Vector response;
	if (dataFlag == 0) {

		response.resize(numDOF);

		const Vector& theResponse = theNode->getTrialDisp();
		for (int j = 0; j < numDOF; j++) {

			if (theTimeSeries != 0) {
				timeSeriesTerm = timeSeriesValues[j];
			}

			int dof = (*theDofs)(j);
			if (theResponse.Size() > dof) {
				response[j] = theResponse(dof) + timeSeriesTerm;
			}
			else {
				response[j] = 0.0 + timeSeriesTerm;
			}
		}
		return response;
	}
	if (dataFlag == 1) {
		response.resize(numDOF);


		const Vector& theResponse = theNode->getTrialVel();
		for (int j = 0; j < numDOF; j++) {

			if (theTimeSeries != 0) {
				timeSeriesTerm = timeSeriesValues[j];
			}

			int dof = (*theDofs)(j);
			if (theResponse.Size() > dof) {
				response[j] = theResponse(dof) + timeSeriesTerm;
			}
			else
				response[j] = 0.0 + timeSeriesTerm;
		}
		return response;

	}
	if (dataFlag == 2) {
		response.resize(numDOF);


		const Vector& theResponse = theNode->getTrialAccel();
		for (int j = 0; j < numDOF; j++) {

			if (theTimeSeries != 0) {
				timeSeriesTerm = timeSeriesValues[j];
			}

			int dof = (*theDofs)(j);
			if (theResponse.Size() > dof) {
				response[j] = theResponse(dof) + timeSeriesTerm;
			}
			else
				response[j] = 0.0 + timeSeriesTerm;

		}
		return response;
	}
	if (dataFlag == 3) {
		response.resize(numDOF);

		const Vector& theResponse = theNode->getIncrDisp();
		for (int j = 0; j < numDOF; j++) {

			if (theTimeSeries != 0) {
				timeSeriesTerm = timeSeriesValues[j];
			}

			int dof = (*theDofs)(j);
			if (theResponse.Size() > dof) {
				response[j] = theResponse(dof);
			}
			else
				response[j] = 0.0;
		}
		return response;
	}
	if (dataFlag == 4) {
		response.resize(numDOF);

		const Vector& theResponse = theNode->getIncrDeltaDisp();
		for (int j = 0; j < numDOF; j++) {
			int dof = (*theDofs)(j);
			if (theResponse.Size() > dof) {
				response[j] = theResponse(dof);
			}
			else
				response[j] = 0.0;
		}
		return response;
	}
	if (dataFlag == 5) {
		response.resize(numDOF);

		const Vector& theResponse = theNode->getUnbalancedLoad();
		for (int j = 0; j < numDOF; j++) {
			int dof = (*theDofs)(j);
			if (theResponse.Size() > dof) {
				response[j] = theResponse(dof);
			}
			else
				response[j] = 0.0;
		}
		return response;
	}
	if (dataFlag == 6) {
		response.resize(numDOF);

		const Vector& theResponse = theNode->getUnbalancedLoadIncInertia();
		for (int j = 0; j < numDOF; j++) {
			int dof = (*theDofs)(j);
			if (theResponse.Size() > dof) {
				response[j] = theResponse(dof);
			}
			else
				response[j] = 0.0;
		}
		return response;
	}
	if (dataFlag == 7 || dataFlag == 8 || dataFlag == 9) {
		response.resize(numDOF);

		const Vector& theResponse = theNode->getReaction();
		for (int j = 0; j < numDOF; j++) {
			int dof = (*theDofs)(j);
			if (theResponse.Size() > dof) {
				response[j] = theResponse(dof);
			}
			else
				response[j] = 0.0;
		}
		return response;
	}
	if (dataFlag == 999997) {
		response.resize(1);

		if (theTimeSeries == 0)
		{
			opserr << "WARNING! NodeRecorder::motionEnergy: the timeSeries tag is missing. Please use the -TimeSeries option\n";
		}
		response[0] = theNode->getMotionEnergy(theTimeSeries);
		return response;
	}
	if (dataFlag == 999998)
	{
		response.resize(1);

		response[0] = theNode->getKineticEnergy(theTimeSeries);
		return response;
	}
	if (dataFlag == 999999)
	{
		response.resize(1);

		response[0] = theNode->getDampEnergy();
		return response;
	}
	return response;
}


int
EnvelopeNodeRecorder::record(int commitTag, double timeStamp)
{
	 if (theDomain == 0) {
		  return 0;
	 }
	 if (initializationDone != true)
		  if (this->initialize() != 0) {
				opserr << "EnvelopeNodeRecorder::record() - failed in initialize()\n";
				return -1;
		  }

	 double timeSeriesTerm = 0.0;

	 if (theTimeSeries != 0) {
		  for (int i = 0; i < numDOF; i++) {
				if (theTimeSeries[i] != 0)
					 timeSeriesValues[i] = theTimeSeries[i]->getFactor(timeStamp);
		  }
	 }

#ifdef _CSS
	 Modified = 0;
	 if (procDataMethod != 0)
	 {
		 Vector* respVecs = new Vector[numValidNodes];
		 for (int i = 0; i < numValidNodes; i++)
		 {
			 Node* theNode = theNodes[i];
			 respVecs[i] = getResponse(theNode);
		 }
		 int cnt = 0;
		  for (int j = 0; j < numDOF; j++) {
				int nProcOuts;
				int nVals = numValidNodes;
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
				int dof = (*theDofs)(j);
				for (int i = 0; i < numValidNodes; i++) {
					 Node* theNode = theNodes[i];
					 val1 = respVecs[i][j];
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
					 cnt = i * numDOF + j;
					 val = &vals[i];
					 (*currentData)(cnt) = *val;
				}
				delete[] vals;
		  }
		  delete[] respVecs;
	 }
	 else
#endif //_CSS
		  for (int i = 0; i < numValidNodes; i++) {
				int cnt = i * numDOF;

				if (dataFlag == 10000)
					 cnt = i;

				Node* theNode = theNodes[i];
				Vector resp = getResponse(theNode);
				for (int j = 0; j < resp.Size(); j++)
					(*currentData)(cnt++) = resp[j];
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
EnvelopeNodeRecorder::restart(void)
{
	 data->Zero();
	 first = true;
	 return 0;
}


int
EnvelopeNodeRecorder::setDomain(Domain& theDom)
{
	 theDomain = &theDom;
	 initializationDone = false;
	 return 0;
}


int
EnvelopeNodeRecorder::sendSelf(int commitTag, Channel& theChannel)
{
	 addColumnInfo = 1;

	 int numDOF = theDofs->Size();

	 if (theChannel.isDatastore() == 1) {
		  opserr << "EnvelopeNodeRecorder::sendSelf() - does not send data to a datastore\n";
		  return -1;
	 }

	 initializationDone = false;
	 static ID idData(7);
	 idData.Zero();

	 if (theDofs != 0)
		  idData(0) = numDOF;
	 if (theNodalTags != 0)
		  idData(1) = theNodalTags->Size();
	 if (theHandler != 0) {
		  idData(2) = theHandler->getClassTag();
	 }

	 idData(3) = dataFlag;

	 if (echoTimeFlag == true)
		  idData(4) = 1;
	 else
		  idData(4) = 0;

	 idData(5) = this->getTag();

	 if (theTimeSeries == 0)
		  idData[6] = 0;
	 else
		  idData[6] = 1;

	 if (theChannel.sendID(0, commitTag, idData) < 0) {
		  opserr << "EnvelopeNodeRecorder::sendSelf() - failed to send idData\n";
		  return -1;
	 }

	 if (theDofs != 0)
		  if (theChannel.sendID(0, commitTag, *theDofs) < 0) {
				opserr << "EnvelopeNodeRecorder::sendSelf() - failed to send dof id's\n";
				return -1;
		  }

	 if (theNodalTags != 0)
		  if (theChannel.sendID(0, commitTag, *theNodalTags) < 0) {
				opserr << "EnvelopeNodeRecorder::sendSelf() - failed to send nodal tags\n";
				return -1;
		  }

	 if (theHandler != 0)

		  if (theHandler->sendSelf(commitTag, theChannel) < 0) {
				opserr << "EnvelopeNodeRecorder::sendSelf() - failed to send the DataOutputHandler\n";
				return -1;
		  }

	 if (theTimeSeries != 0) {
		  ID timeSeriesTags(numDOF);
		  for (int i = 0; i < numDOF; i++) {
				if (theTimeSeries[i] != 0) {
					 timeSeriesTags[i] = theTimeSeries[i]->getClassTag();
				}
				else
					 timeSeriesTags[i] = -1;
		  }
		  if (theChannel.sendID(0, commitTag, timeSeriesTags) < 0) {
				opserr << "EnvelopeNodeRecorder::sendSelf() - failed to send time series tags\n";
				return -1;
		  }
		  for (int i = 0; i < numDOF; i++) {
				if (theTimeSeries[i] != 0) {
					 if (theTimeSeries[i]->sendSelf(commitTag, theChannel) < 0) {
						  opserr << "EnvelopeNodeRecorder::sendSelf() - time series failed in send\n";
						  return -1;

					 }
				}
		  }
	 }

	 return 0;
}

int
EnvelopeNodeRecorder::recvSelf(int commitTag, Channel& theChannel,
	 FEM_ObjectBroker& theBroker)
{
	 addColumnInfo = 1;

	 if (theChannel.isDatastore() == 1) {
		  opserr << "EnvelopeNodeRecorder::sendSelf() - does not send data to a datastore\n";
		  return -1;
	 }

	 static ID idData(7);
	 if (theChannel.recvID(0, commitTag, idData) < 0) {
		  opserr << "EnvelopeNodeRecorder::recvSelf() - failed to send idData\n";
		  return -1;
	 }

	 int numDOFs = idData(0);
	 int numNodes = idData(1);

	 dataFlag = idData(3);

	 this->setTag(idData(5));

	 if (idData(4) == 1)
		  echoTimeFlag = true;
	 else
		  echoTimeFlag = false;


	 //
	 // get the DOF ID data
	 //

	 if (theDofs == 0 || theDofs->Size() != numDOFs) {
		  if (theDofs != 0)
				delete theDofs;

		  if (numDOFs != 0) {
				theDofs = new ID(numDOFs);
				if (theDofs == 0 || theDofs->Size() != numDOFs) {
					 opserr << "EnvelopeNodeRecorder::recvSelf() - out of memory\n";
					 return -1;
				}
		  }
	 }
	 if (theDofs != 0)
		  if (theChannel.recvID(0, commitTag, *theDofs) < 0) {
				opserr << "EnvelopeNodeRecorder::recvSelf() - failed to recv dof data\n";
				return -1;
		  }

	 //
	 // get the NODAL tag data
	 //

	 if (theNodalTags == 0 || theNodalTags->Size() != numNodes) {
		  if (theNodalTags != 0)
				delete theNodalTags;

		  if (numNodes != 0) {
				theNodalTags = new ID(numNodes);
				if (theNodalTags == 0 || theNodalTags->Size() != numNodes) {
					 opserr << "EnvelopeNodeRecorder::recvSelf() - out of memory\n";
					 return -1;
				}
		  }
	 }
	 if (theNodalTags != 0)
		  if (theChannel.recvID(0, commitTag, *theNodalTags) < 0) {
				opserr << "EnvelopeNodeRecorder::recvSelf() - failed to recv dof data\n";
				return -1;
		  }


	 if (theHandler != 0)
		  delete theHandler;

	 theHandler = theBroker.getPtrNewStream(idData(2));
	 if (theHandler == 0) {
		  opserr << "EnvelopeNodeRecorder::sendSelf() - failed to get a data output handler\n";
		  return -1;
	 }

	 if (theHandler->recvSelf(commitTag, theChannel, theBroker) < 0) {
		  delete theHandler;
		  theHandler = 0;
	 }


	 if (idData[6] == 1) {
		  theTimeSeries = new TimeSeries * [numDOFs];
		  ID timeSeriesTags(numDOFs);
		  if (theChannel.recvID(0, commitTag, timeSeriesTags) < 0) {
				opserr << "EnvelopeNodeRecorder::recvSelf() - failed to recv time series tags\n";
				return -1;
		  }
		  for (int i = 0; i < numDOFs; i++) {
				if (timeSeriesTags[i] == -1)
					 theTimeSeries[i] = 0;
				else {
					 theTimeSeries[i] = theBroker.getNewTimeSeries(timeSeriesTags(i));
					 if (theTimeSeries[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
						  opserr << "EnvelopeNodeRecorder::recvSelf() - time series failed in recv\n";
						  return -1;
					 }
				}
		  }
	 }

	 return 0;
}


int
EnvelopeNodeRecorder::initialize(void)
{
#ifdef _CSS
	 if (theDomain == 0) {
		  opserr << "EnvelopeNodeRecorder::initialize() - either nodes or domain has not been set\n";
		  return -1;
	 }
#else
	 if (theDofs == 0 || theDomain == 0) {
		  opserr << "EnvelopeNodeRecorder::initialize() - either nodes, dofs or domain has not been set\n";
		  return -1;
	 }
#endif // _CSS

	 //
	 // create & set nodal array pointer
	 //

	 if (theNodes != 0)
		  delete[] theNodes;

	 numValidNodes = 0;

	 if (theNodalTags != 0) {
		  int numNode = theNodalTags->Size();
		  theNodes = new Node * [numNode];
		  if (theNodes == 0) {
				opserr << "EnvelopeNodeRecorder::domainChanged - out of memory\n";
				return -1;
}

		  for (int i = 0; i < numNode; i++) {
				int nodeTag = (*theNodalTags)(i);
				Node* theNode = theDomain->getNode(nodeTag);
				if (theNode != 0) {
					 theNodes[numValidNodes] = theNode;
					 numValidNodes++;
				}
		  }
	 }
	 else {

		  int numNodes = theDomain->getNumNodes();
		  if (numNodes != 0) {
				theNodes = new Node * [numNodes];

				if (theNodes == 0) {
					 opserr << "NodeRecorder::domainChanged - out of memory\n";
					 return -1;
				}
				NodeIter& theDomainNodes = theDomain->getNodes();
				Node* theNode;
				numValidNodes = 0;
				while (((theNode = theDomainNodes()) != 0) && (numValidNodes < numNodes)) {
					 theNodes[numValidNodes] = theNode;
					 numValidNodes++;
				}
		  }
		  else
				numValidNodes = 0;
	 }

	 //
	 // need to create the data description, i.e. what each column of data is
	 //

	 //
	 // need to create the data description, i.e. what each column of data is
	 //

	 char outputData[32];
	 char dataType[10];

	 if (dataFlag == 0) {
		  strcpy(dataType, "D");
	 }
	 else if (dataFlag == 1) {
		  strcpy(dataType, "V");
	 }
	 else if (dataFlag == 2) {
		  strcpy(dataType, "A");
	 }
	 else if (dataFlag == 3) {
		  strcpy(dataType, "dD");
	 }
	 else if (dataFlag == 4) {
		  strcpy(dataType, "ddD");
	 }
	 else if (dataFlag == 5) {
		  strcpy(dataType, "U");
	 }
	 else if (dataFlag == 6) {
		  strcpy(dataType, "U");
	 }
	 else if (dataFlag == 7) {
		  strcpy(dataType, "R");
	 }
	 else if (dataFlag == 8) {
		  strcpy(dataType, "R");
	 }
	 else if (dataFlag == 10000) {
		  strcpy(dataType, "|D|");
	 }
	 else if (dataFlag > 10) {
		  sprintf(dataType, "E%d", dataFlag - 10);
	 }
	 else
		  strcpy(dataType, "Unknown");


	 //
	 // resize the output matrix
	 //

#ifdef _CSS
	 int nProcOuts;
	 int nVals = numValidNodes;
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
	 int numValidResponse = nProcOuts * numDOF;
#else
	 int numDOF = theDofs->Size();
	 int numValidResponse = numValidNodes * numDOF;
#endif // _CSS

	 if (dataFlag == 10000)
		  numValidResponse = numValidNodes;

	 if (echoTimeFlag == true) {
		  numValidResponse *= 2;
	 }

	 currentData = new Vector(numValidResponse);
	 data = new Matrix(3, numValidResponse);
	 data->Zero();

	 ID dataOrder(numValidResponse);
	 ID xmlOrder(numValidNodes);

	 if (theNodalTags != 0 && addColumnInfo == 1) {

		  int count = 0;
		  int nodeCount = 0;

		  int numNode = theNodalTags->Size();
		  for (int i = 0; i < numNode; i++) {
				int nodeTag = (*theNodalTags)(i);
				Node* theNode = theDomain->getNode(nodeTag);
				if (theNode != 0) {
					 xmlOrder(nodeCount++) = i + 1;
					 for (int j = 0; j < numDOF; j++)
						  dataOrder(count++) = i + 1;
					 if (echoTimeFlag == true) {
						  for (int j = 0; j < numDOF; j++)
								dataOrder(count++) = i + 1;
					 }
				}
		  }
		  if (theHandler != 0)
				theHandler->setOrder(xmlOrder);
	 }

	 if (theHandler != 0)
		  for (int i = 0; i < numValidNodes; i++) {
				int nodeTag = theNodes[i]->getTag();

				theHandler->tag("NodeOutput");
				theHandler->attr("nodeTag", nodeTag);

				for (int j = 0; j < numDOF; j++) {

					 if (echoTimeFlag == true) {
						  theHandler->tag("TimeOutput");
						  theHandler->tag("ResponseType", "time");
						  theHandler->endTag();
					 }

					 sprintf(outputData, "%s%d", dataType, j + 1);
					 theHandler->tag("ResponseType", outputData);
				}

				theHandler->endTag();
		  }

	 if (theHandler != 0)
		  if (theNodalTags != 0 && addColumnInfo == 1) {
				theHandler->setOrder(dataOrder);
		  }

	 initializationDone = true;

	 return 0;
}
//added by SAJalali
double EnvelopeNodeRecorder::getRecordedValue(int clmnId, int rowOffset, bool reset)
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