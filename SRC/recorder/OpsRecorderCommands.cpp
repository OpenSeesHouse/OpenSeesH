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

 // $Revision: 1.58 $
 // $Date: 2010-05-27 17:51:54 $
 // $URL: $


 // Written: fmk 
 // Created: 04/98
 // Revision: AA
 //
 // Modofied SAJalali 1401/02
 // 
 // Description: This file contains the function that is invoked
 // by the interpreter when the comand 'record' is invoked by the 
 // user.
 //
 // What: "@(#) commands.C, revA"




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <EquiSolnAlgo.h>



// recorders
#include <NodeRecorder.h>
#include <EnvelopeNodeRecorder.h>
//#include <PatternRecorder.h>
#include <DriftRecorder.h>
#ifdef _CSS //by SAJalali
#include <ResidDriftRecorder.h>		
#include <ResidElementRecorder.h>
#include <ResidNodeRecorder.h>		
#include <ConditionalNodeRecorder.h>		
#include <ConditionalDriftRecorder.h>		
#include <ConditionalElementRecorder.h>		
#include <../tcl/FeViewer.h>
#endif // _CSS

#include <EnvelopeDriftRecorder.h>
#include <ElementRecorder.h>
#include <EnvelopeElementRecorder.h>
//#include <NormElementRecorder.h>
//#include <NormEnvelopeElementRecorder.h>

//#include <PVDRecorder.h>

//#include <GmshRecorder.h>
//#include <VTK_Recorder.h>

//extern void* OPS_PVDRecorder();

//extern void* OPS_GmshRecorder();
//#ifdef _HDF5
//extern void* OPS_MPCORecorder();
//#endif // _HDF5
//extern void* OPS_VTK_Recorder();
//extern void* OPS_ElementRecorderRMS();
//extern void* OPS_NodeRecorderRMS();


 //#include <NodeIter.h>
 //#include <ElementIter.h>
 //#include <Node.h>
 //#include <Element.h>
//#include <Parameter.h>
 //#include <DamageModel.h>
 //#include <DamageRecorder.h>
#include <MeshRegion.h>
//#include <GSA_Recorder.h>
//#include <RemoveRecorder.h>

//#include <TclModelBuilder.h>

#include <StandardStream.h>
#include <DataFileStream.h>
#include <DataFileStreamAdd.h>
#include <XmlFileStream.h>
#include <BinaryFileStream.h>
#include <DatabaseStream.h>
#include <DummyStream.h>
#include <TCP_Stream.h>

#include <packages.h>
#include <elementAPI.h>


//extern TimeSeries *TclSeriesCommand(ClientData clientData, OPS_Interp * TCL_Char *arg);

//extern const char * getInterpPWD(OPS_Interp *interp);  // commands.cpp

//extern TclModelBuilder *theDamageTclModelBuilder;

typedef struct externalRecorderCommand {
	char* funcName;
	void* (*funcPtr)();
	struct externalRecorderCommand* next;
} ExternalRecorderCommand;

static ExternalRecorderCommand* theExternalRecorderCommands = NULL;
static FeViewer* theFeViewr = 0;

#ifdef _CSS
enum outputMode { No_Output, STANDARD_STREAM, DATA_STREAM, XML_STREAM, DATABASE_STREAM, BINARY_STREAM, DATA_STREAM_CSV, TCP_STREAM, DATA_STREAM_ADD };
#else
enum outputMode { STANDARD_STREAM, DATA_STREAM, XML_STREAM, DATABASE_STREAM, BINARY_STREAM, DATA_STREAM_CSV, TCP_STREAM, DATA_STREAM_ADD };
#endif

#include <SimulationInformation.h>
extern SimulationInformation* OPS_GetSimulationInfo();

void* OPS_NodeRecorder(const char* type)
{
	if (OPS_GetNumRemainingInputArgs() < 5) {
		opserr << "WARING: recorder Node ";
		opserr << "-node <list nodes> -dof <doflist> -file <fileName> -dT <dT> reponse\n";
		opserr << "insufficient number of inputs: " << OPS_GetNumRemainingInputArgs() << "\n";
		return 0;
	}

	const char* responseID = 0;
	OPS_Stream* theOutputStream = 0;
	const char* filename = 0;


#ifdef _CSS
	int procDataMethod = 0;
	int nProcGrp = -1;
	int cntrlRcrdrTag = 0;
	outputMode eMode = No_Output;
#else
	int eMode = STANDARD_STREAM;
#endif // _CSS

	bool echoTimeFlag = false;
	double dT = 0.0;
	bool doScientific = false;

	int precision = 6;

	bool closeOnWrite = false;

	const char* inetAddr = 0;
	int inetPort;

	int gradIndex = -1;

	TimeSeries** theTimeSeries = 0;

	ID nodes(0, 6);
	ID dofs(0, 6);
	ID timeseries(0, 6);
	int numData = 1;
	while (OPS_GetNumRemainingInputArgs() > 0) {

		const char* option = OPS_GetString();
		responseID = option;

		if (strcmp(option, "-time") == 0 || strcmp(option, "-load") == 0) {
			echoTimeFlag = true;
		}
		else if (strcmp(option, "-scientific") == 0) {
			doScientific = true;
		}
		else if (strcmp(option, "-file") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				filename = OPS_GetString();
			}
			eMode = DATA_STREAM;
		}
		else if (strcmp(option, "-fileAdd") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				filename = OPS_GetString();
			}
			eMode = DATA_STREAM_ADD;
		}
		else if (strcmp(option, "-closeOnWrite") == 0) {
			closeOnWrite = true;
		}
		else if (strcmp(option, "-csv") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				filename = OPS_GetString();
			}
			eMode = DATA_STREAM_CSV;
		}
		else if (strcmp(option, "-tcp") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				inetAddr = OPS_GetString();
			}
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetIntInput(&numData, &inetPort) < 0) {
					opserr << "WARNING: failed to read inetPort\n";
					return 0;
				}
			}
			eMode = TCP_STREAM;
		}
		else if (strcmp(option, "-xml") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				filename = OPS_GetString();
			}
			eMode = XML_STREAM;
		}
		else if (strcmp(option, "-binary") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				filename = OPS_GetString();
			}
			eMode = BINARY_STREAM;
		}
		else if (strcmp(option, "-dT") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDoubleInput(&numData, &dT) < 0) {
					opserr << "WARNING: failed to read dT\n";
					return 0;
				}
			}
		}
		else if (strcmp(option, "-rTolDt") == 0) {
			opserr << "WARNING: " << type << ": rTolDt Option is constantly assumed as 1.e-5 in this version\n";
		}
		else if (strcmp(option, "-timeSeries") == 0) {
			int numTimeSeries = 0;
			while (OPS_GetNumRemainingInputArgs() > 0) {

				int ts;
				if (OPS_GetIntInput(&numData, &ts) < 0) {
					break;
				}
				timeseries[numTimeSeries++] = ts;
			}
			theTimeSeries = new TimeSeries * [numTimeSeries];
			for (int j = 0; j < numTimeSeries; j++) {
				if (timeseries(j) != 0 && timeseries(j) != -1)
					theTimeSeries[j] = OPS_getTimeSeries(timeseries(j));
				else
					theTimeSeries[j] = 0;
			}
		}
		else if (strcmp(option, "-precision") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetIntInput(&numData, &precision) < 0) {
					opserr << "WARNING: failed to read precision\n";
					return 0;
				}
			}
		}
		else if (strcmp(option, "-node") == 0 || strcmp(option, "-nodes") == 0) {
			int numNodes = 0;
			while (OPS_GetNumRemainingInputArgs() > 0) {

				int nd;
				if (OPS_GetIntInput(&numData, &nd) < 0) {
					OPS_ResetCurrentInputArg(-1);
					break;
				}
				nodes[numNodes++] = nd;
			}
		}
		else if (strcmp(option, "-nodeRange") == 0) {
			int start = 0, end = 0;
			if (OPS_GetNumRemainingInputArgs() > 0) {

				if (OPS_GetIntInput(&numData, &start) < 0) {
					opserr << "WARNING: failed to read start node\n";
					return 0;
				}
			}
			if (OPS_GetNumRemainingInputArgs() > 0) {

				if (OPS_GetIntInput(&numData, &end) < 0) {
					opserr << "WARNING: failed to read end node\n";
					return 0;
				}
			}
			if (start == 0 || end == 0)
			{
				opserr << "WARNING: nodeRegion: start or end tags are zero\n";
				return 0;
			}

			if (start > end) {
				int swap = end;
				end = start;
				start = swap;
			}
			int numNodes = 0;
			for (int i = start; i <= end; i++)
				nodes[numNodes++] = i;
		}
		else if (strcmp(option, "-region") == 0) {
			int tag = 0;
			if (OPS_GetNumRemainingInputArgs() > 0) {

				if (OPS_GetIntInput(&numData, &tag) < 0) {
					opserr << "WARNING: failed to read region tag\n";
					return 0;
				}
			}
			Domain* domain = OPS_GetDomain();
			MeshRegion* theRegion = domain->getRegion(tag);
			if (theRegion == 0) {
				opserr << "WARNING: region does not exist\n";
				return 0;
			}
			const ID& nodeRegion = theRegion->getNodes();
			int numNodes = 0;
			for (int i = 0; i < nodeRegion.Size(); i++)
				nodes[numNodes++] = nodeRegion(i);
		}
		else if (strcmp(option, "-dof") == 0) {
			int numDOF = 0;
			while (OPS_GetNumRemainingInputArgs() > 0) {

				int dof;
				if (OPS_GetIntInput(&numData, &dof) < 0) {
					OPS_ResetCurrentInputArg(-1);
					break;
				}
				dofs[numDOF++] = dof - 1;
			}
		}
#ifdef _CSS
		else if (strcmp(option, "-controlRecorder") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {

				if (OPS_GetIntInput(&numData, &cntrlRcrdrTag) < 0) {
					opserr << "WARNING recorder Node -controlRecorder $ntrlRcrdrTag - invalid ntrlRcrdrTag " << endln;
					return 0;
				}
			}
		}
		else if (strcmp(option, "-process") == 0) {
			const char* procType = OPS_GetString();
			if (strcmp(procType, "sum") == 0)
				procDataMethod = 1;
			else if (strcmp(procType, "max") == 0)
				procDataMethod = 2;
			else if (strcmp(procType, "min") == 0)
				procDataMethod = 3;
			else if (strcmp(procType, "maxAbs") == 0)
				procDataMethod = 4;
			else if (strcmp(procType, "minAbs") == 0)
				procDataMethod = 5;
			else
			{
				opserr << "unrecognized element result process method: " << procType << endln;
				return 0;
			}
		}
		else if (strcmp(option, "-procGrpNum") == 0) {
			int num = 1;
			if (OPS_GetIntInput(&numData, &nProcGrp) < 0) {
				opserr << "Failed to read $nGrp after -procGrpNum option in recorder" << endln;
				return 0;
			}
		}
#endif // _CSS
	}

	// data handler
	if (eMode == DATA_STREAM && filename != 0)
		theOutputStream = new DataFileStream(filename, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
	else if (eMode == DATA_STREAM_ADD && filename != 0)
		theOutputStream = new DataFileStreamAdd(filename, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
	else if (eMode == DATA_STREAM_CSV && filename != 0)
		theOutputStream = new DataFileStream(filename, OVERWRITE, 2, 1, closeOnWrite, precision, doScientific);
	else if (eMode == XML_STREAM && filename != 0)
		theOutputStream = new XmlFileStream(filename);
	//else if (eMode == DATABASE_STREAM && tableName != 0)
	//    theOutputStream = new DatabaseStream(theDatabase, tableName);
	else if (eMode == BINARY_STREAM && filename != 0)
		theOutputStream = new BinaryFileStream(filename);
	else if (eMode == TCP_STREAM && inetAddr != 0)
		theOutputStream = new TCP_Stream(inetPort, inetAddr);
	if (theOutputStream != 0)
		theOutputStream->setPrecision(precision);

	Domain* domain = OPS_GetDomain();
	if (domain == 0)
		return 0;
	Recorder* recorder = 0;
	if (strcmp(type, "Node") == 0)
		recorder = new NodeRecorder(dofs, &nodes, gradIndex, responseID, *domain, theOutputStream,
			procDataMethod, nProcGrp, dT, echoTimeFlag, theTimeSeries);
	else if (strcmp(type, "EnvelopeNode") == 0)
		recorder = new EnvelopeNodeRecorder(dofs, &nodes, responseID, *domain, theOutputStream,
			procDataMethod, nProcGrp, echoTimeFlag, theTimeSeries);
	else if (strcmp(type, "ResidNode") == 0)
		recorder = new ResidNodeRecorder(dofs, &nodes, responseID, *domain, theOutputStream,
			procDataMethod, nProcGrp, echoTimeFlag, theTimeSeries);
	else if (strcmp(type, "ConditionalNode") == 0)
		recorder = new ConditionalNodeRecorder(dofs, &nodes, responseID, *domain, theOutputStream,
			cntrlRcrdrTag, procDataMethod, nProcGrp, echoTimeFlag, theTimeSeries);

	return recorder;
}

void* OPS_ElementRecorder(const char* type)
{
	if (OPS_GetNumRemainingInputArgs() < 5) {
		opserr << "WARING: recorder Element ";
		opserr << "-ele <list elements> -file <fileName> -dT <dT> reponse";
		return 0;
	}

	const char** data = 0;
	int nargrem = 0;
	OPS_Stream* theOutputStream = 0;
	const char* filename = 0;


	outputMode eMode = No_Output;

	bool echoTimeFlag = false;
	double dT = 0.0;
	double rTolDt = 0.00001;
	bool doScientific = false;

	int precision = 6;

	bool closeOnWrite = false;

	const char* inetAddr = 0;
	int inetPort;

	ID elements(0, 6);
	ID dofs(0, 6);
	int procDataMethod = 0;
	int nProcGrp = -1;
	int cntrlRcrdrTag = 0;

	char** argv = 0;
	int numData = 1;
	while (OPS_GetNumRemainingInputArgs() > 0) {

		const char* option = OPS_GetString();

		if (strcmp(option, "-time") == 0 || strcmp(option, "-load") == 0) {
			echoTimeFlag = true;
		}
		else if (strcmp(option, "-scientific") == 0) {
			doScientific = true;
		}
		else if (strcmp(option, "-file") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				filename = OPS_GetString();
			}
			eMode = DATA_STREAM;
		}
		else if (strcmp(option, "-fileAdd") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				filename = OPS_GetString();
			}
			eMode = DATA_STREAM_ADD;
		}
		else if (strcmp(option, "-closeOnWrite") == 0) {
			closeOnWrite = true;
		}
		else if (strcmp(option, "-csv") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				filename = OPS_GetString();
			}
			eMode = DATA_STREAM_CSV;
		}
		else if (strcmp(option, "-tcp") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				inetAddr = OPS_GetString();
			}
			if (OPS_GetNumRemainingInputArgs() > 0) {

				if (OPS_GetIntInput(&numData, &inetPort) < 0) {
					opserr << "WARNING: failed to read inetPort\n";
					return 0;
				}
			}
			eMode = TCP_STREAM;
		}
		else if (strcmp(option, "-xml") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				filename = OPS_GetString();
			}
			eMode = XML_STREAM;
		}
		else if (strcmp(option, "-binary") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				filename = OPS_GetString();
			}
			eMode = BINARY_STREAM;
		}
		else if (strcmp(option, "-dT") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {

				if (OPS_GetDoubleInput(&numData, &dT) < 0) {
					opserr << "WARNING: failed to read dT\n";
					return 0;
				}
			}
		}
		else if (strcmp(option, "-rTolDt") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {

				if (OPS_GetDoubleInput(&numData, &rTolDt) < 0) {
					opserr << "WARNING: failed to read rTolDt\n";
					return 0;
				}
			}
		}
		else if (strcmp(option, "-precision") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {

				if (OPS_GetIntInput(&numData, &precision) < 0) {
					opserr << "WARNING: failed to read precision\n";
					return 0;
				}
			}
		}
		else if (strcmp(option, "-ele") == 0) {
			int numEle = 0;
			while (OPS_GetNumRemainingInputArgs() > 0) {

				int el;
				if (OPS_GetIntInput(&numData, &el) < 0) {
					OPS_ResetCurrentInputArg(-1);
					break;
				}
				elements[numEle++] = el;
			}
		}
		else if (strcmp(option, "-eleRange") == 0) {
			int start, end;
			if (OPS_GetNumRemainingInputArgs() > 0) {

				if (OPS_GetIntInput(&numData, &start) < 0) {
					opserr << "WARNING: failed to read start element\n";
					return 0;
				}
			}
			if (OPS_GetNumRemainingInputArgs() > 0) {

				if (OPS_GetIntInput(&numData, &end) < 0) {
					opserr << "WARNING: failed to read end element\n";
					return 0;
				}
			}
			if (start == 0 || end == 0)
			{
				opserr << "WARNING! start or end region tags are zero\n";
				return 0;
			}
			if (start > end) {
				int swap = end;
				end = start;
				start = swap;
			}
			int numEle = 0;
			for (int i = start; i <= end; i++)
				elements[numEle++] = i;
		}
		else if (strcmp(option, "-region") == 0) {
			int tag = 0;
			if (OPS_GetNumRemainingInputArgs() > 0) {

				if (OPS_GetIntInput(&numData, &tag) < 0) {
					opserr << "WARNING: failed to read region tag\n";
					return 0;
				}
			}
			Domain* domain = OPS_GetDomain();
			MeshRegion* theRegion = domain->getRegion(tag);
			if (theRegion == 0) {
				opserr << "WARNING: region does not exist\n";
				return 0;
			}
			const ID& eleRegion = theRegion->getElements();
			int numEle = 0;
			for (int i = 0; i < eleRegion.Size(); i++)
				elements[numEle++] = eleRegion(i);
		}
		else if (strcmp(option, "-dof") == 0) {
			int numDOF = 0;
			while (OPS_GetNumRemainingInputArgs() > 0) {

				int dof;
				if (OPS_GetIntInput(&numData, &dof) < 0) {
					OPS_ResetCurrentInputArg(-1);
					break;
				}
				dofs[numDOF++] = dof - 1;
			}
		}
		else if (strcmp(option, "-controlRecorder") == 0) {
			int num = 0;
			if (OPS_GetInt(&numData, &cntrlRcrdrTag) != 0) {
				opserr << "WARNING recorder ConditionalElement -controlRecorder $ntrlRcrdrTag - invalid cntrlRcrdrTag " << option << endln;
				return 0;
			}
		}
		else if (strcmp(option, "-process") == 0) {
			const char* procType = OPS_GetString();
			if (strcmp(procType, "sum") == 0)
				procDataMethod = 1;
			else if (strcmp(procType, "max") == 0)
				procDataMethod = 2;
			else if (strcmp(procType, "min") == 0)
				procDataMethod = 3;
			else if (strcmp(procType, "maxAbs") == 0)
				procDataMethod = 4;
			else if (strcmp(procType, "minAbs") == 0)
				procDataMethod = 5;
			else
				opserr << "unrecognized element result process method: " << procType << endln;
		}
		else if (strcmp(option, "-procGrpNum") == 0) {
			int num = 0;
			if (OPS_GetInt(&numData, &nProcGrp) != 0) {
				opserr << "Failed to read $nGrp after -procGrpNum option in recorder" << option << endln;
				return 0;
			}
		}
		else {
			// first unknown string then is assumed to start 
									// element response request
			nargrem = 1 + OPS_GetNumRemainingInputArgs();
			data = new const char* [nargrem];
			data[0] = option;
			//argv = new char*[nargrem];
			for (int i = 1; i < nargrem; i++) {
				data[i] = new char[128];

				// Turn everything in to a string for setResponse
				OPS_GetStringFromAll((char*)data[i], 128);
			}
		}
	}

	// data handler
	if (eMode == DATA_STREAM && filename != 0)
		theOutputStream = new DataFileStream(filename, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
	else if (eMode == DATA_STREAM_ADD && filename != 0)
		theOutputStream = new DataFileStreamAdd(filename, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
	else if (eMode == DATA_STREAM_CSV && filename != 0)
		theOutputStream = new DataFileStream(filename, OVERWRITE, 2, 1, closeOnWrite, precision, doScientific);
	else if (eMode == XML_STREAM && filename != 0)
		theOutputStream = new XmlFileStream(filename);
	//else if (eMode == DATABASE_STREAM && tableName != 0)
	//    theOutputStream = new DatabaseStream(theDatabase, tableName);
	else if (eMode == BINARY_STREAM && filename != 0)
		theOutputStream = new BinaryFileStream(filename);
	else if (eMode == TCP_STREAM && inetAddr != 0)
		theOutputStream = new TCP_Stream(inetPort, inetAddr);
	if (theOutputStream != 0)
		theOutputStream->setPrecision(precision);

	Domain* domain = OPS_GetDomain();
	if (domain == 0)
		return 0;
	Recorder* recorder = 0;
	if (strcmp(type, "Element") == 0)
		recorder = new ElementRecorder(&elements, data, nargrem, echoTimeFlag, *domain,
			theOutputStream, procDataMethod, nProcGrp, dT, &dofs);
	if (strcmp(type, "EnvelopeElement") == 0)
		recorder = new EnvelopeElementRecorder(&elements, data, nargrem, *domain,
			theOutputStream, procDataMethod, nProcGrp, echoTimeFlag, &dofs);
	if (strcmp(type, "ResidElement") == 0)
		recorder = new ResidElementRecorder(&elements, data, nargrem, *domain,
			theOutputStream, procDataMethod, nProcGrp, echoTimeFlag, &dofs);
	if (strcmp(type, "ConditionalElement") == 0)
		recorder = new ConditionalElementRecorder(&elements, data, nargrem, *domain,
			theOutputStream, cntrlRcrdrTag, procDataMethod, nProcGrp, echoTimeFlag, &dofs);

	if (data != 0) {
		delete[] data;
	}

	return recorder;
}

void* OPS_DriftRecorder(const char* type)
{

	outputMode eMode = No_Output;       // enum found in DataOutputFileHandler.h

	bool echoTimeFlag = false;
	ID iNodes(0, 16);
	ID jNodes(0, 16);
	int dof = 1;
	int perpDirn = 2;
	int pos = 2;
	double dT = 0.0;
	int precision = 6;
	bool doScientific = false;
	bool closeOnWrite = false;
	int cntrlRcrdrTag = 0;
	int procDataMethod = 0;
	int nProcGrp = -1;
	const char* fileName;
	int numData = 1;
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* option = OPS_GetString();

		if (strcmp(option, "-file") == 0) {
			fileName = OPS_GetString();
			eMode = DATA_STREAM;
			const char* pwd = OPS_GetInterpPWD();
			OPS_GetSimulationInfo()->addOutputFile(fileName, pwd);
		}
		else if (strcmp(option, "-xml") == 0) {
			fileName = OPS_GetString();
			eMode = XML_STREAM;
		}
		else if (strcmp(option, "-headings") == 0) {
			eMode = DATA_STREAM;
			fileName = OPS_GetString();
		}
#ifdef _CSS
		else if (strcmp(option, "-controlRecorder") == 0) {
			int num = 0;
			if (OPS_GetInt(&numData, &cntrlRcrdrTag) != 0) {
				opserr << "WARNING recorder ConditionalNode -controlRecorder $ntrlRcrdrTag - invalid ntrlRcrdrTag " << endln;
				return 0;
			}
		}
		else if (strcmp(option, "-process") == 0) {
			const char* procType = OPS_GetString();
			if (strcmp(procType, "sum") == 0)
				procDataMethod = 1;
			else if (strcmp(procType, "max") == 0)
				procDataMethod = 2;
			else if (strcmp(procType, "min") == 0)
				procDataMethod = 3;
			else if (strcmp(procType, "maxAbs") == 0)
				procDataMethod = 4;
			else if (strcmp(procType, "minAbs") == 0)
				procDataMethod = 5;
			else
				opserr << "unrecognized element result process method: " << procType << endln;
			pos += 2;
		}
		else if (strcmp(option, "-procGrpNum") == 0) {
			int num = 0;
			if (OPS_GetInt(&numData, &nProcGrp) != 0) {
				opserr << "Failed to read $nGrp after -procGrpNum option in recorder" << endln;
				return 0;
			}
		}
#endif // _CSS
		else if (strcmp(option, "-fileCSV") == 0) {
			fileName = OPS_GetString();
			eMode = DATA_STREAM_CSV;
			const char* pwd = OPS_GetInterpPWD();
			OPS_GetSimulationInfo()->addOutputFile(fileName, pwd);
			pos += 2;
		}

		else if ((strcmp(option, "-binary") == 0)) {
			// allow user to specify load pattern other than current
			fileName = OPS_GetString();
			const char* pwd = OPS_GetInterpPWD();
			OPS_GetSimulationInfo()->addOutputFile(fileName, pwd);
			eMode = BINARY_STREAM;
			pos += 2;
		}

		else if ((strcmp(option, "-nees") == 0) || (strcmp(option, "-xml") == 0)) {
			// allow user to specify load pattern other than current
			fileName = OPS_GetString();
			const char* pwd = OPS_GetInterpPWD();
			OPS_GetSimulationInfo()->addOutputFile(fileName, pwd);
			eMode = XML_STREAM;
		}

		else if (strcmp(option, "-closeOnWrite") == 0) {
			closeOnWrite = true;
		}

		else if (strcmp(option, "-scientific") == 0) {
			doScientific = true;
		}

		else if (strcmp(option, "-precision") == 0) {
			int num = 0;
			if (OPS_GetInt(&numData, &precision) != 0)
			{
				opserr << "inavlid precision, wants integer number\n";
				return 0;
			}
		}

		else if ((strcmp(option, "-iNode") == 0) || (strcmp(option, "-iNodes") == 0)) {
			int node;
			int numNodes = 0;
			while (OPS_GetNumRemainingInputArgs() > 0)
			{
				if (OPS_GetInt(&numData, &node) != 0) {
					OPS_ResetCurrentInputArg(-1);
					break;
				}
				iNodes[numNodes++] = node;
			}
		}

		else if ((strcmp(option, "-jNode") == 0) || (strcmp(option, "-jNodes") == 0)) {
			int node;
			int numNodes = 0;
			while (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetInt(&numData, &node) != 0) {
					OPS_ResetCurrentInputArg(-1);
					break;
				}
				jNodes[numNodes++] = node;
			}
		}

		else if (strcmp(option, "-dof") == 0) {
			if (OPS_GetInt(&numData, &dof) != 0) {
				opserr << "WARNING: INVALID dof number\n";
				return 0;
			}
		}

		else if (strcmp(option, "-perpDirn") == 0) {
			int num = 0;
			if (OPS_GetInt(&numData, &perpDirn) != 0) {
				opserr << "WARNING: INVALID perpDirn\n";
				return 0;
			}
		}

		else if (strcmp(option, "-time") == 0) {
			echoTimeFlag = true;
		}

		else if (strcmp(option, "-dT") == 0) {
			// allow user to specify time step size for recording
			int num = 0;
			if (OPS_GetDouble(&numData, &dT) != 0) {
				opserr << "WARNING: INVALID perpDirn\n";
				return 0;
			}
		}
	}
	if (iNodes.Size() != jNodes.Size()) {
		opserr << "WARNING recorder Drift - the number of iNodes and jNodes must be the same " << iNodes << " " << jNodes << endln;
		return 0;
	}

	OPS_Stream* theOutputStream = 0;

	// construct the DataHandler
	if (eMode == DATA_STREAM && fileName != 0) {
		theOutputStream = new DataFileStream(fileName, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
	}
	else if (eMode == DATA_STREAM_CSV && fileName != 0) {
		theOutputStream = new DataFileStream(fileName, OVERWRITE, 2, 1, closeOnWrite, precision, doScientific);
	}
	else if (eMode == XML_STREAM && fileName != 0) {
		theOutputStream = new XmlFileStream(fileName);
	} /*else if (eMode == DATABASE_STREAM && tableName != 0) {
theOutputStream = new DatabaseStream(theDatabase, tableName);
	}*/ else if (eMode == BINARY_STREAM) {
		theOutputStream = new BinaryFileStream(fileName);
	}
	if (theOutputStream != 0)
		theOutputStream->setPrecision(precision);
	// Subtract one from dof and perpDirn for C indexing
	Recorder* theRecorder = 0;
	Domain* theDomain = OPS_GetDomain();
	if (strcmp(type, "ResidDrift") == 0)
		theRecorder = new ResidDriftRecorder(iNodes, jNodes, dof - 1, perpDirn - 1,
			*theDomain, theOutputStream, procDataMethod, nProcGrp, echoTimeFlag);
	else if (strcmp(type, "ConditionalDrift") == 0)
		theRecorder = new ConditionalDriftRecorder(iNodes, jNodes, dof - 1, perpDirn - 1,
			*theDomain, theOutputStream, cntrlRcrdrTag, procDataMethod, nProcGrp, echoTimeFlag);
	else if (strcmp(type, "Drift") == 0)
		theRecorder = new DriftRecorder(iNodes, jNodes, dof - 1, perpDirn - 1,
			*theDomain, theOutputStream, procDataMethod, nProcGrp, echoTimeFlag, dT);
	else if (strcmp(type, "EnvelopeDrift") == 0)
		theRecorder = new EnvelopeDriftRecorder(iNodes, jNodes, dof - 1, perpDirn - 1,
			*theDomain, theOutputStream, procDataMethod, nProcGrp, echoTimeFlag);
	else
		opserr << "WARNING! Unrecognized Drift Recorder type: " << type << endln;
	return theRecorder;
}

void* OPS_DisplayRecorder(const char* type = "")
{
	int xLoc, yLoc, width, height;
	int pos = 7;
	int wipeFlag = 0;
	double dT = 0.0;
	double rTolDt = 0.00001;
	int saveToFile = 0;
	const char* fileName = "";
	int argc = OPS_GetNumRemainingInputArgs();
	if (argc < 6) {
		opserr << "WARNING recorder display title xLoc yLoc pixelsX pixelsY <-wipe> <-dT deltaT?> <-file fileName?>";
		return 0;
	}
	const char* title = OPS_GetString();
	int numData = 1;
	if (OPS_GetInt(&numData, &xLoc) != 0)
		return 0;
	if (OPS_GetInt(&numData, &yLoc) != 0)
		return 0;
	if (OPS_GetInt(&numData, &width) != 0)
		return 0;
	if (OPS_GetInt(&numData, &height) != 0)
		return 0;
	const char* option;
	while (OPS_GetNumRemainingInputArgs() > 0) {
		option = OPS_GetString();
		// check if we are to wipe image on each redraw
		if (strcmp(option, "-wipe") == 0) {
			wipeFlag = 1;
		}
		// allow user to specify time step size for recording
		else if (strcmp(option, "-dT") == 0) {
			if (OPS_GetDouble(&numData, &dT) != 0)
				return 0;
		}
		// allow user to specify time step size for recording
		else if (strcmp(option, "-rTolDt") == 0) {
			if (OPS_GetDouble(&numData, &rTolDt) != 0)
				return 0;
		}
		// save images to file
		else if (strcmp(option, "-file") == 0) {
			fileName = OPS_GetString();
			saveToFile = 1;
		}
	}
	Domain* theDomain = OPS_GetDomain();
	Recorder* theRecorder = 0;
	if (!saveToFile)
		theRecorder = new FeViewer(title, xLoc, yLoc, width, height, *theDomain, wipeFlag, dT, rTolDt);
	else
		theRecorder = new FeViewer(title, xLoc, yLoc, width, height, fileName, *theDomain, dT, rTolDt);
	theFeViewr = (FeViewer*)theRecorder;
	return theRecorder;

}

FeViewer* OPS_GetFeViewer()
{
	return theFeViewr;
}
