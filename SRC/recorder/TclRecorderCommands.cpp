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


#include <tcl.h>


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


#include <NodeIter.h>
#include <ElementIter.h>
#include <Node.h>
#include <Element.h>
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

//#include <packages.h>
#include <elementAPI.h>
extern "C" int         OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp* interp, int cArg, int mArg, TCL_Char** argv, Domain* domain);


extern TimeSeries* TclSeriesCommand(ClientData clientData, Tcl_Interp* interp, TCL_Char* arg);

extern const char* getInterpPWD(Tcl_Interp* interp);  // commands.cpp

//extern TclModelBuilder *theDamageTclModelBuilder;

//typedef struct externalRecorderCommand {
//	 char* funcName;
//	 void* (*funcPtr)();
//	 struct externalRecorderCommand* next;
//} ExternalRecorderCommand;
//
//static ExternalRecorderCommand* theExternalRecorderCommands = NULL;

#ifdef _CSS
enum outputMode { No_Output, STANDARD_STREAM, DATA_STREAM, XML_STREAM, DATABASE_STREAM, BINARY_STREAM, DATA_STREAM_CSV, TCP_STREAM, DATA_STREAM_ADD };
#else
enum outputMode { STANDARD_STREAM, DATA_STREAM, XML_STREAM, DATABASE_STREAM, BINARY_STREAM, DATA_STREAM_CSV, TCP_STREAM, DATA_STREAM_ADD };
#endif


#include <FeViewer.h>

#ifdef _NOGRAPHICS

#else

#include <FilePlotter.h>
//#include <AlgorithmIncrements.h>

#endif


#include <SimulationInformation.h>
extern SimulationInformation simulationInfo;
void* OPS_NodeRecorder(const char* type);
void* OPS_ElementRecorder(const char* type);
void* OPS_DriftRecorder(const char* type);

static EquiSolnAlgo* theAlgorithm = 0;
extern FE_Datastore* theDatabase;
extern FEM_ObjectBroker theBroker;
FeViewer* OPS_GetFeViewer();
int
FeViewer_setVRP(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return TCL_OK;
#else
	// check destructor has not been called
	if (OPS_GetFeViewer() == 0)
		return TCL_OK;

	// ensure corrcet num arguments
	if (argc < 4) {
		opserr << "WARNING args incorrect - vrp xloc yloc zloc \n";
		return TCL_ERROR;
	}

	// get the xLoc, yLoc and zLoc
	double xLoc, yLoc, zLoc;
	if (Tcl_GetDouble(interp, argv[1], &xLoc) != TCL_OK) {
		opserr << "WARNING invalid x_loc - vrp x_loc y_loc z_loc\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[2], &yLoc) != TCL_OK) {
		opserr << "WARNING invalid y_loc - vrp x_loc y_loc z_loc\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[3], &zLoc) != TCL_OK) {
		opserr << "WARNING invalid z_loc - vrp x_loc y_loc z_loc\n";
		return TCL_ERROR;
	}

	OPS_GetFeViewer()->setVRP(float(xLoc), float(yLoc), float(zLoc));
	return TCL_OK;
#endif
}

int
FeViewer_setVPN(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return TCL_OK;
#else
	// check destructor has not been called
	if (OPS_GetFeViewer() == 0)
		return TCL_OK;

	// make sure at least one other argument to contain type of system
	if (argc < 4) {
		opserr << "WARNING args incorrect - vpn xdirn ydirn zdirn \n";
		return TCL_ERROR;
	}

	// get the id, x_dirn and y_dirn
	double xDirn, yDirn, zDirn;
	if (Tcl_GetDouble(interp, argv[1], &xDirn) != TCL_OK) {
		opserr << "WARNING invalid x_dirn - vpn x_dirn y_dirn z_dirn\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[2], &yDirn) != TCL_OK) {
		opserr << "WARNING invalid y_dirn - vpn x_dirn y_dirn z_dirn\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[3], &zDirn) != TCL_OK) {
		opserr << "WARNING invalid z_dirn - vpn x_dirn y_dirn z_dirn\n";
		return TCL_ERROR;
	}

	OPS_GetFeViewer()->setVPN(float(xDirn), float(yDirn), float(zDirn));
	return 0;
#endif
}

int
FeViewer_setVUP(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return TCL_OK;
#else
	// check destructor has not been called
	if (OPS_GetFeViewer() == 0)
		return TCL_OK;

	// make sure at least one other argument to contain type of system
	if (argc < 4) {
		opserr << "WARNING args incorrect - vup xdirn ydirn zdirn \n";
		return TCL_ERROR;
	}

	// get the id, x_dirn and y_dirn
	double xDirn, yDirn, zDirn;
	if (Tcl_GetDouble(interp, argv[1], &xDirn) != TCL_OK) {
		opserr << "WARNING invalid x_dirn - vup x_dirn y_dirn z_dirn\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[2], &yDirn) != TCL_OK) {
		opserr << "WARNING invalid y_dirn - vup x_dirn y_dirn z_dirn\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[3], &zDirn) != TCL_OK) {
		opserr << "WARNING invalid z_dirn - vup x_dirn y_dirn z_dirn\n";
		return TCL_ERROR;
	}

	OPS_GetFeViewer()->setVUP(float(xDirn), float(yDirn), float(zDirn));
	return TCL_OK;
#endif
}

int
FeViewer_setViewWindow(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return TCL_OK;
#else
	// check destructor has not been called
	if (OPS_GetFeViewer() == 0)
		return TCL_OK;

	// check num args
	if (argc < 5) {
		opserr << "WARNING args incorrect - prp uMin uMax vMin vMax \n";
		return TCL_ERROR;
	}

	double uMin, uMax, vMin, vMax;
	if (Tcl_GetDouble(interp, argv[1], &uMin) != TCL_OK) {
		opserr << "WARNING invalid uMin - vup uMin uMax vMin vMax\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[2], &uMax) != TCL_OK) {
		opserr << "WARNING invalid uMax - vup uMin uMax vMin vMax\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[3], &vMin) != TCL_OK) {
		opserr << "WARNING invalid vMin - vup uMin uMax vMin vMax\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[4], &vMax) != TCL_OK) {
		opserr << "WARNING invalid vMin - vup uMin uMax vMin vMax\n";
		return TCL_ERROR;
	}

	OPS_GetFeViewer()->setViewWindow(float(uMin), float(uMax), float(vMin), float(vMax));
	return TCL_OK;
#endif
}
int
FeViewer_setPlaneDist(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return TCL_OK;
#else

	// check destructor has not been called
	if (OPS_GetFeViewer() == 0)
		return TCL_OK;

	// check num args
	if (argc < 3) {
		opserr << "WARNING args incorrect - dist near far \n";
		return TCL_ERROR;
	}
	// get distances to near view and far planes
	double anear, afar;
	if (Tcl_GetDouble(interp, argv[1], &anear) != TCL_OK) {
		opserr << "WARNING invalid near - vup near far\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[2], &afar) != TCL_OK) {
		opserr << "WARNING invalid far - vup near far\n";
		return TCL_ERROR;
	}

	OPS_GetFeViewer()->setPlaneDist(float(anear), float(afar));
	return TCL_OK;
#endif
}

int
FeViewer_setProjectionMode(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return TCL_OK;
#else

	// check destructor has not been called
	if (OPS_GetFeViewer() == 0)
		return TCL_OK;

	// ensure corrcet num arguments
	if (argc < 2) {
		opserr << "WARNING args incorrect - projection modeID \n";
		return TCL_ERROR;
	}

	// set the mode
	OPS_GetFeViewer()->setProjectionMode(argv[1]);
	return TCL_OK;
#endif
}

int
FeViewer_setFillMode(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return TCL_OK;
#else

	// check destructor has not been called
	if (OPS_GetFeViewer() == 0)
		return TCL_OK;

	// ensure corrcet num arguments
	if (argc < 2) {
		opserr << "WARNING args incorrect - fill modeID \n";
		return TCL_ERROR;
	}

	// set the mode
	OPS_GetFeViewer()->setFillMode(argv[1]);
	return TCL_OK;
#endif
}

int
FeViewer_setPRP(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return TCL_OK;
#else

	// check destructor has not been called
	if (OPS_GetFeViewer() == 0)
		return TCL_OK;

	// ensure corrcet num arguments
	if (argc < 4) {
		opserr << "WARNING args incorrect - cop xloc yloc zloc \n";
		return TCL_ERROR;
	}

	// get the xLoc, yLoc and zLoc
	double xLoc, yLoc, zLoc;
	if (Tcl_GetDouble(interp, argv[1], &xLoc) != TCL_OK) {
		opserr << "WARNING invalid x_loc - cop x_loc y_loc z_loc\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[2], &yLoc) != TCL_OK) {
		opserr << "WARNING invalid y_loc - cop x_loc y_loc z_loc\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[3], &zLoc) != TCL_OK) {
		opserr << "WARNING invalid z_loc - cop x_loc y_loc z_loc\n";
		return TCL_ERROR;
	}

	OPS_GetFeViewer()->setPRP(float(xLoc), float(yLoc), float(zLoc));
	return TCL_OK;
#endif
}

int
FeViewer_setPortWindow(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return TCL_OK;
#else

	// check destructor has not been called
	if (OPS_GetFeViewer() == 0)
		return TCL_OK;

	// check num args  
	if (argc < 5) {
		opserr << "WARNING args incorrect - prp uMin uMax vMin vMax \n";
		return TCL_ERROR;
	}

	double uMin, uMax, vMin, vMax;
	if (Tcl_GetDouble(interp, argv[1], &uMin) != TCL_OK) {
		opserr << "WARNING invalid uMin - vup uMin uMax vMin vMax\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[2], &uMax) != TCL_OK) {
		opserr << "WARNING invalid uMax - vup uMin uMax vMin vMax\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[3], &vMin) != TCL_OK) {
		opserr << "WARNING invalid vMin - vup uMin uMax vMin vMax\n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[4], &vMax) != TCL_OK) {
		opserr << "WARNING invalid vMin - vup uMin uMax vMin vMax\n";
		return TCL_ERROR;
	}

	OPS_GetFeViewer()->setPortWindow(float(uMin), float(uMax), float(vMin), float(vMax));
	return TCL_OK;
#endif
}

int
FeViewer_displayModel(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return TCL_OK;
#else

	// check destructor has not been called
	if (OPS_GetFeViewer() == 0)
		return TCL_OK;

	// check number of args  
	if (argc < 3 || argc > 5) {
		opserr << "WARNING args incorrect - display eleMode <nodeMode> displayFact <lineWidth>\n";
		return TCL_ERROR;
	}

	if (argc == 3) {
		int displayMode;
		double displayFact;
		if (Tcl_GetInt(interp, argv[1], &displayMode) != TCL_OK) {
			opserr << "WARNING invalid displayMode - display displayMode displayFact\n";
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[2], &displayFact) != TCL_OK) {
			opserr << "WARNING invalid displayFact - display displayMode displayFact\n";
			return TCL_ERROR;
		}

		OPS_GetFeViewer()->displayModel(displayMode, -1, float(displayFact));
		return TCL_OK;
	}
	else {
		int eleMode, nodeMode;
		double displayFact;
		if (Tcl_GetInt(interp, argv[1], &eleMode) != TCL_OK) {
			opserr << "WARNING invalid eleMode - display eleMode nodeMode displayFact\n";
			return TCL_ERROR;
		}
		if (Tcl_GetInt(interp, argv[2], &nodeMode) != TCL_OK) {
			opserr << "WARNING invalid nodeMode - display eleMode nodeMode displayFact\n";
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[3], &displayFact) != TCL_OK) {
			opserr << "WARNING invalid displayFact - display eleMode nodeMode displayFact\n";
			return TCL_ERROR;
		}
		// line width
		if (argc == 5) {
			int lineWidth;
			if (Tcl_GetInt(interp, argv[4], &lineWidth) != TCL_OK) {
				opserr << "WARNING invalid lineWidth - display eleMode nodeMode displayFact lineWidth\n";
				return TCL_ERROR;
			}
			OPS_GetFeViewer()->displayModel(eleMode, nodeMode, float(displayFact), lineWidth);
			return TCL_OK;
		}
		OPS_GetFeViewer()->displayModel(eleMode, nodeMode, float(displayFact));
		return TCL_OK;
	}
#endif
}

int
FeViewer_clearImage(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return TCL_OK;
#else

	// check destructor has not been called
	if (OPS_GetFeViewer() == 0)
		return TCL_OK;

	OPS_GetFeViewer()->clearImage();
	return TCL_OK;
#endif
}

int
FeViewer_saveImage(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return TCL_OK;
#else

	// check destructor has not been called
	if (OPS_GetFeViewer() == 0)
		return TCL_OK;

	// check number of args  
	if (argc != 3 && argc != 5) {
		opserr << "WARNING args incorrect - saveImage -image imageName? -file fileName?\n";
		return TCL_ERROR;
	}

	int loc = 1;
	const char* imageName = 0;
	const char* fileName = 0;
	while (loc < argc) {
		if (strcmp(argv[loc], "-image") == 0)
			imageName = argv[loc + 1];
		else if (strcmp(argv[loc], "-file") == 0)
			fileName = argv[loc + 1];
		else {
			opserr << "WARNING invalid option: " << argv[loc] << " - saveImage -image imageName? -file fileName?\n";
			return TCL_ERROR;
		}
		loc += 2;
	}

	if (imageName == 0 && fileName != 0)
		OPS_GetFeViewer()->saveImage(fileName);
	else if (imageName != 0 && fileName != 0)
		OPS_GetFeViewer()->saveImage(imageName, fileName);
	else {
		opserr << "WARNING saveImage - need a fileName\n";
		return TCL_ERROR;
	}

	return TCL_OK;
#endif
}

int
TclCreateRecorder(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv, Domain& theDomain, Recorder** theRecorder)
{
	// make sure at least one other argument to contain integrator
	if (argc < 2) {
		opserr << "WARNING need to specify a Recorder type\n";
		return TCL_ERROR;
	}
	FE_Datastore* theRecorderDatabase = 0;
	OPS_Stream* theOutputStream = 0;

	TCL_Char* fileName = 0;
	TCL_Char* tableName = 0;
	// an Element Recorder or ElementEnvelope Recorder
	if ((strcmp(argv[1], "Element") == 0) || (strcmp(argv[1], "EnvelopeElement") == 0)
		|| (strcmp(argv[1], "NormElement") == 0) || (strcmp(argv[1], "NormEnvelopeElement") == 0)
		|| strcmp(argv[1], "ResidElement") == 0 || strcmp(argv[1], "ConditionalElement") == 0)
	{
		(*theRecorder) = (Recorder*)OPS_ElementRecorder(argv[1]);
	}

	//else if ((strcmp(argv[1], "Damage") == 0) || (strcmp(argv[1], "ElementDamage") == 0) ||
	//	(strcmp(argv[1], "damage") == 0) || (strcmp(argv[1], "elementDamage") == 0)) {


	//	if (argc < 7) {
	//		opserr << "WARNING recorder ElementDamage eleID? <-time> "
	//			<< "<-file fileName?> <-section secID1? secID2? ...> <-dof dofID?> <-damage dmgID?>";
	//		return TCL_ERROR;
	//	}

	//	double dT = 0.0;
	//	double rTolDt = 0.00001;
	//	bool echoTime = false;
	//	int loc = 2;
	//	int eleID;

	//	if (Tcl_GetInt(interp, argv[loc], &eleID) != TCL_OK) {
	//		opserr << "WARNING recorder ElementDamage: No element tag specified ";
	//		return TCL_ERROR;
	//	}
	//	loc++;

	//	if ((strcmp(argv[loc], "-time") == 0) || (strcmp(argv[loc], "-load") == 0)) {
	//		// allow user to specify const load
	//		echoTime = true;
	//		loc++;
	//	}
	//	else if (strcmp(argv[loc], "-dT") == 0) {
	//		// allow user to specify time step size for recording
	//		loc++;
	//		if (Tcl_GetDouble(interp, argv[loc], &dT) != TCL_OK)
	//			return TCL_ERROR;
	//		loc++;
	//	}
	//	else if (strcmp(argv[loc], "-rTolDt") == 0) {
	//		// allow user to specify time step tolerance for recording
	//		loc++;
	//		if (Tcl_GetDouble(interp, argv[loc], &rTolDt) != TCL_OK)
	//			return TCL_ERROR;
	//		loc++;
	//	}


	//	if (strcmp(argv[loc], "-file") == 0) {
	//		// allow user to specify load pattern other than current
	//		loc++;
	//		fileName = argv[loc];
	//		loc++;
	//	}

	//	if (strcmp(argv[loc], "-section") != 0 && strcmp(argv[loc], "section") != 0) {
	//		opserr << "WARNING recorder ElementDamage: Section keyword not specified ";
	//		return TCL_ERROR;
	//	}
	//	loc++;

	//	int secID;
	//	int endSecIDs = loc;
	//	int numSec = 0;
	//	while (Tcl_GetInt(interp, argv[endSecIDs], &secID) == TCL_OK) {
	//		endSecIDs++;
	//	}

	//	numSec = endSecIDs - loc;
	//	// create an ID to hold section/material tags
	//	ID secIDs(numSec);

	//	// read in the sec tags to the ID
	//	for (int i = loc; i < endSecIDs; i++) {
	//		if (Tcl_GetInt(interp, argv[i], &secID) != TCL_OK)
	//			return TCL_ERROR;
	//		secIDs[loc - i] = secID;
	//	}

	//	loc = endSecIDs;

	//	int dofID = 0;
	//	if (strcmp(argv[loc], "-dof") == 0 || strcmp(argv[loc], "dof") == 0 ||
	//		strcmp(argv[loc], "-DOF") == 0 || strcmp(argv[loc], "DOF") == 0) {
	//		loc++;
	//		if (Tcl_GetInt(interp, argv[loc], &dofID) != TCL_OK) {
	//			opserr << "WARNING recorder ElementDamage: No dof tag specified ";
	//			return TCL_ERROR;
	//		}
	//		loc++;
	//	}

	//	if (strcmp(argv[loc], "-damage") != 0 && strcmp(argv[loc], "damage") != 0) {
	//		opserr << "WARNING recorder ElementDamage: No damege tag specified ";
	//		return TCL_ERROR;
	//	}
	//	loc++;

	//	int dmgID;
	//	if (Tcl_GetInt(interp, argv[loc], &dmgID) != TCL_OK) {
	//		opserr << "WARNING recorder ElementDamage: No damege tag specified ";
	//		return TCL_ERROR;
	//	}

	//	DamageModel* dmgPTR;
	//	dmgPTR = OPS_getDamageModel(dmgID);

	//	if (dmgPTR == NULL)
	//	{
	//		opserr << "WARNING recorder ElementDamage: specified damage model not found\n";
	//		return TCL_ERROR;
	//	}


	//	//	const char **data = new const char *[argc-eleData];

	//	OPS_Stream* theOutput = new DataFileStream(fileName);

	//	// now construct the recorder
	//	(*theRecorder) = new DamageRecorder(eleID, secIDs, dofID, dmgPTR, theDomain, echoTime, dT, rTolDt, *theOutput);

	//}
	//////////////////////End of ElementDamage recorder////////////////////////////

	//else if ((strcmp(argv[1], "Remove") == 0) || (strcmp(argv[1], "ElementRemoval") == 0) ||
	//	(strcmp(argv[1], "NodeRemoval") == 0) || (strcmp(argv[1], "Collapse") == 0)) {


	//	if (argc < 4) {
	//		opserr << "WARNING recorder Collapse -ele eleID <eleID2? ...>  -node nodeID <-time> <-file fileName?> ? "
	//			<< "\n or recorder Collapse -ele eleID1 <eleID2? ...>? <-sec secID1? secID2? ...> -crit crit1? value1?"
	//			<< " <-crit crit2? value2?> <-time> <-file fileName?> <-mass mass1? mass2? ...> <-g gAcc gDir? gPat?>?" << endln;
	//		return TCL_ERROR;
	//	}


	//	int maxNumCriteria = 2;	// current maximum number of either-or criteria for an element
	//	double dT = 0.0;
	//	double rTolDt = 0.00001;
	//	bool echoTime = false;
	//	int endEleIDs = 2;
	//	int numEle = endEleIDs - 2;
	//	int loc = endEleIDs;
	//	int flags = 0;
	//	int eleData = 0;
	//	int nodeTag = 0;
	//	//  new
	//	int nTagbotn = 0;
	//	int nTagmidn = 0;
	//	int nTagtopn = 0;
	//	int globgrav = 0;
	//	const char* fileNameinf = 0;

	//	//  end of new
	//	int numSecondaryEle = 0;

	//	// create an ID to hold ele tags
	//	//ID eleIDs(numEle, numEle+1); 
	//	ID eleIDs;
	//	eleIDs = ID(1);
	//	ID secondaryEleIDs = ID(1);
	//	secondaryEleIDs[0] = 0;
	//	bool secondaryFlag = false;
	//	ID secIDs = 0;
	//	//	secIDs = 0;

	//	// optional mass and weight definition
	//	Vector eleMass(1);
	//	eleMass.Zero();
	//	double gAcc = 0;
	//	int gDir = 0, gPat = 0;


	//	Vector remCriteria(2 * maxNumCriteria);
	//	remCriteria.Zero();
	//	int numCrit = 0;

	//	//   DataOutputHandler *theDataOutputHandler = 0;
	//	//    TCL_Char *fileName = 0;


	//	while (flags == 0 && loc < argc) {

	//		if (strcmp(argv[loc], "-node") == 0) {

	//			if (Tcl_GetInt(interp, argv[loc + 1], &nodeTag) != TCL_OK) {
	//				opserr << "WARNING recorder Collapse -node - invalid node tag " << argv[loc + 1] << endln;
	//				return TCL_ERROR;
	//			}

	//			Node* theNode = theDomain.getNode(nodeTag);
	//			if (theNode == 0) {
	//				opserr << "WARNING recorder Collapse -node - invalid node " << argv[loc + 1] << endln;
	//				return TCL_ERROR;
	//			}
	//			loc += 2;
	//		}
	//		// new

	//		else if (strcmp(argv[loc], "-file_infill") == 0) {
	//			fileNameinf = argv[loc + 1];
	//			loc += 2;
	//		}

	//		else if (strcmp(argv[loc], "-checknodes") == 0) {
	//			if (Tcl_GetInt(interp, argv[loc + 1], &nTagbotn) != TCL_OK) {
	//				opserr << "WARNING recorder Collapse -node - invalid node tag " << argv[loc + 1] << endln;
	//				return TCL_ERROR;
	//			}
	//			if (Tcl_GetInt(interp, argv[loc + 2], &nTagmidn) != TCL_OK) {
	//				opserr << "WARNING recorder Collapse -node - invalid node tag " << argv[loc + 1] << endln;
	//				return TCL_ERROR;
	//			}
	//			if (Tcl_GetInt(interp, argv[loc + 3], &nTagtopn) != TCL_OK) {
	//				opserr << "WARNING recorder Collapse -node - invalid node tag " << argv[loc + 1] << endln;
	//				return TCL_ERROR;
	//			}
	//			loc += 4;
	//		}

	//		else if (strcmp(argv[loc], "-global_gravaxis") == 0) {
	//			if (Tcl_GetInt(interp, argv[loc + 1], &globgrav) != TCL_OK) {
	//				opserr << "WARNING recorder Collapse -global_gravaxis - invalid global axis for gravity " << argv[loc + 1] << endln;
	//				return TCL_ERROR;
	//			}
	//			loc += 2;
	//		}


	//		//    end of new

	//		else if ((strcmp(argv[loc], "-slave") == 0) || (strcmp(argv[loc], "-secondary") == 0)) {
	//			secondaryFlag = true;
	//			loc++;
	//		}

	//		else if ((strcmp(argv[loc], "-ele") == 0) ||
	//			(strcmp(argv[loc], "-eles") == 0) ||
	//			(strcmp(argv[loc], "-element") == 0)) {

	//			// ensure no segmentation fault if user messes up
	//			if (argc < loc + 2) {
	//				opserr << "WARNING recorder Collapse .. -ele tag1? .. - no ele tags specified\n";
	//				return TCL_ERROR;
	//			}

	//			//
	//			// read in a list of ele until end of command or other flag
	//			//
	//			loc++;
	//			int eleTag;
	//			while (loc < argc && Tcl_GetInt(interp, argv[loc], &eleTag) == TCL_OK) {
	//				if (secondaryFlag == false)
	//					eleIDs[numEle++] = eleTag;
	//				else
	//					secondaryEleIDs[numSecondaryEle++] = eleTag;
	//				secondaryFlag = false;
	//				loc++;
	//			}

	//			Tcl_ResetResult(interp);

	//			//	    if (loc == argc) {
	//			//	      opserr << "ERROR: No response type specified for element recorder. " << endln;
	//			//	      return TCL_ERROR;
	//			//	    }

	//			if (strcmp(argv[loc], "all") == 0) {
	//				ElementIter& theEleIter = theDomain.getElements();
	//				Element* theEle;
	//				while ((theEle = theEleIter()) != 0)
	//					eleIDs[numEle++] = theEle->getTag();
	//				loc++;
	//			}

	//		}

	//		else if (strcmp(argv[loc], "-eleRange") == 0) {

	//			// ensure no segmentation fault if user messes up
	//			if (argc < loc + 3) {
	//				opserr << "WARNING recorder Element .. -eleRange start? end?  .. - no ele tags specified\n";
	//				return TCL_ERROR;
	//			}

	//			//
	//			// read in start and end tags of two elements & add set [start,end]
	//			//

	//			int start, end;
	//			if (Tcl_GetInt(interp, argv[loc + 1], &start) != TCL_OK) {
	//				opserr << "WARNING recorder Element -eleRange start? end? - invalid start " << argv[loc + 1] << endln;
	//				return TCL_ERROR;
	//			}
	//			if (Tcl_GetInt(interp, argv[loc + 2], &end) != TCL_OK) {
	//				opserr << "WARNING recorder Element -eleRange start? end? - invalid end " << argv[loc + 2] << endln;
	//				return TCL_ERROR;
	//			}
	//			if (start > end) {
	//				int swap = end;
	//				end = start;
	//				start = swap;
	//			}

	//			for (int i = start; i <= end; i++)
	//				if (secondaryFlag == false)
	//					eleIDs[numEle++] = i;
	//				else
	//					secondaryEleIDs[numSecondaryEle++] = i;

	//			secondaryFlag = false;
	//			loc += 3;
	//		}

	//		else if (strcmp(argv[loc], "-region") == 0) {
	//			// allow user to specif elements via a region

	//			if (argc < loc + 2) {
	//				opserr << "WARNING recorder Element .. -region tag?  .. - no region specified\n";
	//				return TCL_ERROR;
	//			}
	//			int tag;
	//			if (Tcl_GetInt(interp, argv[loc + 1], &tag) != TCL_OK) {
	//				opserr << "WARNING recorder Element -region tag? - invalid tag " << argv[loc + 1] << endln;
	//				return TCL_ERROR;
	//			}
	//			MeshRegion* theRegion = theDomain.getRegion(tag);
	//			if (theRegion == 0) {
	//				opserr << "WARNING recorder Element -region " << tag << " - region does not exist" << endln;
	//				return TCL_OK;
	//			}
	//			const ID& eleRegion = theRegion->getElements();
	//			for (int i = 0; i < eleRegion.Size(); i++)
	//				if (secondaryFlag == false)
	//					eleIDs[numEle++] = eleRegion(i);
	//				else
	//					secondaryEleIDs[numSecondaryEle++] = eleRegion(i);

	//			secondaryFlag = false;
	//			loc += 2;
	//		}

	//		else if ((strcmp(argv[loc], "-time") == 0) || (strcmp(argv[loc], "-load") == 0)) {
	//			// allow user to specify const load
	//			echoTime = true;
	//			loc++;
	//		}

	//		else if (strcmp(argv[loc], "-dT") == 0) {
	//			// allow user to specify time step size for recording
	//			loc++;
	//			if (Tcl_GetDouble(interp, argv[loc], &dT) != TCL_OK)
	//				return TCL_ERROR;
	//			loc++;
	//		}
	//		else if (strcmp(argv[loc], "-rTolDt") == 0) {
	//			// allow user to specify time step tolerance for recording
	//			loc++;
	//			if (Tcl_GetDouble(interp, argv[loc], &rTolDt) != TCL_OK)
	//				return TCL_ERROR;
	//			loc++;
	//		}

	//		else if (strcmp(argv[loc], "-file") == 0) {
	//			fileName = argv[loc + 1];
	//			//	    simulationInfo.addWriteFile(fileName);
	//			loc += 2;
	//		}

	//		else if ((strcmp(argv[loc], "-mass") == 0)) {
	//			eleMass.resize(numEle);
	//			eleMass.Zero();
	//			loc++;
	//			double eleM = 0;

	//			if (loc < argc && Tcl_GetDouble(interp, argv[loc + 1], &eleM) != TCL_OK) {
	//				Tcl_GetDouble(interp, argv[loc], &eleM);
	//				for (int i = 0; i < numEle; i++)
	//					eleMass(i) = eleM;
	//				loc++;
	//			}
	//			else {
	//				int i = 0;
	//				while (loc < argc && Tcl_GetDouble(interp, argv[loc], &eleM) == TCL_OK) {
	//					eleMass(i) = eleM;
	//					loc++;
	//					i++;
	//				}
	//			}

	//			//Tcl_ResetResult(interp);
	//		}

	//		else if ((strcmp(argv[loc], "-g") == 0)) {
	//			loc++;

	//			if (Tcl_GetDouble(interp, argv[loc], &gAcc) != TCL_OK) {
	//				opserr << "WARNING recorder Remove -g gValue? gDir? gPat?... invalid gValue ";
	//				opserr << argv[loc] << endln;
	//				return TCL_ERROR;
	//			}

	//			if (Tcl_GetInt(interp, argv[loc + 1], &gDir) != TCL_OK) {
	//				opserr << "WARNING recorder Remove -g gValue? gDir? gPat?... invalid gDir ";
	//				opserr << argv[loc + 1] << endln;
	//				return TCL_ERROR;
	//			}

	//			if (Tcl_GetInt(interp, argv[loc + 2], &gPat) != TCL_OK) {
	//				opserr << "WARNING recorder Remove -g gValue? gDir? gPat?... invalid gPat ";
	//				opserr << argv[loc + 1] << endln;
	//				return TCL_ERROR;
	//			}
	//			loc += 3;
	//		}


	//		else if (strcmp(argv[loc], "-database") == 0) {
	//			theRecorderDatabase = theDatabase;
	//			if (theRecorderDatabase != 0)
	//				tableName = argv[loc + 1];
	//			else {
	//				opserr << "WARNING recorder Node .. -database &lt;fileName&gt; - NO CURRENT DATABASE, results to File instead\n";
	//				fileName = argv[loc + 1];
	//			}

	//			loc += 2;
	//		}

	//		else if ((strcmp(argv[loc], "-section") == 0) || (strcmp(argv[loc], "-sec") == 0)
	//			|| (strcmp(argv[loc], "-comp") == 0)) {

	//			loc++;

	//			int secID;
	//			int endSecIDs = loc;
	//			int numSec = 0;
	//			while (Tcl_GetInt(interp, argv[endSecIDs], &secID) == TCL_OK) {
	//				endSecIDs++;
	//			}

	//			numSec = endSecIDs - loc;
	//			// create an ID to hold section/material tags
	//			secIDs = ID(numSec);

	//			// read in the sec tags to the ID
	//			for (int i = loc; i < endSecIDs; i++) {
	//				if (Tcl_GetInt(interp, argv[i], &secID) != TCL_OK)
	//					return TCL_ERROR;
	//				secIDs[i - loc] = secID;
	//			}


	//			loc = endSecIDs;
	//		}

	//		else if (strcmp(argv[loc], "-criteria") == 0 || strcmp(argv[loc], "-crit") == 0) {

	//			int critTag = 0;
	//			double critValue = 0;

	//			if (strcmp(argv[loc + 1], "minStrain") == 0 || strcmp(argv[loc + 1], "1") == 0)
	//				critTag = 1;
	//			else if (strcmp(argv[loc + 1], "maxStrain") == 0 || strcmp(argv[loc + 1], "2") == 0)
	//				critTag = 2;
	//			else if (strcmp(argv[loc + 1], "axialDI") == 0 || strcmp(argv[loc + 1], "3") == 0)
	//				critTag = 3;
	//			else if (strcmp(argv[loc + 1], "flexureDI") == 0 || strcmp(argv[loc + 1], "4") == 0)
	//				critTag = 4;
	//			else if (strcmp(argv[loc + 1], "axialLS") == 0 || strcmp(argv[loc + 1], "5") == 0)
	//				critTag = 5;
	//			else if (strcmp(argv[loc + 1], "shearLS") == 0 || strcmp(argv[loc + 1], "6") == 0)
	//				critTag = 6;
	//			//      new
	//			else if (strcmp(argv[loc + 1], "INFILLWALL") == 0 || strcmp(argv[loc + 1], "7") == 0)
	//				critTag = 7;
	//			//
	//			else {
	//				opserr << "Error: RemoveRecorder - Removal Criteria " << argv[loc + 1] << " not recognized" << endln;
	//				return TCL_ERROR;
	//			}
	//			//     new
	//			if (critTag != 7) {
	//				if (Tcl_GetDouble(interp, argv[loc + 2], &critValue) != TCL_OK) {
	//					opserr << "WARNING recorder Remove -crit critTag? critValue?... invalid critValue ";
	//					opserr << argv[loc + 1] << endln;
	//					return TCL_ERROR;
	//				}
	//			}

	//			remCriteria[2 * numCrit] = critTag;
	//			//      new
	//			if (critTag != 7) {
	//				remCriteria[2 * numCrit + 1] = critValue;
	//			}
	//			else {
	//				remCriteria[2 * numCrit + 1] = 100.0;
	//			}
	//			numCrit++;
	//			if (critTag != 7) {
	//				loc += 3;
	//			}
	//			else {
	//				loc += 2;
	//				secIDs = ID(1);
	//				secIDs[1] = 1;   // this is not used directly, gets rid of the "-sec" for infillwall 	 
	//			}
	//		}

	//		else {
	//			// first unknown string then is assumed to start 
	//			// element response request starts
	//			eleData = loc;
	//			flags = 1;
	//		}
	//	}

	//	// if user has specified no element tags lets assume he wants them all
	//	if (numEle == 0) {
	//		ElementIter& theEleIter = theDomain.getElements();
	//		Element* theEle;
	//		while ((theEle = theEleIter()) != 0)
	//			eleIDs[numEle++] = theEle->getTag();
	//	}

	//	theOutputStream = new DummyStream();

	//	// now construct the recorder
	//	(*theRecorder) = new RemoveRecorder(nodeTag,
	//		eleIDs,
	//		secIDs,
	//		secondaryEleIDs,
	//		remCriteria,
	//		theDomain,
	//		*theOutputStream,
	//		echoTime,
	//		dT,
	//		rTolDt,
	//		fileName,
	//		eleMass,
	//		gAcc,
	//		gDir,
	//		gPat,
	//		nTagbotn,
	//		nTagmidn,
	//		nTagtopn,
	//		globgrav,
	//		fileNameinf);
	//	//      
	//}

	//////////////////////End of Component Remove recorder////////////////////////////

	// create a recorder to write nodal displacement quantities to a file
	else if ((strcmp(argv[1], "Node") == 0) || (strcmp(argv[1], "EnvelopeNode") == 0)
		|| (strcmp(argv[1], "NodeEnvelope") == 0) || (strcmp(argv[1], "ResidNode") == 0)
		|| (strcmp(argv[1], "ConditionalNode") == 0))
	{
		(*theRecorder) = (Recorder*)OPS_NodeRecorder(argv[1]);
	}

	//else if (strcmp(argv[1], "Pattern") == 0) {
	//	if (argc < 4) {
	//		opserr << "WARNING recorder Pattern filename? <startFlag> patternTag?";
	//		return TCL_ERROR;
	//	}

	//	int flag = 0;
	//	if (strcmp(argv[3], "-time") == 0)
	//		flag = 1;
	//	if (strcmp(argv[3], "-load") == 0)
	//		flag = 2;

	//	int pos = 3;
	//	if (flag != 0) pos = 4;

	//	int patternTag;

	//	if (Tcl_GetInt(interp, argv[pos++], &patternTag) != TCL_OK)
	//		return TCL_ERROR;

	//	(*theRecorder) = new PatternRecorder(patternTag, theDomain, argv[2], 0.0, flag);
	//}

	// Create a recorder to write nodal drifts to a file
	else if ((strcmp(argv[1], "Drift") == 0) || (strcmp(argv[1], "EnvelopeDrift") == 0)
		|| (strcmp(argv[1], "ResidDrift") == 0) || (strcmp(argv[1], "ConditionalDrift") == 0))
	{
		(*theRecorder) = (Recorder*)OPS_DriftRecorder(argv[1]);
	}

	// a recorder for the graphical display of the domain
	else if (strcmp(argv[1], "display") == 0) {

		int xLoc, yLoc, width, height;
		int pos = 7;
		int wipeFlag = 0;
		double dT = 0.0;
		double rTolDt = 0.00001;
		int saveToFile = 0;

		if (argc < 7) {
			opserr << "WARNING recorder display title xLoc yLoc pixelsX pixelsY <-wipe> <-dT deltaT?> <-file fileName?>";
			return TCL_ERROR;
		}
		if (Tcl_GetInt(interp, argv[3], &xLoc) != TCL_OK)
			return TCL_ERROR;
		if (Tcl_GetInt(interp, argv[4], &yLoc) != TCL_OK)
			return TCL_ERROR;
		if (Tcl_GetInt(interp, argv[5], &width) != TCL_OK)
			return TCL_ERROR;
		if (Tcl_GetInt(interp, argv[6], &height) != TCL_OK)
			return TCL_ERROR;

		while (pos < argc) {

			// check if we are to wipe image on each redraw
			if (strcmp(argv[pos], "-wipe") == 0) {
				wipeFlag = 1;
				pos++;
			}
			// allow user to specify time step size for recording
			else if (strcmp(argv[pos], "-dT") == 0) {
				if (Tcl_GetDouble(interp, argv[pos + 1], &dT) != TCL_OK)
					return TCL_ERROR;
				pos += 2;
			}
			// allow user to specify time step size for recording
			else if (strcmp(argv[pos], "-rTolDt") == 0) {
				if (Tcl_GetDouble(interp, argv[pos + 1], &rTolDt) != TCL_OK)
					return TCL_ERROR;
				pos += 2;
			}
			// save images to file
			else if (strcmp(argv[pos], "-file") == 0) {
				fileName = argv[pos + 1];
				saveToFile = 1;
				pos += 2;
			}
		}
		if (!saveToFile)
			(*theRecorder) = new FeViewer(argv[2], xLoc, yLoc, width, height, theDomain, wipeFlag, dT, rTolDt);
		else
			(*theRecorder) = new FeViewer(argv[2], xLoc, yLoc, width, height, fileName, theDomain, dT, rTolDt);

		Tcl_CreateCommand(interp, "vrp", FeViewer_setVRP,
			(ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

		Tcl_CreateCommand(interp, "vpn", FeViewer_setVPN,
			(ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

		Tcl_CreateCommand(interp, "vup", FeViewer_setVUP,
			(ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

		Tcl_CreateCommand(interp, "viewWindow", FeViewer_setViewWindow,
			(ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

		Tcl_CreateCommand(interp, "plane", FeViewer_setPlaneDist,
			(ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

		Tcl_CreateCommand(interp, "projection", FeViewer_setProjectionMode,
			(ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

		Tcl_CreateCommand(interp, "fill", FeViewer_setFillMode,
			(ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

		Tcl_CreateCommand(interp, "prp", FeViewer_setPRP,
			(ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

		Tcl_CreateCommand(interp, "port", FeViewer_setPortWindow,
			(ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

		Tcl_CreateCommand(interp, "display", FeViewer_displayModel,
			(ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

		Tcl_CreateCommand(interp, "clearImage", FeViewer_clearImage,
			(ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

		Tcl_CreateCommand(interp, "saveImage", FeViewer_saveImage,
			(ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

	}

	else if (strcmp(argv[1], "plot") == 0) {

		int xLoc, yLoc, width, height;
		if (argc < 9) {
			opserr << "WARNING recorder plot fileName? windowTitle? xLoc yLoc xPixels yPixels -columns colX1 colY1 -columns colX2 ...";
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[4], &xLoc) != TCL_OK)
			return TCL_ERROR;
		if (Tcl_GetInt(interp, argv[5], &yLoc) != TCL_OK)
			return TCL_ERROR;
		if (Tcl_GetInt(interp, argv[6], &width) != TCL_OK)
			return TCL_ERROR;
		if (Tcl_GetInt(interp, argv[7], &height) != TCL_OK)
			return TCL_ERROR;

		int loc = 8;

		double dT = 0.0;
		double rTolDt = 0.00001;
		loc = 8;
		ID cols(0, 16);
		int numCols = 0;
		while (loc < argc) {
			if ((strcmp(argv[loc], "-columns") == 0) ||
				(strcmp(argv[loc], "-cols") == 0) ||
				(strcmp(argv[loc], "-col") == 0)) {
				if (argc < loc + 2)
					return TCL_ERROR;

				int colX, colY;
				if (Tcl_GetInt(interp, argv[loc + 1], &colX) != TCL_OK)
					return TCL_ERROR;

				if (Tcl_GetInt(interp, argv[loc + 2], &colY) != TCL_OK)
					return TCL_ERROR;

				cols[numCols++] = colX;
				cols[numCols++] = colY;
				loc += 3;
			}
			else if (strcmp(argv[loc], "-dT") == 0) {

				if (Tcl_GetDouble(interp, argv[loc + 1], &dT) != TCL_OK)
					return TCL_ERROR;
				loc += 2;
			}
			else if (strcmp(argv[loc], "-rTolDt") == 0) {

				if (Tcl_GetDouble(interp, argv[loc + 1], &rTolDt) != TCL_OK)
					return TCL_ERROR;
				loc += 2;
			}
			else
				loc++;
		}

#ifdef _NOGRAPHICS
		return TCL_OK;
#else
		FilePlotter* thePlotter = new FilePlotter(argv[2], argv[3], xLoc, yLoc, width, height, dT, rTolDt);
		(*theRecorder) = thePlotter;
		thePlotter->setCol(cols);
#endif
	}

	//	else if (strcmp(argv[1], "plotDifferent") == 0) {
	//
	//		int xLoc, yLoc, width, height;
	//		if (argc < 10) {
	//			opserr << "WARNING recorder display fileName? windowTitle? xLoc yLoc pixelsX pixelsY -columns colX1 colY1 -columns colX2 ...";
	//			return TCL_ERROR;
	//		}
	//		if (Tcl_GetInt(interp, argv[5], &xLoc) != TCL_OK)
	//			return TCL_ERROR;
	//		if (Tcl_GetInt(interp, argv[6], &yLoc) != TCL_OK)
	//			return TCL_ERROR;
	//		if (Tcl_GetInt(interp, argv[7], &width) != TCL_OK)
	//			return TCL_ERROR;
	//		if (Tcl_GetInt(interp, argv[8], &height) != TCL_OK)
	//			return TCL_ERROR;
	//
	//		int loc = 9;
	//
	//		double dT = 0.0;
	//		double rTolDt = 0.00001;
	//		loc = 0;
	//		ID cols(0, 16);
	//		int numCols = 0;
	//		while (loc < argc) {
	//			if ((strcmp(argv[loc], "-columns") == 0) ||
	//				(strcmp(argv[loc], "-cols") == 0) ||
	//				(strcmp(argv[loc], "-col") == 0)) {
	//				if (argc < loc + 2)
	//					return TCL_ERROR;
	//
	//				int colX, colY;
	//				if (Tcl_GetInt(interp, argv[loc + 1], &colX) != TCL_OK)
	//					return TCL_ERROR;
	//
	//				if (Tcl_GetInt(interp, argv[loc + 2], &colY) != TCL_OK)
	//					return TCL_ERROR;
	//
	//				cols[numCols++] = colX;
	//				cols[numCols++] = colY;
	//				loc += 3;
	//			}
	//			else if (strcmp(argv[loc], "-dT") == 0) {
	//
	//				if (Tcl_GetDouble(interp, argv[loc + 1], &dT) != TCL_OK)
	//					return TCL_ERROR;
	//				loc += 2;
	//			}
	//			else if (strcmp(argv[loc], "-rTolDt") == 0) {
	//
	//				if (Tcl_GetDouble(interp, argv[loc + 1], &rTolDt) != TCL_OK)
	//					return TCL_ERROR;
	//				loc += 2;
	//			}
	//			else
	//				loc++;
	//		}
	//
	//#ifdef _NOGRAPHICS
	//		return TCL_OK;
	//#else
	//		FilePlotter* thePlotter = new FilePlotter(argv[2], argv[3], argv[4], xLoc, yLoc, width, height, dT, rTolDt);
	//		(*theRecorder) = thePlotter;
	//		thePlotter->setCol(cols);
	//#endif
	//	}


	//	else if (strcmp(argv[1], "increments") == 0) {
	//
	//		int xLoc, yLoc, width, height;
	//
	//		if (theAlgorithm == 0) {
	//			opserr << "WARNING recorder increments - only allowed as algorithmRecorder";
	//			return TCL_ERROR;
	//		}
	//		if (argc < 7) {
	//			opserr << "WARNING recorder display windowTitle? xLoc yLoc pixelsX pixelsY ";
	//			return TCL_ERROR;
	//		}
	//		if (Tcl_GetInt(interp, argv[3], &xLoc) != TCL_OK)
	//			return TCL_ERROR;
	//		if (Tcl_GetInt(interp, argv[4], &yLoc) != TCL_OK)
	//			return TCL_ERROR;
	//		if (Tcl_GetInt(interp, argv[5], &width) != TCL_OK)
	//			return TCL_ERROR;
	//		if (Tcl_GetInt(interp, argv[6], &height) != TCL_OK)
	//			return TCL_ERROR;
	//
	//		TCL_Char* fileName = 0;
	//		bool displayRecord = false;
	//		int loc = 7;
	//		while (loc < argc) {
	//			if ((strcmp(argv[loc], "-file") == 0) ||
	//				(strcmp(argv[loc], "-file") == 0)) {
	//
	//				if (argc < loc + 1)
	//					return TCL_ERROR;
	//				loc++;
	//				fileName = argv[loc];
	//				loc++;
	//			}
	//			else if (strcmp(argv[loc], "-display") == 0) {
	//				displayRecord = true;
	//				loc++;
	//			}
	//			else
	//				loc++;
	//		}
	//
	//#ifdef _NOGRAPHICS
	//		return TCL_OK;
	//#else
	//		AlgorithmIncrements* thePlotter = new AlgorithmIncrements(theAlgorithm,
	//			argv[2], xLoc, yLoc, width, height,
	//			displayRecord, fileName);
	//		(*theRecorder) = thePlotter;
	//#endif // _NOGRAPHICS
	//	}
		//else if (strcmp(argv[1], "pvd") == 0 || strcmp(argv[1], "PVD") == 0) {
		//	OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);
		//	(*theRecorder) = (Recorder*)OPS_PVDRecorder();
		//}

		//else if (strcmp(argv[1], "vtk") == 0 || strcmp(argv[1], "VTK") == 0) {
		//	OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);
		//	(*theRecorder) = (Recorder*)OPS_VTK_Recorder();
		//}
		//else if (strcmp(argv[1], "ElementRMS") == 0) {
		//	OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);
		//	(*theRecorder) = (Recorder*)OPS_ElementRecorderRMS();
		//}
		//else if (strcmp(argv[1], "NodeRMS") == 0) {
		//	OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);
		//	(*theRecorder) = (Recorder*)OPS_NodeRecorderRMS();
		//}
	//#ifdef _HDF5
	//	else if (strcmp(argv[1], "mpco") == 0) {
	//		OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);
	//		(*theRecorder) = (Recorder*)OPS_MPCORecorder();
	//	}
	//#endif // _HDF5
		//else if (strcmp(argv[1], "gmsh") == 0 || strcmp(argv[1], "GMSH") == 0) {
		//	OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);
		//	(*theRecorder) = (Recorder*)OPS_GmshRecorder();
		//}
		// else if (strcmp(argv[1],"gmshparallel") == 0 || strcmp(argv[1],"GMSHPARALLEL") == 0) {
		//  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);
		//  (*theRecorder) = (Recorder*) OPS_GmshRecorderParallel();
		//  }

		 /* *****************************************
		 else if (strcmp(argv[1],"GSA") == 0) {
		if (argc < 3) {
		 opserr << argc;
		 opserr << "WARNING recorder GSA -file fileName? -dT deltaT? - not enough arguments\n";
		 return TCL_ERROR;
		}
		TCL_Char *fileName = 0;
		TCL_Char *title1 =0;
		TCL_Char *title2 =0;
		TCL_Char *title3 =0;
		TCL_Char *jobno =0;
		TCL_Char *initials =0;
		TCL_Char *spec =0;
		TCL_Char *currency =0;
		TCL_Char *length =0;
		TCL_Char *force =0;
		TCL_Char *temp =0;
		double dT = 0.0;
		int loc = 2;

		while (loc < argc) {
		 if ((strcmp(argv[loc],"-file") == 0) ||
				(strcmp(argv[loc],"-file") == 0)) {
			 fileName = argv[loc+1];
			 loc += 2;
		 } else if ((strcmp(argv[loc],"-title1") == 0) ||
				(strcmp(argv[loc],"-Title1e") == 0)) {
			 title1 = argv[loc+1];
			 loc += 2;
		 } else if ((strcmp(argv[loc],"-title2") == 0) ||
				(strcmp(argv[loc],"-Title2e") == 0)) {
			 title2 = argv[loc+1];
			 loc += 2;
		 } else if ((strcmp(argv[loc],"-title3") == 0) ||
				(strcmp(argv[loc],"-Title3e") == 0)) {
			 title3 = argv[loc+1];
			 loc += 2;
		 } else if ((strcmp(argv[loc],"-jobno") == 0) ||
				(strcmp(argv[loc],"-JobNo") == 0)) {
			 jobno = argv[loc+1];
			 loc += 2;
		 } else if ((strcmp(argv[loc],"-initials") == 0) ||
				(strcmp(argv[loc],"-Initials") == 0)) {
			 initials = argv[loc+1];
			 loc += 2;
		 } else if ((strcmp(argv[loc],"-spec") == 0) ||
				(strcmp(argv[loc],"-Spec") == 0)) {
			 spec = argv[loc+1];
			 loc += 2;
		 } else if ((strcmp(argv[loc],"-currency") == 0) ||
				(strcmp(argv[loc],"-Currency") == 0)) {
			 currency = argv[loc+1];
			 loc += 2;
		 } else if ((strcmp(argv[loc],"-length") == 0) ||
				(strcmp(argv[loc],"-Length") == 0)) {
			 length = argv[loc+1];
			 loc += 2;
		 } else if ((strcmp(argv[loc],"-force") == 0) ||
				(strcmp(argv[loc],"-Force") == 0)) {
			 force = argv[loc+1];
			 loc += 2;
		 } else if ((strcmp(argv[loc],"-temp") == 0) ||
				(strcmp(argv[loc],"-Temp") == 0)) {
			 temp = argv[loc+1];
			 loc += 2;
		 }
		 else if (strcmp(argv[loc],"-dT") == 0) {
			 if (Tcl_GetDouble(interp, argv[loc+1], &dT) != TCL_OK)
				return TCL_ERROR;
			 loc += 2;
		 }
		 else
			 loc++;
		}

		GSA_Recorder *theR = new GSA_Recorder(theDomain, fileName, title1, title2, title3,
								jobno, initials, spec, currency, length, force,
								temp, dT);
		(*theRecorder) = theR;
		 }
		 ************************************************* */

		 //else {

		 //	// try existing loaded packages
		 //	ExternalRecorderCommand* recorderCommands = theExternalRecorderCommands;
		 //	bool found = false;

		 //	while (recorderCommands != NULL && found == false) {
		 //		if (strcmp(argv[1], recorderCommands->funcName) == 0) {

		 //			OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);
		 //			void* theRes = (*(recorderCommands->funcPtr))();
		 //			if (theRes != 0) {
		 //				*theRecorder = (Recorder*)theRes;
		 //				found = true;
		 //			}
		 //		}
		 //		else
		 //			recorderCommands = recorderCommands->next;
		 //	}

		 //	//
		 //	// if not there try loading package
		 //	//

		 //	if (found == false) {

		 //		void* libHandle;
		 //		void* (*funcPtr)();
		 //		int recorderNameLength = strlen(argv[1]);
		 //		char* tclFuncName = new char[recorderNameLength + 5];
		 //		strcpy(tclFuncName, "OPS_");
		 //		strcpy(&tclFuncName[4], argv[1]);

		 //		int res = getLibraryFunction(argv[1], tclFuncName, &libHandle, (void**)&funcPtr);

		 //		delete[] tclFuncName;

		 //		if (res == 0) {
		 //			opserr << "commands.cpp - : " << argv[1] << " FOUND & LOADED\n";
		 //			char* recorderName = new char[recorderNameLength + 1];
		 //			strcpy(recorderName, argv[1]);
		 //			ExternalRecorderCommand* theRecorderCommand = new ExternalRecorderCommand;
		 //			theRecorderCommand->funcPtr = funcPtr;
		 //			theRecorderCommand->funcName = recorderName;
		 //			theRecorderCommand->next = theExternalRecorderCommands;
		 //			theExternalRecorderCommands = theRecorderCommand;

		 //			OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);

		 //			void* theRes = (*funcPtr)();
		 //			if (theRes != 0) {
		 //				*theRecorder = (Recorder*)theRes;
		 //			}
		 //		}
		 //		else {
		 //			opserr << "commands.cpp - : " << argv[1] << " NOT FOUND\n";
		 //		}
		 //	}

		 //	if (*theRecorder == 0) {
		 //		opserr << "WARNING No recorder type exists ";
		 //		opserr << "for recorder of type:" << argv[1];

		 //		return TCL_ERROR;
		 //	}
		 //}

		 // check we instantiated a recorder .. if not ran out of memory
	if ((*theRecorder) == 0) {
		opserr << "WARNING could not create the recorder: " << argv[1] << endln;
		return TCL_ERROR;
	}

	// operation successfull
	return TCL_OK;
}


int
TclAddRecorder(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv, Domain& theDomain)
{
	Recorder* theRecorder = 0;
	OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);
	TclCreateRecorder(clientData, interp, argc, argv, theDomain, &theRecorder);

	if (theRecorder == 0) {

		char buffer[] = "-1";
		Tcl_SetResult(interp, buffer, TCL_VOLATILE);
		//sprintf(interp->result,"-1");
		return TCL_ERROR;
	}

	if ((theDomain.addRecorder(*theRecorder)) < 0) {
		opserr << "WARNING could not add to domain - recorder " << argv[1] << endln;
		delete theRecorder;
		char buffer[] = "-1";
		Tcl_SetResult(interp, buffer, TCL_VOLATILE);
		//    sprintf(interp->result,"-1");
		return TCL_ERROR;
	}

	int recorderTag = theRecorder->getTag();
	char buffer[30];
	sprintf(buffer, "%d", recorderTag);
	Tcl_SetResult(interp, buffer, TCL_VOLATILE);
	//  sprintf(interp->result,"%d",recorderTag);

	return TCL_OK;
}


int
TclAddAlgorithmRecorder(ClientData clientData, Tcl_Interp* interp, int argc,
	TCL_Char** argv, Domain& theDomain, EquiSolnAlgo* theAlgo)
{
	Recorder* theRecorder = 0;
	theAlgorithm = theAlgo;

	TclCreateRecorder(clientData, interp, argc, argv, theDomain, &theRecorder);

	if (theRecorder == 0) {
		char buffer[] = "-1";
		Tcl_SetResult(interp, buffer, TCL_VOLATILE);
		//    sprintf(interp->result,"-1");
		return TCL_ERROR;
	}

	// add the recorder to the domain, 
	// NOTE: will not be called with theALgo == 0
	// see ~/g3/SRC/tcl/commands.C file
	if (theAlgorithm != 0) {
		if ((theAlgorithm->addRecorder(*theRecorder)) < 0) {
			opserr << "WARNING could not add to domain - recorder " << argv[1] << endln;
			delete theRecorder;
			return TCL_ERROR;
		}
	}

	int recorderTag = theRecorder->getTag();
	char buffer[30];
	sprintf(buffer, "%d", recorderTag);
	Tcl_SetResult(interp, buffer, TCL_VOLATILE);

	return TCL_OK;
}

