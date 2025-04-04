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
// $Date: 2009-08-07 20:01:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/TclFourNodeQuadCommand.cpp,v $

// File: ~/element/TclFourNodeQuadCommand.C
// 
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the implementation of the TclModelBuilder_addFourNodeQuad()
// command. 
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <FourNodeQuad.h>
//#include <FourNodeQuadWithSensitivity.h>
//#include <ConstantPressureVolumeQuad.h>
//#include <EnhancedQuad.h>
//#include <NineNodeMixedQuad.h>
#include <NineNodeQuad.h>
#include <EightNodeQuad.h>
//#include <SixNodeTri.h>

#include <TclModelBuilder.h>
extern NDMaterial* OPS_getNDMaterial(int tag);

extern void printCommand(int argc, TCL_Char** argv);

/*  *****************************************************************************

		R E G U L A R    Q U A D

		***************************************************************************** */

int
TclModelBuilder_addFourNodeQuad(ClientData clientData, Tcl_Interp* interp,
	int argc,
	TCL_Char** argv,
	Domain* theTclDomain,
	TclModelBuilder* theTclBuilder)
{
	// ensure the destructor has not been called - 
	if (theTclBuilder == 0) {
		opserr << "WARNING builder has been destroyed\n";
		return TCL_ERROR;
	}
#ifndef _CSS
	if (theTclBuilder->getNDM() != 2 || theTclBuilder->getNDF() != 2) {
		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
		return TCL_ERROR;
	}
#else
	if (theTclBuilder->getNDM() != 2) {
		opserr << "WARNING -- element can only be used in 2d models\n";
		return TCL_ERROR;
	}
#endif // !_CSS

	// check the number of arguments is correct
	int argStart = 2;

	if ((argc - argStart) < 8) {
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		opserr << "Want: element FourNodeQuad eleTag? iNode? jNode? kNode? lNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
		return TCL_ERROR;
	}

	// get the id and end nodes 
	int FourNodeQuadId, iNode, jNode, kNode, lNode, matID;
	double thickness = 1.0;
	double p = 0.0;		// uniform normal traction (pressure)
	double rho = 0.0;		// mass density
	double b1 = 0.0;
	double b2 = 0.0;

	if (Tcl_GetInt(interp, argv[argStart], &FourNodeQuadId) != TCL_OK) {
		opserr << "WARNING invalid FourNodeQuad eleTag" << endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
		opserr << "WARNING invalid iNode\n";
		opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
		opserr << "WARNING invalid jNode\n";
		opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
		opserr << "WARNING invalid kNode\n";
		opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
		opserr << "WARNING invalid lNode\n";
		opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[5 + argStart], &thickness) != TCL_OK) {
		opserr << "WARNING invalid thickness\n";
		opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
		return TCL_ERROR;
	}

	TCL_Char* type = argv[6 + argStart];

	if (Tcl_GetInt(interp, argv[7 + argStart], &matID) != TCL_OK) {
		opserr << "WARNING invalid matID\n";
		opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
		return TCL_ERROR;
	}

	int dampingTag = 0;
	Damping* theDamping = 0;

	if ((argc - argStart) == 10) {
		if (strcmp(argv[8 + argStart], "-damp") == 0 && Tcl_GetInt(interp, argv[9 + argStart], &dampingTag) == TCL_OK) {
			theDamping = OPS_getDamping(dampingTag);
			if (theDamping == 0) {
				opserr << "WARNING damping not found\n";
				opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
				return TCL_ERROR;
			}
		}
	}
	else if ((argc - argStart) == 12) {
		if (Tcl_GetDouble(interp, argv[8 + argStart], &p) != TCL_OK) {
			opserr << "WARNING invalid pressure\n";
			opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[9 + argStart], &rho) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[10 + argStart], &b1) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[11 + argStart], &b2) != TCL_OK) {
			opserr << "WARNING invalid b2\n";
			opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
			return TCL_ERROR;
		}
	}
	else if ((argc - argStart) == 14) {
		int argPress = 8;
		if (strcmp(argv[8 + argStart], "-damp") == 0 && Tcl_GetInt(interp, argv[9 + argStart], &dampingTag) == TCL_OK) {
			theDamping = OPS_getDamping(dampingTag);
			if (theDamping == 0) {
				opserr << "WARNING damping not found\n";
				opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
				return TCL_ERROR;
			}
			argPress = 10;
		}
		else if (strcmp(argv[12 + argStart], "-damp") == 0 && Tcl_GetInt(interp, argv[13 + argStart], &dampingTag) == TCL_OK) {
			theDamping = OPS_getDamping(dampingTag);
			if (theDamping == 0) {
				opserr << "WARNING damping not found\n";
				opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
				return TCL_ERROR;
			}
		}
		if (Tcl_GetDouble(interp, argv[argPress + argStart], &p) != TCL_OK) {
			opserr << "WARNING invalid pressure\n";
			opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[1 + argPress + argStart], &rho) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[2 + argPress + argStart], &b1) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[3 + argPress + argStart], &b2) != TCL_OK) {
			opserr << "WARNING invalid b2\n";
			opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
			return TCL_ERROR;
		}
	}

	NDMaterial* theMaterial = OPS_getNDMaterial(matID);

	if (theMaterial == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << matID;
		opserr << "\nFourNodeQuad element: " << FourNodeQuadId << endln;
		return TCL_ERROR;
	}

	// now create the FourNodeQuad and add it to the Domain
	FourNodeQuad* theFourNodeQuad =
		new FourNodeQuad(FourNodeQuadId, iNode, jNode, kNode, lNode,
			*theMaterial, type, thickness, p, rho, b1, b2, theDamping);
	if (theFourNodeQuad == 0) {
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (theTclDomain->addElement(theFourNodeQuad) == false) {
		opserr << "WARNING could not add element to the domain\n";
		opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
		delete theFourNodeQuad;
		return TCL_ERROR;
	}

	// if get here we have successfully created the element and added it to the domain
	return TCL_OK;
}


/*  *****************************************************************************

		C O N S T A N T    P R E S S U R E    V O L U M E    Q U A D

		***************************************************************************** */


//int
//TclModelBuilder_addConstantPressureVolumeQuad(ClientData clientData, Tcl_Interp* interp,
//	int argc,
//	TCL_Char** argv,
//	Domain* theTclDomain,
//	TclModelBuilder* theTclBuilder)
//{
//	// ensure the destructor has not been called - 
//	if (theTclBuilder == 0) {
//		opserr << "WARNING builder has been destroyed\n";
//		return TCL_ERROR;
//	}
//
//	if (theTclBuilder->getNDM() != 2 || theTclBuilder->getNDF() != 2) {
//		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
//		return TCL_ERROR;
//	}
//
//	// check the number of arguments is correct
//	int argStart = 2;
//
//	if ((argc - argStart) < 7) {
//		opserr << "WARNING insufficient arguments\n";
//		printCommand(argc, argv);
//		opserr << "Want: element ConstantPressureVolumeQuad eleTag? iNode? jNode? kNode? lNode? thk? matTag?\n";
//		return TCL_ERROR;
//	}
//
//	// get the id and end nodes 
//	int ConstantPressureVolumeQuadId, iNode, jNode, kNode, lNode, matID;
//	double thickness = 1.0;
//
//	if (Tcl_GetInt(interp, argv[argStart], &ConstantPressureVolumeQuadId) != TCL_OK) {
//		opserr << "WARNING invalid ConstantPressureVolumeQuad eleTag" << endln;
//		return TCL_ERROR;
//	}
//	if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
//		opserr << "WARNING invalid iNode\n";
//		opserr << "ConstantPressureVolumeQuad element: " << ConstantPressureVolumeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
//		opserr << "WARNING invalid jNode\n";
//		opserr << "ConstantPressureVolumeQuad element: " << ConstantPressureVolumeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
//		opserr << "WARNING invalid kNode\n";
//		opserr << "ConstantPressureVolumeQuad element: " << ConstantPressureVolumeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
//		opserr << "WARNING invalid lNode\n";
//		opserr << "ConstantPressureVolumeQuad element: " << ConstantPressureVolumeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetDouble(interp, argv[5 + argStart], &thickness) != TCL_OK) {
//		opserr << "WARNING invalid thickness\n";
//		opserr << "ConstantPressureVolumeQuad element: " << ConstantPressureVolumeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[6 + argStart], &matID) != TCL_OK) {
//		opserr << "WARNING invalid matID\n";
//		opserr << "ConstantPressureVolumeQuad element: " << ConstantPressureVolumeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	NDMaterial* theMaterial = OPS_getNDMaterial(matID);
//
//	if (theMaterial == 0) {
//		opserr << "WARNING material not found\n";
//		opserr << "Material: " << matID;
//		opserr << "\nConstantPressureVolumeQuad element: " << ConstantPressureVolumeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	// now create the ConstantPressureVolumeQuad and add it to the Domain
//	ConstantPressureVolumeQuad* theConstantPressureVolumeQuad =
//		new ConstantPressureVolumeQuad(ConstantPressureVolumeQuadId, iNode, jNode, kNode, lNode,
//			*theMaterial, thickness);
//	if (theConstantPressureVolumeQuad == 0) {
//		opserr << "WARNING ran out of memory creating element\n";
//		opserr << "ConstantPressureVolumeQuad element: " << ConstantPressureVolumeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (theTclDomain->addElement(theConstantPressureVolumeQuad) == false) {
//		opserr << "WARNING could not add element to the domain\n";
//		opserr << "ConstantPressureVolumeQuad element: " << ConstantPressureVolumeQuadId << endln;
//		delete theConstantPressureVolumeQuad;
//		return TCL_ERROR;
//	}
//
//	// if get here we have successfully created the element and added it to the domain
//	return TCL_OK;
//}


/*  *****************************************************************************

		E N H A N C E D    Q U A D

		***************************************************************************** */
//int
//TclModelBuilder_addEnhancedQuad(ClientData clientData, Tcl_Interp* interp,
//	int argc,
//	TCL_Char** argv,
//	Domain* theTclDomain,
//	TclModelBuilder* theTclBuilder)
//{
//	// ensure the destructor has not been called - 
//	if (theTclBuilder == 0) {
//		opserr << "WARNING builder has been destroyed\n";
//		return TCL_ERROR;
//	}
//
//	if (theTclBuilder->getNDM() != 2 || theTclBuilder->getNDF() != 2) {
//		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
//		return TCL_ERROR;
//	}
//
//	// check the number of arguments is correct
//	int argStart = 2;
//
//	if ((argc - argStart) < 8) {
//		opserr << "WARNING insufficient arguments\n";
//		printCommand(argc, argv);
//		opserr << "Want: element EnhancedQuad eleTag? iNode? jNode? kNode? lNode? thk? type? matTag? \n";
//		return TCL_ERROR;
//	}
//
//	// get the id and end nodes 
//	int EnhancedQuadId, iNode, jNode, kNode, lNode, matID;
//	double thickness = 1.0;
//
//	if (Tcl_GetInt(interp, argv[argStart], &EnhancedQuadId) != TCL_OK) {
//		opserr << "WARNING invalid EnhancedQuad eleTag" << endln;
//		return TCL_ERROR;
//	}
//	if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
//		opserr << "WARNING invalid iNode\n";
//		opserr << "EnhancedQuad element: " << EnhancedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
//		opserr << "WARNING invalid jNode\n";
//		opserr << "EnhancedQuad element: " << EnhancedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
//		opserr << "WARNING invalid kNode\n";
//		opserr << "EnhancedQuad element: " << EnhancedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
//		opserr << "WARNING invalid lNode\n";
//		opserr << "EnhancedQuad element: " << EnhancedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetDouble(interp, argv[5 + argStart], &thickness) != TCL_OK) {
//		opserr << "WARNING invalid thickness\n";
//		opserr << "EnhancedQuad element: " << EnhancedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	TCL_Char* type = argv[6 + argStart];
//
//	if (Tcl_GetInt(interp, argv[7 + argStart], &matID) != TCL_OK) {
//		opserr << "WARNING invalid matID\n";
//		opserr << "EnhancedQuad element: " << EnhancedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	NDMaterial* theMaterial = OPS_getNDMaterial(matID);
//
//	if (theMaterial == 0) {
//		opserr << "WARNING material not found\n";
//		opserr << "Material: " << matID;
//		opserr << "\nEnhancedQuad element: " << EnhancedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	// now create the EnhancedQuad and add it to the Domain
//	EnhancedQuad* theEnhancedQuad =
//		new EnhancedQuad(EnhancedQuadId, iNode, jNode, kNode, lNode,
//			*theMaterial, type, thickness);
//	if (theEnhancedQuad == 0) {
//		opserr << "WARNING ran out of memory creating element\n";
//		opserr << "EnhancedQuad element: " << EnhancedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (theTclDomain->addElement(theEnhancedQuad) == false) {
//		opserr << "WARNING could not add element to the domain\n";
//		opserr << "EnhancedQuad element: " << EnhancedQuadId << endln;
//		delete theEnhancedQuad;
//		return TCL_ERROR;
//	}
//
//	// if get here we have successfully created the element and added it to the domain
//	return TCL_OK;
//}
//


/*  *****************************************************************************

		N I N E   N O D E   M I X E D  Q U A D

		***************************************************************************** */


//int
//TclModelBuilder_addNineNodeMixedQuad(ClientData clientData, Tcl_Interp* interp,
//	int argc,
//	TCL_Char** argv,
//	Domain* theTclDomain,
//	TclModelBuilder* theTclBuilder)
//{
//	// ensure the destructor has not been called - 
//	if (theTclBuilder == 0) {
//		opserr << "WARNING builder has been destroyed\n";
//		return TCL_ERROR;
//	}
//
//	if (theTclBuilder->getNDM() != 2 || theTclBuilder->getNDF() != 2) {
//		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
//		return TCL_ERROR;
//	}
//
//	// check the number of arguments is correct
//	int argStart = 2;
//
//	if ((argc - argStart) < 11) {
//		opserr << "WARNING insufficient arguments\n";
//		printCommand(argc, argv);
//		opserr << "Want: element NineNodeMixedQuad  eleTag?"
//			<< " iNode? jNode? kNode? lNode? mNode, nNode, pNode, qNode, centerNode "
//			<< " matTag?\n";
//		return TCL_ERROR;
//	}
//
//	// get the id and end nodes 
//	int NineNodeMixedQuadId, iNode, jNode, kNode, lNode;
//	int                  mNode, nNode, pNode, qNode;
//	int centerNode;
//	int matID;
//
//	if (Tcl_GetInt(interp, argv[argStart], &NineNodeMixedQuadId) != TCL_OK) {
//		opserr << "WARNING invalid NineNodeMixedQuad eleTag" << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
//		opserr << "WARNING invalid iNode\n";
//		opserr << "NineNodeMixedQuad element: " << NineNodeMixedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
//		opserr << "WARNING invalid jNode\n";
//		opserr << "NineNodeMixedQuad element: " << NineNodeMixedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
//		opserr << "WARNING invalid kNode\n";
//		opserr << "NineNodeMixedQuad element: " << NineNodeMixedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
//		opserr << "WARNING invalid lNode\n";
//		opserr << "NineNodeMixedQuad element: " << NineNodeMixedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[5 + argStart], &mNode) != TCL_OK) {
//		opserr << "WARNING invalid mNode\n";
//		opserr << "NineNodeMixedQuad element: " << NineNodeMixedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[6 + argStart], &nNode) != TCL_OK) {
//		opserr << "WARNING invalid nNode\n";
//		opserr << "NineNodeMixedQuad element: " << NineNodeMixedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[7 + argStart], &pNode) != TCL_OK) {
//		opserr << "WARNING invalid pNode\n";
//		opserr << "NineNodeMixedQuad element: " << NineNodeMixedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[8 + argStart], &qNode) != TCL_OK) {
//		opserr << "WARNING invalid qNode\n";
//		opserr << "NineNodeMixedQuad element: " << NineNodeMixedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[9 + argStart], &centerNode) != TCL_OK) {
//		opserr << "WARNING invalid centerNode\n";
//		opserr << "NineNodeMixedQuad element: " << NineNodeMixedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//
//	if (Tcl_GetInt(interp, argv[10 + argStart], &matID) != TCL_OK) {
//		opserr << "WARNING invalid matID\n";
//		opserr << "NineNodeMixedQuad element: " << NineNodeMixedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	NDMaterial* theMaterial = OPS_getNDMaterial(matID);
//
//	if (theMaterial == 0) {
//		opserr << "WARNING material not found\n";
//		opserr << "Material: " << matID;
//		opserr << "\nNineNodeMixedQuad element: " << NineNodeMixedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	// now create the NineNodeMixedQuad and add it to the Domain
//	NineNodeMixedQuad* theNineNodeMixed =
//		new NineNodeMixedQuad(NineNodeMixedQuadId, iNode, jNode, kNode, lNode,
//			mNode, nNode, pNode, qNode, centerNode, *theMaterial);
//
//	if (theNineNodeMixed == 0) {
//		opserr << "WARNING ran out of memory creating element\n";
//		opserr << "NineNodeMixedQuad element: " << NineNodeMixedQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (theTclDomain->addElement(theNineNodeMixed) == false) {
//		opserr << "WARNING could not add element to the domain\n";
//		opserr << "NineNodeMixedQuad element: " << NineNodeMixedQuadId << endln;
//		delete theNineNodeMixed;
//		return TCL_ERROR;
//	}
//
//	// if get here we have successfully created the element and added it to the domain
//	return TCL_OK;
//}


//int
//TclModelBuilder_addFourNodeQuadWithSensitivity(ClientData clientData, Tcl_Interp* interp,
//	int argc,
//	TCL_Char** argv,
//	Domain* theTclDomain,
//	TclModelBuilder* theTclBuilder)
//{
//	// ensure the destructor has not been called - 
//	if (theTclBuilder == 0) {
//		opserr << "WARNING builder has been destroyed\n";
//		return TCL_ERROR;
//	}
//
//	if (theTclBuilder->getNDM() != 2 || theTclBuilder->getNDF() != 2) {
//		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
//		return TCL_ERROR;
//	}
//
//	// check the number of arguments is correct
//	int argStart = 2;
//
//	if ((argc - argStart) < 8) {
//		opserr << "WARNING insufficient arguments\n";
//		printCommand(argc, argv);
//		opserr << "Want: element FourNodeQuad eleTag? iNode? jNode? kNode? lNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
//		return TCL_ERROR;
//	}
//
//	// get the id and end nodes 
//	int FourNodeQuadId, iNode, jNode, kNode, lNode, matID;
//	double thickness = 1.0;
//	double p = 0.0;		// uniform normal traction (pressure)
//	double r = 0.0;		// mass density
//	double b1 = 0.0;
//	double b2 = 0.0;
//
//	if (Tcl_GetInt(interp, argv[argStart], &FourNodeQuadId) != TCL_OK) {
//		opserr << "WARNING invalid FourNodeQuadWithSensitivity eleTag" << endln;
//		return TCL_ERROR;
//	}
//	if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
//		opserr << "WARNING invalid iNode\n";
//		opserr << "FourNodeQuadWithSensitivity element: " << FourNodeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
//		opserr << "WARNING invalid jNode\n";
//		opserr << "FourNodeQuadWithSensitivity element: " << FourNodeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
//		opserr << "WARNING invalid kNode\n";
//		opserr << "FourNodeQuadWithSensitivity element: " << FourNodeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
//		opserr << "WARNING invalid lNode\n";
//		opserr << "FourNodeQuadWithSensitivity element: " << FourNodeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetDouble(interp, argv[5 + argStart], &thickness) != TCL_OK) {
//		opserr << "WARNING invalid thickness\n";
//		opserr << "FourNodeQuadWithSensitivity element: " << FourNodeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	TCL_Char* type = argv[6 + argStart];
//
//	if (Tcl_GetInt(interp, argv[7 + argStart], &matID) != TCL_OK) {
//		opserr << "WARNING invalid matID\n";
//		opserr << "FourNodeQuadWithSensitivity element: " << FourNodeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if ((argc - argStart) > 11) {
//		if (Tcl_GetDouble(interp, argv[8 + argStart], &p) != TCL_OK) {
//			opserr << "WARNING invalid pressure\n";
//			opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
//			return TCL_ERROR;
//		}
//		if (Tcl_GetDouble(interp, argv[9 + argStart], &r) != TCL_OK) {
//			opserr << "WARNING invalid rho\n";
//			opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
//			return TCL_ERROR;
//		}
//		if (Tcl_GetDouble(interp, argv[10 + argStart], &b1) != TCL_OK) {
//			opserr << "WARNING invalid b1\n";
//			opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
//			return TCL_ERROR;
//		}
//		if (Tcl_GetDouble(interp, argv[11 + argStart], &b2) != TCL_OK) {
//			opserr << "WARNING invalid b2\n";
//			opserr << "FourNodeQuad element: " << FourNodeQuadId << endln;
//			return TCL_ERROR;
//		}
//
//	}
//
//	NDMaterial* theMaterial = OPS_getNDMaterial(matID);
//
//	if (theMaterial == 0) {
//		opserr << "WARNING material not found\n";
//		opserr << "Material: " << matID;
//		opserr << "\nFourNodeQuad element: " << FourNodeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	// now create the FourNodeQuad and add it to the Domain
//	FourNodeQuadWithSensitivity* theFourNodeQuadWithSensitivity =
//		new FourNodeQuadWithSensitivity(FourNodeQuadId, iNode, jNode, kNode, lNode,
//			*theMaterial, type, thickness, p, r, b1, b2);
//	if (theFourNodeQuadWithSensitivity == 0) {
//		opserr << "WARNING ran out of memory creating element\n";
//		opserr << "FourNodeQuadWithSensitivity element: " << FourNodeQuadId << endln;
//		return TCL_ERROR;
//	}
//
//	if (theTclDomain->addElement(theFourNodeQuadWithSensitivity) == false) {
//		opserr << "WARNING could not add element to the domain\n";
//		opserr << "FourNodeQuadWithSensitivity element: " << FourNodeQuadId << endln;
//		delete theFourNodeQuadWithSensitivity;
//		return TCL_ERROR;
//	}
//
//	// if get here we have successfully created the element and added it to the domain
//	return TCL_OK;
//}

// Regular nine node quad

int
TclModelBuilder_addNineNodeQuad(ClientData clientData, Tcl_Interp* interp,
	int argc,
	TCL_Char** argv,
	Domain* theTclDomain,
	TclModelBuilder* theTclBuilder)
{
	// ensure the destructor has not been called - 
	if (theTclBuilder == 0) {
		opserr << "WARNING builder has been destroyed\n";
		return TCL_ERROR;
	}

	if (theTclBuilder->getNDM() != 2 || theTclBuilder->getNDF() != 2) {
		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
		return TCL_ERROR;
	}

	// check the number of arguments is correct
	int argStart = 2;

	if ((argc - argStart) < 13) {
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		opserr << "Want: element NineNodeQuad eleTag? iNode? jNode? kNode? lNode? nNode? mNode? pNode? qNode? cNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
		return TCL_ERROR;
	}

	// get the id and end nodes
	int NineNodeQuadId, iNode, jNode, kNode, lNode;
	int nNode, mNode, pNode, qNode, cNode, matID;
	double thickness = 1.0;
	double p = 0.0;		// uniform normal traction (pressure)
	double rho = 0.0;		// mass density
	double b1 = 0.0;
	double b2 = 0.0;

	if (Tcl_GetInt(interp, argv[argStart], &NineNodeQuadId) != TCL_OK) {
		opserr << "WARNING invalid NineNodeQuad eleTag" << endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
		opserr << "WARNING invalid iNode\n";
		opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
		opserr << "WARNING invalid jNode\n";
		opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
		opserr << "WARNING invalid kNode\n";
		opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
		opserr << "WARNING invalid lNode\n";
		opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[5 + argStart], &nNode) != TCL_OK) {
		opserr << "WARNING invalid nNode\n";
		opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[6 + argStart], &mNode) != TCL_OK) {
		opserr << "WARNING invalid mNode\n";
		opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[7 + argStart], &pNode) != TCL_OK) {
		opserr << "WARNING invalid pNode\n";
		opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[8 + argStart], &qNode) != TCL_OK) {
		opserr << "WARNING invalid qNode\n";
		opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[9 + argStart], &cNode) != TCL_OK) {
		opserr << "WARNING invalid cNode\n";
		opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[10 + argStart], &thickness) != TCL_OK) {
		opserr << "WARNING invalid thickness\n";
		opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
		return TCL_ERROR;
	}

	TCL_Char* type = argv[11 + argStart];

	if (Tcl_GetInt(interp, argv[12 + argStart], &matID) != TCL_OK) {
		opserr << "WARNING invalid matID\n";
		opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
		return TCL_ERROR;
	}

	if ((argc - argStart) > 16) {
		if (Tcl_GetDouble(interp, argv[13 + argStart], &p) != TCL_OK) {
			opserr << "WARNING invalid pressure\n";
			opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[14 + argStart], &rho) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[15 + argStart], &b1) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[16 + argStart], &b2) != TCL_OK) {
			opserr << "WARNING invalid b2\n";
			opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
			return TCL_ERROR;
		}
	}

	NDMaterial* theMaterial = OPS_getNDMaterial(matID);

	if (theMaterial == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << matID;
		opserr << "\nNineNodeQuad element: " << NineNodeQuadId << endln;
		return TCL_ERROR;
	}

	// now create the NineNodeQuad and add it to the Domain
	NineNodeQuad* theNineNodeQuad =
		new NineNodeQuad(NineNodeQuadId, iNode, jNode, kNode, lNode,
			nNode, mNode, pNode, qNode, cNode,
			*theMaterial, type, thickness, p, rho, b1, b2);
	if (theNineNodeQuad == 0) {
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (theTclDomain->addElement(theNineNodeQuad) == false) {
		opserr << "WARNING could not add element to the domain\n";
		opserr << "NineNodeQuad element: " << NineNodeQuadId << endln;
		delete theNineNodeQuad;
		return TCL_ERROR;
	}

	// if get here we have successfully created the element and added it to the domain
	return TCL_OK;
}


// Regular eight node quad

int
TclModelBuilder_addEightNodeQuad(ClientData clientData, Tcl_Interp* interp,
	int argc,
	TCL_Char** argv,
	Domain* theTclDomain,
	TclModelBuilder* theTclBuilder)
{
	// ensure the destructor has not been called - 
	if (theTclBuilder == 0) {
		opserr << "WARNING builder has been destroyed\n";
		return TCL_ERROR;
	}

	if (theTclBuilder->getNDM() != 2 || theTclBuilder->getNDF() != 2) {
		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
		return TCL_ERROR;
	}

	// check the number of arguments is correct
	int argStart = 2;

	if ((argc - argStart) < 12) {
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		opserr << "Want: element EightNodeQuad eleTag? iNode? jNode? kNode? lNode? nNode? mNode? pNode? qNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
		return TCL_ERROR;
	}

	// get the id and end nodes
	int EightNodeQuadId, iNode, jNode, kNode, lNode;
	int nNode, mNode, pNode, qNode, matID;
	double thickness = 1.0;
	double p = 0.0;		// uniform normal traction (pressure)
	double rho = 0.0;		// mass density
	double b1 = 0.0;
	double b2 = 0.0;

	if (Tcl_GetInt(interp, argv[argStart], &EightNodeQuadId) != TCL_OK) {
		opserr << "WARNING invalid EightNodeQuad eleTag" << endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
		opserr << "WARNING invalid iNode\n";
		opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
		opserr << "WARNING invalid jNode\n";
		opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
		opserr << "WARNING invalid kNode\n";
		opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
		opserr << "WARNING invalid lNode\n";
		opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[5 + argStart], &nNode) != TCL_OK) {
		opserr << "WARNING invalid nNode\n";
		opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[6 + argStart], &mNode) != TCL_OK) {
		opserr << "WARNING invalid mNode\n";
		opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[7 + argStart], &pNode) != TCL_OK) {
		opserr << "WARNING invalid pNode\n";
		opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[8 + argStart], &qNode) != TCL_OK) {
		opserr << "WARNING invalid qNode\n";
		opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[9 + argStart], &thickness) != TCL_OK) {
		opserr << "WARNING invalid thickness\n";
		opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
		return TCL_ERROR;
	}

	TCL_Char* type = argv[10 + argStart];

	if (Tcl_GetInt(interp, argv[11 + argStart], &matID) != TCL_OK) {
		opserr << "WARNING invalid matID\n";
		opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
		return TCL_ERROR;
	}

	if ((argc - argStart) > 15) {
		if (Tcl_GetDouble(interp, argv[12 + argStart], &p) != TCL_OK) {
			opserr << "WARNING invalid pressure\n";
			opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[13 + argStart], &rho) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[14 + argStart], &b1) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[15 + argStart], &b2) != TCL_OK) {
			opserr << "WARNING invalid b2\n";
			opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
			return TCL_ERROR;
		}
	}

	NDMaterial* theMaterial = OPS_getNDMaterial(matID);

	if (theMaterial == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << matID;
		opserr << "\nEightNodeQuad element: " << EightNodeQuadId << endln;
		return TCL_ERROR;
	}

	// now create the EightNodeQuad and add it to the Domain
	EightNodeQuad* theEightNodeQuad =
		new EightNodeQuad(EightNodeQuadId, iNode, jNode, kNode, lNode,
			nNode, mNode, pNode, qNode,
			*theMaterial, type, thickness, p, rho, b1, b2);
	if (theEightNodeQuad == 0) {
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
		return TCL_ERROR;
	}

	if (theTclDomain->addElement(theEightNodeQuad) == false) {
		opserr << "WARNING could not add element to the domain\n";
		opserr << "EightNodeQuad element: " << EightNodeQuadId << endln;
		delete theEightNodeQuad;
		return TCL_ERROR;
	}

	// if get here we have successfully created the element and added it to the domain
	return TCL_OK;
}

// Six node tri

//int
//TclModelBuilder_addSixNodeTri(ClientData clientData, Tcl_Interp* interp,
//	int argc,
//	TCL_Char** argv,
//	Domain* theTclDomain,
//	TclModelBuilder* theTclBuilder)
//{
//	// ensure the destructor has not been called - 
//	if (theTclBuilder == 0) {
//		opserr << "WARNING builder has been destroyed\n";
//		return TCL_ERROR;
//	}
//
//	if (theTclBuilder->getNDM() != 2 || theTclBuilder->getNDF() != 2) {
//		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
//		return TCL_ERROR;
//	}
//
//	// check the number of arguments is correct
//	int argStart = 2;
//
//	if ((argc - argStart) < 10) {
//		opserr << "WARNING insufficient arguments\n";
//		printCommand(argc, argv);
//		opserr << "Want: element SixNodeTri eleTag? iNode? jNode? kNode? lNode? nNode? mNode? pNode? qNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
//		return TCL_ERROR;
//	}
//
//	// get the id and end nodes
//	int SixNodeTriId, iNode, jNode, kNode, lNode;
//	int nNode, mNode, matID;
//	double thickness = 1.0;
//	double p = 0.0;		// uniform normal traction (pressure)
//	double rho = 0.0;		// mass density
//	double b1 = 0.0;
//	double b2 = 0.0;
//
//	if (Tcl_GetInt(interp, argv[argStart], &SixNodeTriId) != TCL_OK) {
//		opserr << "WARNING invalid SixNodeTri eleTag" << endln;
//		return TCL_ERROR;
//	}
//	if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
//		opserr << "WARNING invalid iNode\n";
//		opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
//		opserr << "WARNING invalid jNode\n";
//		opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
//		opserr << "WARNING invalid kNode\n";
//		opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
//		opserr << "WARNING invalid lNode\n";
//		opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[5 + argStart], &nNode) != TCL_OK) {
//		opserr << "WARNING invalid nNode\n";
//		opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetInt(interp, argv[6 + argStart], &mNode) != TCL_OK) {
//		opserr << "WARNING invalid mNode\n";
//		opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//		return TCL_ERROR;
//	}
//
//	if (Tcl_GetDouble(interp, argv[7 + argStart], &thickness) != TCL_OK) {
//		opserr << "WARNING invalid thickness\n";
//		opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//		return TCL_ERROR;
//	}
//
//	TCL_Char* type = argv[8 + argStart];
//
//	if (Tcl_GetInt(interp, argv[9 + argStart], &matID) != TCL_OK) {
//		opserr << "WARNING invalid matID\n";
//		opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//		return TCL_ERROR;
//	}
//
//	if ((argc - argStart) > 13) {
//		if (Tcl_GetDouble(interp, argv[10 + argStart], &p) != TCL_OK) {
//			opserr << "WARNING invalid pressure\n";
//			opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//			return TCL_ERROR;
//		}
//		if (Tcl_GetDouble(interp, argv[11 + argStart], &rho) != TCL_OK) {
//			opserr << "WARNING invalid b1\n";
//			opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//			return TCL_ERROR;
//		}
//		if (Tcl_GetDouble(interp, argv[12 + argStart], &b1) != TCL_OK) {
//			opserr << "WARNING invalid b1\n";
//			opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//			return TCL_ERROR;
//		}
//		if (Tcl_GetDouble(interp, argv[13 + argStart], &b2) != TCL_OK) {
//			opserr << "WARNING invalid b2\n";
//			opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//			return TCL_ERROR;
//		}
//	}
//
//	NDMaterial* theMaterial = OPS_getNDMaterial(matID);
//
//	if (theMaterial == 0) {
//		opserr << "WARNING material not found\n";
//		opserr << "Material: " << matID;
//		opserr << "\nSixNodeTri element: " << SixNodeTriId << endln;
//		return TCL_ERROR;
//	}
//
//	// now create the SixNodeTri and add it to the Domain
//	SixNodeTri* theSixNodeTri =
//		new SixNodeTri(SixNodeTriId, iNode, jNode, kNode, lNode,
//			nNode, mNode,
//			*theMaterial, type, thickness, p, rho, b1, b2);
//	if (theSixNodeTri == 0) {
//		opserr << "WARNING ran out of memory creating element\n";
//		opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//		return TCL_ERROR;
//	}
//
//	if (theTclDomain->addElement(theSixNodeTri) == false) {
//		opserr << "WARNING could not add element to the domain\n";
//		opserr << "SixNodeTri element: " << SixNodeTriId << endln;
//		delete theSixNodeTri;
//		return TCL_ERROR;
//	}
//
//	// if get here we have successfully created the element and added it to the domain
//	return TCL_OK;
//}
