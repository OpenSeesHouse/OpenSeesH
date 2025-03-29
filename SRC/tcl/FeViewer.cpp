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

// $Revision: 1.12 $
// $Date: 2008-02-15 23:46:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/FeViewer.cpp,v $

// Written: fmk 
// Created: 04/98
// Revision: A
//
// Description: This file contains the class implementation for FeViewer
//
// What: "@(#) FeViewer.C, revA"


#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <NodeIter.h>

#include <Renderer.h>

#ifndef _NOGRAPHICS

#ifdef _WGL
#include <OpenGLRenderer.h>
#elif _GLX
#include <OpenGLRenderer.h>
#elif _AGL
#include <OpenGLRenderer.h>
#else
#include <X11Renderer.h>
#endif

#endif

#include <PlainMap.h>

#include "FeViewer.h"

//
// some static variables used in the functions
//

// 
// the functions that will be invoked by the interpreter while building the model
//

//
// the class constructor, destructor and methods


//

FeViewer::FeViewer()
	:Recorder(RECORDER_TAGS_FeViewer),
	theMap(0), theRenderer(0), theDomain(0),
	theEleMode(-1), theNodeMode(-1), theDisplayFact(1),
	deltaT(0.0), relDeltaTTol(0.00001), nextTimeStampToRecord(0.0), wipeFlag(0),
	vrpSet(0), vpwindowSet(0), clippingPlaneDistancesSet(0)
{
	theMap = 0;
}

FeViewer::FeViewer(const char* title, int xLoc, int yLoc, int width, int height,
	Domain& _theDomain, int WipeFlag, double dT, double rTolDt)
	:Recorder(RECORDER_TAGS_FeViewer),
	theMap(0), theRenderer(0), theDomain(&_theDomain),
	theEleMode(-1), theNodeMode(-1), theDisplayFact(1),
	deltaT(dT), relDeltaTTol(rTolDt), nextTimeStampToRecord(0.0), wipeFlag(WipeFlag),
	vrpSet(0), vpwindowSet(0), clippingPlaneDistancesSet(0)
{

	// set the static pointer used in the class
	theMap = new PlainMap();

#ifndef _NOGRAPHICS

#ifdef _WGL
	theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap);
#elif _GLX
	theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap);
#elif _AGL
	theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap);
#else
	theRenderer = new X11Renderer(title, xLoc, yLoc, width, height, *theMap);
#endif

#endif

	// Call Tcl_CreateCommand for class specific commands
}

FeViewer::FeViewer(const char* title, int xLoc, int yLoc, int width, int height,
	const char* fileName, Domain& _theDomain, double dT, double rTolDt)
	:Recorder(RECORDER_TAGS_FeViewer),
	theMap(0), theRenderer(0), theDomain(&_theDomain),
	theEleMode(-1), theNodeMode(-1), theDisplayFact(1),
	deltaT(dT), relDeltaTTol(rTolDt), nextTimeStampToRecord(0.0), wipeFlag(1),
	vrpSet(0), vpwindowSet(0), clippingPlaneDistancesSet(0)
{

	// set the static pointer used in the class
	theMap = new PlainMap();

#ifndef _NOGRAPHICS

#ifdef _WGL
	theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap, 0, fileName);
#elif _GLX
	theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap, fileName, 0);
#elif _AGL
	theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap, fileName, 0);
#else
	theRenderer = new X11Renderer(title, xLoc, yLoc, width, height, *theMap, fileName);
#endif

#endif

}

FeViewer::~FeViewer()
{
	// may possibly invoke Tcl_DeleteCommand() later
	// for moment just invoke destructor on map and renderer
	// and set pointer to NULL
	if (theMap != 0)
		delete theMap;

	if (theRenderer != 0)
		delete theRenderer;

}

int
FeViewer::record(int cTag, double timeStamp)
{

#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return 0;
#else

	if (theRenderer == 0)
		return 0;


	//  theRenderer->displayModel(thetheEleMode, theNodeMode, theDisplayFact);

	// loop over the elements getting each to display itself
	// using theRenderer and displayTag as arguments.
	// first clear the image
	int res = 0;
	// where relDeltaTTol is the maximum reliable ratio between analysis time step and deltaT
	// and provides tolerance for floating point precision (see floating-point-tolerance-for-recorder-time-step.md)
	if (deltaT == 0.0 || timeStamp - nextTimeStampToRecord >= -deltaT * relDeltaTTol) {

		if (deltaT != 0.0)
			nextTimeStampToRecord = timeStamp + deltaT;

		if (wipeFlag == 1)
			theRenderer->clearImage();

		// set some quantities if not set
		if (vrpSet == 0 || vpwindowSet == 0 || clippingPlaneDistancesSet == 0) {
			const Vector& theBounds = theDomain->getPhysicalBounds();
			double xAvg = (theBounds(0) + theBounds(3)) / 2.0;
			double yAvg = (theBounds(1) + theBounds(4)) / 2.0;
			double zAvg = (theBounds(2) + theBounds(5)) / 2.0;

			if (vrpSet == 0)
				this->setVRP(float(xAvg), float(yAvg), float(zAvg));

			double diff, xDiff, yDiff, zDiff;
			xDiff = (theBounds(3) - theBounds(0));
			yDiff = (theBounds(4) - theBounds(1));
			zDiff = (theBounds(5) - theBounds(2));
			diff = xDiff;
			if (yDiff > diff)
				diff = yDiff;
			if (zDiff > diff)
				diff = zDiff;

			diff *= 1.25 * 0.5;

			if (vpwindowSet == 0)
				this->setViewWindow(float(-diff), float(diff), float(-diff), float(diff));

			if (clippingPlaneDistancesSet == 0) {
				diff = sqrt(xDiff * xDiff + yDiff * yDiff + zDiff * zDiff);
				this->setPlaneDist(float(diff), float(-diff));
			}
		}

		theRenderer->startImage();

		if (theEleMode != 0) {
			ElementIter& theElements = theDomain->getElements();
			Element* theEle;
			while ((theEle = theElements()) != 0) {
				res = theEle->displaySelf(*theRenderer, theEleMode, float(theDisplayFact));
				if (res < 0) {
					opserr << "Renderer::displayModel() - Element: \n";
					opserr << theEle->getTag() << " failed to display itself\n";
				}
			}
		}

		if (theNodeMode != 0) {
			NodeIter& theNodes = theDomain->getNodes();
			Node* theNode;
			while ((theNode = theNodes()) != 0) {
				res = theNode->displaySelf(*theRenderer, theEleMode, theNodeMode, float(theDisplayFact));
				if (res < 0) {
					opserr << "Renderer::displayModel() - Node: ";
					opserr << theNode->getTag() << " failed to display itself\n";
				}
			}
		}

		// now mark the image has having been completed
		theRenderer->doneImage();
	}

	return res;
#endif
}

int
FeViewer::playback(int cTag)
{
	return 0;
}

int
FeViewer::restart(void)
{
	return 0;
}


int
FeViewer::setVRP(float xLoc, float yLoc, float zLoc)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return 0;
#else

	int ok = theRenderer->setVRP(xLoc, yLoc, zLoc);
	if (ok == 0)
		vrpSet = 1;

	return ok;
#endif
}

int
FeViewer::setVPN(float xdirn, float ydirn, float zdirn)
{
#ifdef _NOGRAPHICS

	// if no graphics .. just return 0
	return 0;

#else

	// view plane normal
	return theRenderer->setVPN(xdirn, ydirn, zdirn);

#endif
}

int
FeViewer::setVUP(float xdirn, float ydirn, float zdirn)
{
#ifdef _NOGRAPHICS

	// if no graphics .. just return 0
	return 0;

#else

	return theRenderer->setVUP(xdirn, ydirn, zdirn);

#endif
}

int
FeViewer::setViewWindow(float uMin, float uMax, float vMin, float vMax)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return 0;
#else

	int ok = theRenderer->setViewWindow(uMin, uMax, vMin, vMax);
	if (ok == 0)
		vpwindowSet = 1;

	return ok;

#endif
}

int
FeViewer::setPlaneDist(float anear, float afar)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return 0;
#else

	int ok = theRenderer->setPlaneDist(anear, afar);
	if (ok == 0)
		clippingPlaneDistancesSet = 1;
	return ok;

#endif
}

int
FeViewer::setProjectionMode(const char* mode)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return 0;
#else
	return theRenderer->setProjectionMode(mode);
#endif
}

int
FeViewer::setFillMode(const char* mode)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return 0;
#else
	return theRenderer->setFillMode(mode);
#endif
}

int
FeViewer::setPRP(float uLoc, float vLoc, float nLoc)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return 0;
#else
	return theRenderer->setPRP(uLoc, vLoc, nLoc);
#endif
}

int
FeViewer::setPortWindow(float left, float right, float bottom, float top)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return 0;
#else
	return theRenderer->setPortWindow(left, right, bottom, top);
#endif
}

int
FeViewer::displayModel(int eleFlag, int nodeFlag, float displayFact, int lineWidth)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return 0;
#else
	// methods invoked on the FE_Viewer
	theEleMode = eleFlag;
	theNodeMode = nodeFlag;
	theDisplayFact = displayFact;
	theRenderer->setLineWidth(lineWidth);
	return this->record(0, 0.0);
#endif
}

int
FeViewer::clearImage(void)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return 0;
#else
	return theRenderer->clearImage();
#endif
}

int
FeViewer::saveImage(const char* fileName)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return 0;
#else
	return theRenderer->saveImage(fileName);
#endif
}

int
FeViewer::saveImage(const char* imageName, const char* fileName)
{
#ifdef _NOGRAPHICS
	// if no graphics .. just return 0
	return 0;
#else
	return theRenderer->saveImage(imageName, fileName);
#endif
}

int
FeViewer::sendSelf(int commitTag, Channel& theChannel)
{
	return 0;
}

int
FeViewer::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	return 0;
}
