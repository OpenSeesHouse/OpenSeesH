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
                                                                        
// $Revision: 1.10 $
// $Date: 2009-04-14 21:14:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/EnvelopeElementRecorder.h,v $
                                                                        
#ifndef EnvelopeElementRecorder_h
#define EnvelopeElementRecorder_h

// Written: fmk 
//
// What: "@(#) EnvelopeElementRecorder.h, revA"

#include <Recorder.h>
#include <Information.h>
#include <OPS_Globals.h>
#include <ID.h>


class Domain;
class Vector;
class Matrix;
class Element;
class Response;
class FE_Datastore;

class EnvelopeElementRecorder: public Recorder
{
  public:
    EnvelopeElementRecorder();
    EnvelopeElementRecorder(const ID *eleID, 
			    const char **argv, 
			    int argc,
			    Domain &theDomain, 
			    OPS_Stream *theOutputHandler,
          int procMethod, int procGrpNum,
			    bool echoTimeFlag,
			    const ID *dof); 


    ~EnvelopeElementRecorder();

    int record(int commitTag, double timeStamp);
    int restart(void);    
    int flush(void);    

    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
	virtual double getRecordedValue(int clmnId, int rowOffset, bool reset); //added by SAJalali
#ifdef _CSS
	virtual int removeComponentResponse(int compTag);
   int procDataMethod;  //flag to indicate element group processing method:
                        //0: no processing, print separate results
                        //1: summate results
                        //2: maximize results
                        //3: minimize results
   int procGrpNum;
   virtual int getModified() { return Modified; }
   int Modified;
#endif // _CSS

  protected:
    
  private:	
    int initialize(void);

    int numEle;
    int numDOF;
    
    ID *eleID;
    ID *dof;

    Response **theResponses;

    Domain *theDomain;
    OPS_Stream *theHandler;

    Matrix *data;
    Vector *currentData;
    bool first;

    bool initializationDone;
    char **responseArgs;
    int numArgs;

    bool echoTimeFlag; 

    int addColumnInfo;
    bool closeOnWrite;
};


#endif
