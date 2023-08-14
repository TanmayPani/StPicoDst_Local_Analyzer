#ifndef TTreeEvent_h
#define TTreeEvent_h

#include "TObject.h"

class TTreeEvent : public TObject {
public:
    TTreeEvent(){}
    virtual ~TTreeEvent(){}

    unsigned int runId = 0;
    unsigned int eventId = 0;
    double centrality = -1.0;
    double primaryVertexZ = -99;
    double genWeight = 1.0;
    double refMultWeight = 1.0;

    double nDetectorJets = 0;
    double nGenJets = 0;

    double raw_Qx_2   = 0;                  
    double raw_Qy_2   = 0;
    double raw_Qx_A_2 = 0;
    double raw_Qy_A_2 = 0;
    double raw_Qx_B_2 = 0;
    double raw_Qy_B_2 = 0;
    double raw_psi_2 = 0;
    double raw_psi_A_2 = 0;
    double raw_psi_B_2 = 0;
    double eventPlaneWeight      = 0;    
    double subEventPlaneWeight_A = 0;    
    double subEventPlaneWeight_B = 0;
    unsigned int eventPlaneMult  = 0;    
    unsigned int subEventPlaneMult_A  = 0;
    unsigned int subEventPlaneMult_B  = 0;    
};

#endif