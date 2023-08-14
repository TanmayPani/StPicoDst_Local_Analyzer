#ifndef TTreeJet_h
#define TTreeJet_h

#include "TObject.h"

class TTreeJet : public TObject{
public:
    TTreeJet(){}
    virtual ~TTreeJet(){}

    double Pt = 0;
    double Eta = -99;
    double Phi = -99;
    double NEF = -1;
    double Area = 0;  
    double NNeutral = 0;  
    double NCharged = 0;
    double JetPtD = 0;
    double JetGirth = 0;
    double JetLeSub = 0;

    ClassDef(TTreeJet, 1)
};

#endif