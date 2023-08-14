#ifndef BEMCLocator_H
#define BEMCLocator_H

#include "TMath.h"
#include "TVector3.h"

#include <iostream>
#include <cassert>
#include <vector>

class BEMCLocator {
public:
    BEMCLocator();
    virtual ~BEMCLocator() {}

    void getBin(const int softId, int &m, int &e, int &s) const;
    double getEta(const int m, const int e) const;
    double getEta(const int softId) const {return towerEta[softId-1];}
    double getPhi(const int m, const int s) const;
    double getPhi(const int softId) const {return towerPhi[softId-1];}
    void setEtaPhiArrays();
    TVector3 getVector(const int softId) const;
    TVector3 getTowerPosition(const int softId, TVector3& vtx) const {return getVector(softId)-vtx;}

private:
    const float pi = TMath::Pi();

    static const int mNModule = 120; 
    static const int mNEta = 20;
    static const int mNSub = 2;
    static const int mNes = 40; //mNEta * mNSub
    static const int mNRaw = 4800; //mNModule * mNEta * mNSub

    const float mRadius = 225.405;
    const float mYWidth = 11.174;
    const float mEtaMax = 0.984;
    const float mEtaMin = 0.0035;
    const float mPhiStepHalf = pi/60.0; //(2pi/nModule)

    double mPhiOffset[2] = {(72./180.)*pi, (108./180.)*pi};
    double mPhiStep[2] = {-pi/30.0, pi/30.0};
    double mPhiBound[2] = {(75./180.)*pi, (105./180.)*pi};
    double mYlocal[2] = {-mYWidth/2.0, mYWidth/2.0};
    double mPhi[2] = {std::atan2(mYlocal[0], mRadius), std::atan2(mYlocal[1], mRadius)};
    double mYB[3] = {mYlocal[0], 0.0, mYlocal[1]};
    double mPhiB[3] = {mPhi[0], 0.0, mPhi[1]};

    const int mMaxAdc = 4096;

    double mZlocal[mNEta];
    double mEta[mNEta];
    double mEtaB[mNEta+1];
    double mPhiModule[mNModule];

    double towerEta[mNRaw];
    double towerPhi[mNRaw];

    ClassDef(BEMCLocator, 1)
};

#endif