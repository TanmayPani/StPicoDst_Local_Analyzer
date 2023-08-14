#define BEMCLocator_CXX

#include "BEMCLocator.h"

using namespace std;

ClassImp(BEMCLocator)

BEMCLocator::BEMCLocator(){
    int mw;
    for(int i = 0; i < mNModule/2; i++){
        double phiW = mPhiOffset[0] + i*mPhiStep[0];
        while(phiW >= pi) phiW -= 2*pi;
        while(phiW < -pi) phiW += 2*pi;
        if(phiW > (pi-0.0001)) phiW = -pi;// -pi<=phi<phi
        mPhiModule[i] = phiW;
        double phi_w = mPhiBound[1] - phiW;
        if(phi_w < 0.0) phi_w += 2*pi;
        if(phi_w < 0.0 || phi_w > 2*pi) cout<<"phi_w = "<<phi_w<<endl;
        mw = 120 - int(phi_w/mPhiStep[1]);
        mPhiModule[mw] = phiW;
    }

    for(int i=0; i<mNEta; i++){
        mEtaB[i] = 0.05*i;
    } 
    mEtaB[mNEta]=mEtaMax; 
    mEtaB[0]=mEtaMin;

    for(int i=0; i< mNEta; i++){
        mEta[i]    = (mEtaB[i+1] + mEtaB[i])/2.;
        mZlocal[i] = mRadius * sinh(mEta[i]);  // z=r/tan(theta) => 1./tan(theta) = sinh(eta)
    }

    setEtaPhiArrays();

}

void BEMCLocator::getBin(const int softId, int &m, int &e, int &s) const {
    assert(softId >= 1 && softId <= mNRaw);
    int wid = softId - 1;
    m = wid / mNes + 1;
    int j = wid - (m-1)*mNes;
    s = j / mNEta + 1;
    e = j%mNEta + 1;
}

double BEMCLocator::getEta(const int m, const int e) const {
    assert(m >= 1 && m <= mNModule);
    assert(e >= 1 && e <= mNEta);
    return (m <= mNModule/2) ? mEta[e-1] : -mEta[e-1];
}

double BEMCLocator::getPhi(const int m, const int s) const {
    assert(m >= 1 && m <= mNModule);
    assert(s >= 1 && s <= mNSub);
    int iphi, im;
    double phiW;
    if(m <= mNModule/2){
        iphi = m-1;
        im = 0;
        phiW = -mPhi[s-1];
    }else{
        iphi = m-mNModule/2 - 1;
        im = 1;
        phiW = mPhi[s-1];
    }

    phiW += mPhiOffset[im] + iphi*mPhiStep[im];

    while(phiW >= pi) phiW -= 2*pi;
    while(phiW < -pi) phiW += 2*pi;
    if(phiW > (pi-0.0001)) phiW = -pi;// -pi<=phi<phi

    return phiW;
}

void BEMCLocator::setEtaPhiArrays(){
    for(int i = 0; i < mNRaw; i++){
        int m, e, s;
        getBin(i+1, m, e, s);
        towerEta[i] = getEta(m, e);
        towerPhi[i] = getPhi(m, s);
        //cout<<"i = "<<i<<", m = "<<m<<", e = "<<e<<", s = "<<s<<", eta = "<<towerEta[i]<<", phi = "<<towerPhi[i]<<endl;
    }
}

TVector3 BEMCLocator::getVector(const int softId) const {
    assert(softId >= 1 && softId <= mNRaw);
    double x, y, z;
    double phiW = towerPhi[softId-1];
    x = mRadius * cos(phiW);
    y = mRadius * sin(phiW);
    int m, e, s;
    getBin(softId, m, e, s);
    if(m <= mNModule/2) z = mZlocal[e-1];
    else z = -mZlocal[e-1];
    return TVector3(x, y, z); 
}
