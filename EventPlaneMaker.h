#ifndef EventPlaneMaker_H
#define EventPlaneMaker_H

#include <map>
#include <vector>
#include <memory>
#include <string>
#include <functional>

#include "StPicoTrack.h"

class JetVector;
class TProfile2D;
class TFile;

class EventPlaneMaker {
public:
    EventPlaneMaker(unsigned int n = 2);
    virtual ~EventPlaneMaker();

    void clear();
    void finish();
    void declareTProfile2Ds(std::string var1name, int nVar1Bins, const double* var1Bins, std::string var2name, int nVar2Bins, const double* var2Bins);
    void addTrack(StPicoTrack& trk){trackVector.emplace_back(trk);}
    void calculateEventPlane(double v1, double v2, double weight = 1.0);

    double getQx(){return Qx_raw;}
    double getQy(){return Qy_raw;}
    double getQx_A(){return Qx_raw_A;}
    double getQy_A(){return Qy_raw_A;}
    double getQx_B(){return Qx_raw_B;}
    double getQy_B(){return Qy_raw_B;}

    double getEventWeight(){return eventWeight;}
    double getSubEventWeight_A(){return subEventWeight_A;}
    double getSubEventWeight_B(){return subEventWeight_B;}

    unsigned int getEventMult(){return eventMultiplicity;}
    unsigned int getSubEventMult_A(){return subEventMultiplicity_A;}
    unsigned int getSubEventMult_B(){return subEventMultiplicity_B;}

    double getPsi(){return psi_raw;}
    double getPsi_A(){return psi_raw_A;}
    double getPsi_B(){return psi_raw_B;}

    double getEPResolution(unsigned int m){return cos(m * (psi_raw_A - psi_raw_B));}

    void setN(unsigned int n){N = n;}
    void setMaxTrackPt(double pt){maxTrackPt = pt;}
    void setOutFileName(std::string name){outFileName = name;}

    void setLeadingJet(JetVector& jet);
    void setSubLeadingJet(JetVector& jet);

    void setRemoveLeadingEtaStrip(bool remove){removeLeadingEtaStrip = remove; removeLeadingEtaPhiCone = !remove;}
    void setRemoveSubLeadingEtaStrip(bool remove){removeSubLeadingEtaStrip = remove; removeSubLeadingEtaPhiCone = !remove;}
    void setRemoveLeadingEtaPhiCone(bool remove){removeLeadingEtaPhiCone = remove; removeLeadingEtaStrip = !remove;}
    void setRemoveSubLeadingEtaPhiCone(bool remove){removeSubLeadingEtaPhiCone = remove; removeSubLeadingEtaStrip = !remove;}

private:
    bool removeLeadingEtaStrip = true;
    bool removeSubLeadingEtaStrip = false;
    bool removeLeadingEtaPhiCone = false;
    bool removeSubLeadingEtaPhiCone = false;

    double maxTrackPt = 1.0;

    std::string outFileName = "";
    TFile* outFile = nullptr;

    std::vector<StPicoTrack> trackVector;

    std::unique_ptr<JetVector> leadingJet;
    std::unique_ptr<JetVector> subLeadingJet;

    unsigned int N = 2;

    double Qx_raw = 0, Qy_raw = 0;
    double Qx_raw_A = 0, Qy_raw_A = 0;
    double Qx_raw_B = 0, Qy_raw_B = 0;

    double psi_raw = 0;
    double psi_raw_A = 0;
    double psi_raw_B = 0;

    double eventWeight = 0.0;
    double subEventWeight_A = 0.0;
    double subEventWeight_B = 0.0;

    unsigned int eventMultiplicity = 0;
    unsigned int subEventMultiplicity_A = 0;
    unsigned int subEventMultiplicity_B = 0;

    std::map<std::string, TProfile2D*> epProf;

    static std::map<std::string, std::function<double(EventPlaneMaker&)>> epVars;
    
};

#endif