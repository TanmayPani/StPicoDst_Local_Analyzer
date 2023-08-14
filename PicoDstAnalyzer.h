#ifndef PicoDstAnalyzer_H
#define PicoDstAnalyzer_H

#include "TVector3.h"

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>

class TClonesArray;
class TTree;
class TFile;

class StPicoDst;
class StPicoDstReader;
class StPicoEvent;
class StPicoTrack;
class StPicoMcTrack;

class StRefMultCorr;
class BEMCLocator;

class JetMaker;
class JetBackgroundMaker;
class JetVector;

class TTreeEvent;

class EventPlaneMaker;

class TH1D;
class TH2D;
class TProfile;
class TProfile2D;

class PicoDstAnalyzer {
public:
    PicoDstAnalyzer(std::string infileName, long nEv = -1, std::string outfileName = "test.root", double WtFactor = 1.0);
    virtual ~PicoDstAnalyzer();

    void run(){init(); eventLoop(); finish();}
    void init();
    void finish();
    void eventLoop();

    StPicoDstReader* getPicoReader();
    JetMaker* getFjWrapper(); 
    JetMaker* getGenFjWrapper();

    TTreeEvent* treeEvent = nullptr;

    EventPlaneMaker* getEPMaker();

    void setAbsZVtxMax(double zVtxMax){absZVtxMax = zVtxMax;}
    void setPtMin(double pt){ptMin = pt;}
    void setPtMax(double pt){ptMax = pt;}
    void setAbsEtaMax(double eta){absEtaMax = eta;}
    void setNHitsFitMin(int nHitsFit){nHitsFitMin = nHitsFit;}
    void setNHitsRatioMin(double nHitsRatio){nHitsRatioMin = nHitsRatio;}
    void setTrackDCAMax(double dca){trkDCAMax = dca;}

    void addHist1D(std::string name, std::string title, int nBins, double xMin, double xMax);
    void addHist2D(std::string name, std::string title, int nBinsX, double xMin, double xMax, int nBinsY, double yMin, double yMax);

private:
    void clear();
    void makeTree();
    void trackLoop();
    void towerLoop();
    void genTrackLoop();
    void jetLoop();
    void genJetLoop();
    
    void declareEventPlaneHistos();
    void makeEventPlane();

    void fillHist1D(std::string name, double x, double w = 1.0);
    void fillHist2D(std::string name, double x, double y, double w = 1.0);
    void fillTrackHistos(StPicoTrack* trk);
    void fillTowerHistos(double towEt, TVector3& towPos);
    void fillGenTrackHistos(StPicoMcTrack* trk);
    void fillJetHistos(JetVector& jet);
    void fillGenJetHistos(JetVector& jet);

    double pi0mass = 0.13957;

    std::unique_ptr<StPicoDstReader> picoReader;
    StPicoDst* picoDst;
    StPicoEvent* picoEvent;

    std::unique_ptr<StRefMultCorr> refMultCorr;
    std::unique_ptr<BEMCLocator> bemcLoc;

    std::unique_ptr<JetMaker> fjMaker;
    std::unique_ptr<JetMaker> fjGenMaker;

    std::unique_ptr<EventPlaneMaker> epMaker;

    TClonesArray* eventTreeArray = nullptr;
    TClonesArray* jetTreeArray = nullptr;
    TClonesArray* genJetTreeArray = nullptr;

    TTree* outTree = nullptr;
    TFile* outFile = nullptr;
    TFile* histOutFile = nullptr;
    TFile* eventPlaneOutFile = nullptr;

    std::string inFileName = "";
    std::string outFileName = "";
    std::string histOutFileName = "";
    std::string eventPlaneOutFileName = "";

    std::vector<double> towerHadCorrSum;
    std::vector<unsigned int> towerNTracksMatched;

    long nEvents = 10;
    TVector3 pVtx;
    double pVtx_Z = -999;
    double absZVtxMax = 30.0;

    int centbin9 = -1;
    int centbin16 = -1;
    int ref9 = -1;
    int ref16 = -1;
    float centrality = -1;

    const int nCentBins9 = 9;
    double centBins9[10] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

    double genWeight = 1.0;
    double weight = 1.0;

    double ptMin = 0.2;
    double ptMax = 30.0;
    double absEtaMax = 1.0;
    int nHitsFitMin = 15;
    double nHitsRatioMin = 0.52;
    double trkDCAMax = 3.0;

    static std::map<std::string, std::function<double(JetVector&)>> jetVars;
    static std::map<std::string, std::function<double(StPicoTrack*)>> trackVars;
    static std::map<std::string, std::function<double(StPicoMcTrack*)>> genTrackVars;
    static std::map<std::string, std::function<double(double, TVector3&)>> towerVars;

    std::map<std::string, TH1D*> hist1D;
    std::map<std::string, TH2D*> hist2D;

    TProfile *pRes22;
    TProfile *pRes24;
};

#endif