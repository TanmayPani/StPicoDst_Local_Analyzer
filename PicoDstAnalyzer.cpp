// This is needed for calling standalone classes (not needed on RACF)
#include "PicoDstAnalyzer.h"
#include "StPicoDstReader.h"
#include "StPicoDst.h"
#include "StPicoEvent.h"
#include "StPicoTrack.h"
#include "StPicoBTowHit.h"
#include "StPicoMcTrack.h"

#include "StRefMultCorr.h"
#include "CentralityMaker.h"
#include "BEMCLocator.h"
#include "TTreeJet.h"
#include "TTreeEvent.h"
#include "EventPlaneMaker.h"

#include "JetMaker.h"
#include "JetBackgroundMaker.h"
#include "JetVector.h"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

using namespace fastjet;

using namespace std;

map<string, function<double(JetVector&)>> PicoDstAnalyzer::jetVars = {
    {"Pt", [](JetVector& vec){return vec.pt(); }},
    {"Eta", [](JetVector& vec){return vec.eta(); }},
    {"Phi", [](JetVector& vec){return vec.phi(); }},
    {"NEF", [](JetVector& vec){return vec.getNeutralPtFraction(); }},
    {"LeSub", [](JetVector& vec){return (vec.nChargedConstituents() > 1) ? vec.getChargedConstituent(0).perp() - vec.getChargedConstituent(1).perp() : -1; }},
    {"PtD", [](JetVector& vec){return (vec.nChargedConstituents() > 1) ? vec.getChargedAngularity(2, 0, true) : -1; }},
    {"Girth", [](JetVector& vec){return (vec.nChargedConstituents() > 1) ? vec.getChargedAngularity(1, 1, false) : -1;}}
};

map<string, function<double(StPicoTrack*)>> PicoDstAnalyzer::trackVars = {
    {"Pt",  [](StPicoTrack* vec){return vec->pMom().Pt(); }},
    {"Eta", [](StPicoTrack* vec){return vec->pMom().Eta(); }},
    {"Phi", [](StPicoTrack* vec){return vec->pMom().Phi(); }},
    {"Charge", [](StPicoTrack* vec){return vec->charge(); }}
};

map<string, function<double(StPicoMcTrack*)>> PicoDstAnalyzer::genTrackVars = {
    {"Pt",  [](StPicoMcTrack* vec){return vec->fourMomentum().Pt(); }},
    {"Eta", [](StPicoMcTrack* vec){return vec->fourMomentum().Eta(); }},
    {"Phi", [](StPicoMcTrack* vec){return vec->fourMomentum().Phi(); }},
    {"E", [](StPicoMcTrack* vec){return vec->fourMomentum().E(); }},
    {"Charge", [](StPicoMcTrack* vec){return vec->charge(); }}
};

map<string, function<double(double, TVector3&)>> PicoDstAnalyzer::towerVars = {
    {"Et",  [](double Et, TVector3& vec){return Et; }},
    {"Eta", [](double Et, TVector3& vec){return vec.Eta(); }},
    {"Phi", [](double Et, TVector3& vec){return vec.Phi(); }}
};


PicoDstAnalyzer::PicoDstAnalyzer(string infile, long nEv, string outfile, double WtFactor){
    inFileName = infile;
    assert(outfile.find(".root") != string::npos);
    outFileName = outfile;
    outFileName.insert(outFileName.find(".root"), ".tree");
    histOutFileName = outfile;
    histOutFileName.insert(histOutFileName.find(".root"), ".hist");
    nEvents = nEv;
    genWeight = WtFactor;

    towerHadCorrSum.resize(4800, 0.0);
    towerNTracksMatched.resize(4800, 0);

    eventTreeArray = new TClonesArray("TTreeEvent", 1);
    jetTreeArray = new TClonesArray("TTreeJet", 20);
    genJetTreeArray = new TClonesArray("TTreeJet", 20);
}

PicoDstAnalyzer::~PicoDstAnalyzer(){
    if(outTree) delete outTree;
    if(outFile) delete outFile;
    if(picoDst) delete picoDst;
    if(picoEvent) delete picoEvent;
    if(outTree) delete outTree;
    if(outFile) delete outFile;
}

void PicoDstAnalyzer::init(){
    if(fjMaker){
        fjMaker->init();
        cout<<"Initialized detector-level JetMaker..."<<endl;
        fjMaker->printDescription();
    }
    if(fjGenMaker){
        fjGenMaker->init();
        cout<<"Initialized particle-level JetMaker..."<<endl;
        fjGenMaker->printDescription();
    }

    if(!picoReader){
        cout<<"No picoReader found. Creating a new one..."<<endl;
        picoReader.reset(new StPicoDstReader(inFileName.c_str()));
    }
    picoReader->Init();

    if(!epMaker){
        cout<<"No epMaker found. Creating a new one..."<<endl;
        epMaker.reset(new EventPlaneMaker());
    }

    declareEventPlaneHistos();

    refMultCorr.reset(CentralityMaker::instance()->getgRefMultCorr_P18ih_VpdMB30_AllLumi());
    cout<<"Set up grefmultCorr..."<<endl;
    refMultCorr->print();

    bemcLoc.reset(new BEMCLocator());

    outFile = new TFile(outFileName.c_str(), "RECREATE");
    outFile->cd();
    outTree = new TTree("JetTree", "JetTree");
    outTree->SetDirectory(gDirectory);
    outTree->Branch("Event", &eventTreeArray);
    outTree->Branch("Jets", &jetTreeArray);
    outTree->Branch("GenJets", &genJetTreeArray);

    if( !picoReader->chain() ) {cout << "No chain has been found." << endl; return;}
    unsigned long events2read = picoReader->chain()->GetEntries();
    cout << "Number of events to read: " << events2read << endl;
    if(nEvents <= 0 || nEvents > events2read) nEvents = events2read;
    cout << "Will read " << nEvents << " events." << endl;

    genWeight = genWeight/(double)events2read;
}

void PicoDstAnalyzer::clear(){
   // if(picoDst)picoDst.reset();
    if(fjMaker)fjMaker->clear();
    if(fjGenMaker)fjGenMaker->clear();
    epMaker->clear();
    eventTreeArray->Clear();
    jetTreeArray->Clear();
    genJetTreeArray->Clear();
    towerHadCorrSum.assign(towerHadCorrSum.size(), 0.0);
    towerNTracksMatched.assign(towerNTracksMatched.size(), 0);
}

void PicoDstAnalyzer::finish(){
    outFile->Write();
    outFile->Close();

    histOutFile = new TFile(histOutFileName.c_str(), "RECREATE");
    histOutFile->cd();
    for(auto& hist : hist1D){
        hist.second->Write();
    }
    for(auto& hist : hist2D){
        hist.second->Write();
    }

    pRes22->Write();
    pRes24->Write();

    histOutFile->Write();
    histOutFile->Close();

    epMaker->finish();
}

void PicoDstAnalyzer::eventLoop(){
    for(unsigned int i = 0; i < nEvents; i++){
        if(i%1000 == 0) cout << "Event " << i << endl;
        bool readEvent = picoReader->readPicoEvent(i);
        if( !readEvent ) {
            cout << "Something went wrong! Nothing to analyze..." << endl;
            break;
        }
        clear();

        picoDst = picoReader->picoDst();
        picoEvent = picoDst->event();
        if( !picoEvent ) {
            cout << "Something went wrong! PicoEvent not found..." << endl;
            break;
        }
        
        pVtx = picoEvent->primaryVertex();
        if(fabs(pVtx.z()) > absZVtxMax) continue;
        pVtx_Z = pVtx.z();
        //cout<<"Z vertex: "<<pVtx_Z<<" bin: "<<zVtxBin<<endl;  

        refMultCorr->init(picoEvent->runId());
        refMultCorr->initEvent(picoEvent->grefMult(), pVtx.z(), picoEvent->ZDCx());
        centbin16 = refMultCorr->getCentralityBin16();
        centbin9 = refMultCorr->getCentralityBin9();

        if(centbin16 < 0 || centbin9 < 0)continue;


        ref16 = 15-centbin16; 
        ref9 = 8-centbin9;

        centrality = 5.0*ref16 + 2.5;

        double refWeight = refMultCorr->getWeight();

        weight = genWeight*refWeight;

        hist1D["hCentrality"]->Fill(centrality, weight);
        hist1D["hRefMult"]->Fill(refMultCorr->getRefMultCorr(picoEvent->grefMult(), pVtx.z(), picoEvent->ZDCx(), 2), weight);

        treeEvent = static_cast<TTreeEvent*>(eventTreeArray->ConstructedAt(0));
        treeEvent->runId = picoEvent->runId();
        treeEvent->eventId = picoEvent->eventId();
        treeEvent->centrality = centrality;
        treeEvent->primaryVertexZ = pVtx_Z;
        treeEvent->genWeight = genWeight;
        treeEvent->refMultWeight = refWeight;

        trackLoop();
        towerLoop();
        if(fjMaker)jetLoop();
        genTrackLoop();
        if(fjGenMaker)genJetLoop();
        //cout<<"Going to make event plane..."<<endl;
        if((treeEvent->nDetectorJets < 1) && (treeEvent->nGenJets < 1))continue;
        makeEventPlane();

        outTree->Fill();
    }
}

void PicoDstAnalyzer::trackLoop(){
    for(unsigned int itrk = 0; itrk < picoDst->numberOfTracks(); itrk++){
        StPicoTrack* trk = picoDst->track(itrk);
        if(!trk) continue;
        if(!(trk->isPrimary())) continue;
        if(trk->nHitsFit() < nHitsFitMin) continue;
        if((trk->nHitsFit()/(double)trk->nHitsMax()) < nHitsRatioMin) continue;
        if(trk->gDCA(pVtx).Mag() > trkDCAMax) continue;

        TVector3 trkMom = trk->pMom();
        if(trkMom.Pt() < ptMin) continue;
        if(trkMom.Pt() > ptMax) continue;
        if(fabs(trkMom.Eta()) > absEtaMax) continue;

        double E = sqrt(trkMom.Mag2() + pi0mass*pi0mass);

        int towerMatched = trk->bemcTowerIndex();
        if(towerMatched >= 0){
            towerNTracksMatched[towerMatched]++;   
            towerHadCorrSum[towerMatched] += E;
        }

        epMaker->addTrack(*trk);

        fillTrackHistos(trk);

        if(!fjMaker)continue;
        fjMaker->inputForClustering(itrk, trkMom.Px(), trkMom.Py(), trkMom.Pz(), E);
    }
}

void PicoDstAnalyzer::towerLoop(){
    for(unsigned int itow = 0; itow < picoDst->numberOfBTowHits(); itow++){
        StPicoBTowHit* tow = picoDst->btowHit(itow);
        if(!tow) continue;
        if(tow->energy() < ptMin) continue;

        double E = tow->energy();
        if(towerNTracksMatched[itow] > 0){
            E -= towerHadCorrSum[itow]/towerNTracksMatched[itow];
        }
        if(E < ptMin) continue;

        TVector3 towPos = bemcLoc->getTowerPosition(itow+1, pVtx);
        double towEta = towPos.Eta();
        if(fabs(towEta) > absEtaMax) continue;

        double Et = E/cosh(towEta);
        if(Et < ptMin) continue;
        if(Et > ptMax) continue;

        double towMom = sqrt(E*E - pi0mass*pi0mass);
        towPos.SetMag(towMom);

        fillTowerHistos(Et, towPos);

        if(!fjMaker)continue;
        fjMaker->inputForClustering(-itow-2, towPos.Px(), towPos.Py(), towPos.Pz(), E);
    }
}

void PicoDstAnalyzer::genTrackLoop(){
    for(unsigned int igen = 0; igen < picoDst->numberOfMcTracks(); igen++){
        StPicoMcTrack* genTrk = picoDst->mcTrack(igen);
        if(!genTrk) continue;

        unsigned int idVtxStop = genTrk->idVtxStop();
        if(idVtxStop > 0) continue;

        TLorentzVector genTrkMom = genTrk->fourMomentum();
        if(genTrkMom.Pt() < ptMin) continue;
        if(fabs(genTrkMom.Eta()) > absEtaMax) continue;

        fillGenTrackHistos(genTrk);

        if(!fjGenMaker)continue;
        if(genTrk->charge() != 0)
            fjGenMaker->inputForClustering(igen, genTrkMom.Px(), genTrkMom.Py(), genTrkMom.Pz(), genTrkMom.E());
        else
            fjGenMaker->inputForClustering(-igen-2, genTrkMom.Px(), genTrkMom.Py(), genTrkMom.Pz(), genTrkMom.E());
    }
}

void PicoDstAnalyzer::jetLoop(){
    vector<JetVector> Jets = fjMaker->getFullJets();
    unsigned int NJets = Jets.size();
    if(NJets < 1) return;

    treeEvent->nDetectorJets = NJets;

    epMaker->setLeadingJet(Jets[0]);
    //epMaker->setSubLeadingJet(Jets[1]);

    hist1D["hNJets"]->Fill(NJets, weight);
    for(JetVector& jet : Jets){
        fillJetHistos(jet);
        TTreeJet* treeJet = static_cast<TTreeJet*>(jetTreeArray->ConstructedAt(jetTreeArray->GetEntriesFast()));
        treeJet->Pt = jet.pt();
        treeJet->Eta = jet.eta();
        treeJet->Phi = jet.phi();
        treeJet->NEF = jet.getNeutralPtFraction();
        if(jet.has_area())treeJet->Area = jet.area();
        treeJet->NNeutral = jet.nNeutralConstituents();
        treeJet->NCharged = jet.nChargedConstituents();
        treeJet->JetPtD = (jet.nChargedConstituents() > 1) ? jet.getChargedAngularity(2, 0, true) : -1;
        treeJet->JetGirth = (jet.nChargedConstituents() > 1) ? jet.getChargedAngularity(1, 1, false) : -1;
        treeJet->JetLeSub = (jet.nChargedConstituents() > 1) ? jet.getChargedConstituent(0).perp() - jet.getChargedConstituent(1).perp() : -1;
    }
}

void PicoDstAnalyzer::genJetLoop(){
    vector<JetVector> Jets = fjGenMaker->getFullJets();
    unsigned int NJets = Jets.size();
    if(NJets < 1) return;

    treeEvent->nGenJets = NJets;

    hist1D["hNGenJets"]->Fill(NJets, weight);
    for(JetVector& jet : Jets){
        fillGenJetHistos(jet);
        TTreeJet* genTreeJet = static_cast<TTreeJet*>(genJetTreeArray->ConstructedAt(genJetTreeArray->GetEntriesFast()));
        genTreeJet->Pt = jet.pt();
        genTreeJet->Eta = jet.eta();
        genTreeJet->Phi = jet.phi();
        genTreeJet->NEF = jet.getNeutralPtFraction();
        if(jet.has_area())genTreeJet->Area = jet.area();
        genTreeJet->NNeutral = jet.nNeutralConstituents();
        genTreeJet->NCharged = jet.nChargedConstituents();
        genTreeJet->JetPtD = (jet.nChargedConstituents() > 1) ? jet.getChargedAngularity(2, 0, true) : -1;
        genTreeJet->JetGirth = (jet.nChargedConstituents() > 1) ? jet.getChargedAngularity(1, 1, false) : -1;
        genTreeJet->JetLeSub = (jet.nChargedConstituents() > 1) ? jet.getChargedConstituent(0).perp() - jet.getChargedConstituent(1).perp() : -1;
    }
}

void PicoDstAnalyzer::declareEventPlaneHistos(){
    pRes22 = new TProfile("profRes22", "<cos(2(#Psi_{2, A}^{Raw} - #Psi_{2, B}^{Raw}))>", nCentBins9, centBins9);
    pRes22->Sumw2();
    pRes24 = new TProfile("profRes24", "<cos(4(#Psi_{2, A}^{Raw} - #Psi_{2, B}^{Raw}))>", nCentBins9, centBins9);
    pRes24->Sumw2();
}

void PicoDstAnalyzer::makeEventPlane(){
    //cout << "PicoDstAnalyzer::makeEventPlane" << endl;
    epMaker->calculateEventPlane(pVtx_Z, centrality);

    treeEvent->raw_Qx_2 = epMaker->getQx();
    treeEvent->raw_Qx_A_2 = epMaker->getQx_A();
    treeEvent->raw_Qx_B_2 = epMaker->getQx_B();
    treeEvent->raw_Qy_2 = epMaker->getQy();
    treeEvent->raw_Qy_A_2 = epMaker->getQy_A();
    treeEvent->raw_Qy_B_2 = epMaker->getQy_B();
    treeEvent->raw_psi_2 = epMaker->getPsi();
    treeEvent->raw_psi_A_2 = epMaker->getPsi_A();
    treeEvent->raw_psi_B_2 = epMaker->getPsi_B();
    treeEvent->eventPlaneWeight      = epMaker->getEventWeight();
    treeEvent->subEventPlaneWeight_A = epMaker->getSubEventWeight_A();
    treeEvent->subEventPlaneWeight_B = epMaker->getSubEventWeight_B();
    treeEvent->eventPlaneMult = epMaker->getEventMult();
    treeEvent->subEventPlaneMult_A = epMaker->getSubEventMult_A();
    treeEvent->subEventPlaneMult_B = epMaker->getSubEventMult_B();

    pRes22->Fill(centrality, epMaker->getEPResolution(2));
    pRes24->Fill(centrality, epMaker->getEPResolution(4));
}

StPicoDstReader* PicoDstAnalyzer::getPicoReader() { 
    if(!picoReader)picoReader.reset(new StPicoDstReader(inFileName.c_str()));
    return picoReader.get(); 
}

JetMaker* PicoDstAnalyzer::getFjWrapper() { 
    if(!fjMaker)fjMaker.reset(new JetMaker()); 
    return fjMaker.get(); 
}

JetMaker* PicoDstAnalyzer::getGenFjWrapper() {
    if(!fjGenMaker)fjGenMaker.reset(new JetMaker()); 
    return fjGenMaker.get(); 
}

EventPlaneMaker* PicoDstAnalyzer::getEPMaker() { 
    if(!epMaker)epMaker.reset(new EventPlaneMaker()); 
    return epMaker.get(); 
}

void PicoDstAnalyzer::addHist1D(string name, string title, int nBins, double xMin, double xMax){
    hist1D[name] = new TH1D(name.c_str(), title.c_str(), nBins, xMin, xMax);
    hist1D[name]->Sumw2(); 
}

void PicoDstAnalyzer::addHist2D(string name, string title, int nBinsX, double xMin, double xMax, int nBinsY, double yMin, double yMax){
    hist2D[name] = new TH2D(name.c_str(), title.c_str(), nBinsX, xMin, xMax, nBinsY, yMin, yMax);
    hist2D[name]->Sumw2();
}

void PicoDstAnalyzer::fillHist1D(string name, double x, double wt){
    if(hist1D.find(name) == hist1D.end()) return;
    hist1D[name]->Fill(x, wt);
}

void PicoDstAnalyzer::fillHist2D(string name, double x, double y, double wt){
    if(hist2D.find(name) == hist2D.end()) return;
    hist2D[name]->Fill(x, y, wt);
}

void PicoDstAnalyzer::fillTrackHistos(StPicoTrack* trk){
    if(!trk) return;
    for(auto& var : trackVars){
        fillHist1D("hTrack" + var.first, var.second(trk), weight);
    }
}

void PicoDstAnalyzer::fillTowerHistos(double towEt, TVector3& towPos){
    for(auto& var : towerVars){
        fillHist1D("hTower" + var.first, var.second(towEt, towPos), weight);
    }
}

void PicoDstAnalyzer::fillGenTrackHistos(StPicoMcTrack* trk){
    if(!trk) return;
    for(auto& var : genTrackVars){
        fillHist1D("hGenTrack" + var.first, var.second(trk), genWeight);
    }

}

void PicoDstAnalyzer::fillJetHistos(JetVector& jet){
    for(auto& var : jetVars){
        fillHist1D("hJet" + var.first, var.second(jet), weight);
        for(auto& var2 : jetVars){
            fillHist2D("h2Jet" + var.first + "v" + var2.first, var.second(jet), var2.second(jet), weight);
        }
    }
}

void PicoDstAnalyzer::fillGenJetHistos(JetVector& jet){
    for(auto& var : jetVars){
        fillHist1D("hGenJet" + var.first, var.second(jet), genWeight);
        for(auto& var2 : jetVars){
            fillHist2D("h2GenJet" + var.first + "v" + var2.first, var.second(jet), var2.second(jet), genWeight);
        }
    }
}









