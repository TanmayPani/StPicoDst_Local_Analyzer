#define EventPlaneMaker_cxx

#include "EventPlaneMaker.h"

#include "JetVector.h"

#include "TFile.h"
#include "TProfile2D.h"

#include <iostream>

using namespace std;

map<string, function<double(EventPlaneMaker&)>> EventPlaneMaker::epVars = {
    {"Qx",                    [](EventPlaneMaker& ep){return ep.getQx();              }},
    {"Qy",                    [](EventPlaneMaker& ep){return ep.getQy();              }},
    {"Qx_A",                  [](EventPlaneMaker& ep){return ep.getQx_A();            }},
    {"Qy_A",                  [](EventPlaneMaker& ep){return ep.getQy_A();            }},
    {"Qx_B",                  [](EventPlaneMaker& ep){return ep.getQx_B();            }},
    {"Qy_B",                  [](EventPlaneMaker& ep){return ep.getQy_B();            }},
    {"psi",                   [](EventPlaneMaker& ep){return ep.getPsi();             }},
    {"psi_A",                 [](EventPlaneMaker& ep){return ep.getPsi_A();           }},
    {"psi_B",                 [](EventPlaneMaker& ep){return ep.getPsi_B();           }},
    {"eventPlaneWeight",      [](EventPlaneMaker& ep){return ep.getEventWeight();     }},
    {"subEventPlaneWeight_A", [](EventPlaneMaker& ep){return ep.getSubEventWeight_A();}},
    {"subEventPlaneWeight_B", [](EventPlaneMaker& ep){return ep.getSubEventWeight_B();}},
    {"eventPlaneMult",        [](EventPlaneMaker& ep){return ep.getEventMult();       }},
    {"subEventPlaneMult_A",   [](EventPlaneMaker& ep){return ep.getSubEventMult_A();  }},
    {"subEventPlaneMult_B",   [](EventPlaneMaker& ep){return ep.getSubEventMult_B();  }}
};

EventPlaneMaker::EventPlaneMaker(unsigned int n){
    N = n;
    cout<<"EventPlaneMaker::EventPlaneMaker() N = "<<N<<endl;
}

EventPlaneMaker::~EventPlaneMaker(){

}

void EventPlaneMaker::clear(){
    trackVector.clear();
    leadingJet.reset();
    subLeadingJet.reset();
}

void EventPlaneMaker::finish(){
    if(outFileName == "") outFileName = "EventPlaneMaker.root";
    outFile = new TFile(outFileName.c_str(), "RECREATE");
    outFile->cd();
    for(auto& p : epProf){
        p.second->Write();
    }
    outFile->Write();
    outFile->Close();
}

void EventPlaneMaker::declareTProfile2Ds(string var1name, int nVar1Bins, const double* var1Bins, string var2name, int nVar2Bins, const double* var2Bins){
    string hname, htitle;
    for(auto& epVar : epVars){
        hname = "p2"+epVar.first+"_raw";
        htitle = "<" + epVar.first + "> raw" + " vs " + var1name + " vs " + var2name;
        epProf[hname]   = new TProfile2D(hname.c_str(), htitle.c_str(), nVar1Bins, var1Bins, nVar2Bins, var2Bins);
        epProf[hname]->Sumw2();
    }
}

void EventPlaneMaker::setLeadingJet(JetVector& jet){
    leadingJet.reset(new JetVector(jet));
}

void EventPlaneMaker::setSubLeadingJet(JetVector& jet){
    subLeadingJet.reset(new JetVector(jet));
}

void EventPlaneMaker::calculateEventPlane(double var1, double var2, double weight){
    Qx_raw   = 0 ; Qy_raw   = 0 ;
    Qx_raw_A = 0 ; Qy_raw_A = 0 ;
    Qx_raw_B = 0 ; Qy_raw_B = 0 ;

    eventWeight = 0.0;
    subEventWeight_A = 0.0;
    subEventWeight_B = 0.0;

    eventMultiplicity = 0;
    subEventMultiplicity_A = 0;
    subEventMultiplicity_B = 0;

    //cout<<"EventPlaneMaker::calculateEventPlane()"<<N<<endl;

    for(auto& trk : trackVector){
        if(!trk.isPrimary()) continue;
        TVector3 mom = trk.pMom();
        double trkPt = mom.Perp();
        //cout<<"trkPt = "<<trkPt<<endl;
        if(trkPt > maxTrackPt) continue;
        double trkEta = mom.Eta();
        double trkPhi = mom.Phi();
        if(trkPhi < 0) trkPhi += 2*TMath::Pi();

        if(leadingJet){
            if(removeLeadingEtaStrip){
                if(fabs(trkEta - leadingJet->eta()) < leadingJet->getRadius()) continue;
            }
            if(removeLeadingEtaPhiCone){
                if(leadingJet->getDeltaR(trkEta, trkPhi) < leadingJet->getRadius()) continue;
            }
        }

        if(subLeadingJet){
            if(removeSubLeadingEtaStrip){
                if(fabs(trkEta - subLeadingJet->eta()) < subLeadingJet->getRadius()) continue;
            }
            if(removeSubLeadingEtaPhiCone){
                if(subLeadingJet->getDeltaR(trkEta, trkPhi) < subLeadingJet->getRadius()) continue;
            }
        }

        double x = trkPt * cos(N * trkPhi);
        double y = trkPt * sin(N * trkPhi);

        epProf["p2Qx_raw"]->Fill(var1, var2, x, weight);
        epProf["p2Qy_raw"]->Fill(var1, var2, y, weight);

        Qx_raw += x;
        Qy_raw += y;
        eventWeight += trkPt;
        eventMultiplicity++;

        if(trkEta > 0){
            epProf["p2Qx_A_raw"]->Fill(var1, var2, x, weight);
            epProf["p2Qy_A_raw"]->Fill(var1, var2, y, weight);
            Qx_raw_A += x;
            Qy_raw_A += y;
            subEventWeight_A += trkPt;
            subEventMultiplicity_A++;
        }
        else{
            epProf["p2Qx_B_raw"]->Fill(var1, var2, x, weight);
            epProf["p2Qy_B_raw"]->Fill(var1, var2, y, weight);
            Qx_raw_B += x;
            Qy_raw_B += y;
            subEventWeight_B += trkPt;
            subEventMultiplicity_B++;
        }
    }

    psi_raw = atan2(Qy_raw, Qx_raw) / (float)N;
    psi_raw_A = atan2(Qy_raw_A, Qx_raw_A) / (float)N;
    psi_raw_B = atan2(Qy_raw_B, Qx_raw_B) / (float)N;

    if(psi_raw < -0.5*TMath::Pi()) psi_raw += TMath::Pi();
    if(psi_raw >  0.5*TMath::Pi()) psi_raw -= TMath::Pi();

    if(psi_raw_A < -0.5*TMath::Pi()) psi_raw_A += TMath::Pi();
    if(psi_raw_A >  0.5*TMath::Pi()) psi_raw_A -= TMath::Pi();

    if(psi_raw_B < -0.5*TMath::Pi()) psi_raw_B += TMath::Pi();
    if(psi_raw_B >  0.5*TMath::Pi()) psi_raw_B -= TMath::Pi();
    //cout<<"EventPlaneMaker::calculateEventPlane() psi_raw = "<<psi_raw<<endl;

    for(auto& epVar : epVars){
        string hname = "p2"+epVar.first+"_raw";
        if(hname.find("p2Q") != string::npos)continue;
        assert(epProf.find(hname) != epProf.end());
        epProf[hname]->Fill(var1, var2, epVar.second(*this), weight);
    }
}
