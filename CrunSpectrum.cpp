#include "CrunSpectrum.h"
#include "TFile.h"

void runSpectrum::initSpectrum(){
    mHSpectrum = new TH1D("mHSpectrum" + mName, "Spectrum;p_{T} (GeV/c)",100, 0, 10);
    mNEvents = new TH1D("mNEvents" + mName, "Number of events",1, 0.5, 1.5);
}
//_________________________________________________________________________
void runSpectrum::process(){
    for(int i=0; i < mVect.size(); ++i){
        const particleCand& p = mVect[i];
        for(int j=0; j < p.pdgOptions.size(); ++j){
            if(std::abs(p.q[j].Eta()) < 1){
                if(p.pdgOptions[j] == mPDG){
                    mHSpectrum->Fill(p.q[j].Pt());
                }
            }
        }
    }
    mNEvents->Fill(1);
}
//_________________________________________________________________________
void runSpectrum::write(){
  gFile->mkdir("runSpectrum" + mName);
  gFile->cd("runSpectrum" + mName);
  mHSpectrum->Write();
  mNEvents->Write();
  gFile->cd();
}
//_________________________________________________________________________
void runSpectrum::setPDG(int pdg){
    mPDG = pdg;
}