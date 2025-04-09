#include "CrunSpectrum.h"
#include "TFile.h"

void runSpectrum::initHistos()
{
    mHProtonSpectrum = new TH1D("mHProtonSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 62, 0.4, 3.5);
    mHNeutronSpectrum = new TH1D("mHNeutronSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 62, 0.4, 3.5);
    mHDeuteronSpectrum = new TH1D("mHDeuteronSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 31, 0.4, 3.5);
}
//_________________________________________________________________________
void runSpectrum::initEventsHisto()
{
    mEvents = new TH1D("mEvents" + mName, "Number of events", 5, 0, 5);
    mEvents->Fill("Number of Events", 0);
    mEvents->Fill("Selected Particles", 0);
    mEvents->Fill("Selected Protons", 0);
    mEvents->Fill("Selected Neutrons", 0);
    mEvents->Fill("Selected Deuterons", 0);
}
//_________________________________________________________________________
void runSpectrum::process()
{
    for (int i = 0; i < mVect.size(); ++i)
    {
        const particleCand &p = mVect[i];
        int ipdg = selectP(p);
        if (ipdg == -1)
                continue;
        if (std::abs(p.q[ipdg].Rapidity()) < 0.5)
        {
            mEvents->Fill("Selected Particles", 1);
            if (p.pdgOptions[ipdg] == mPDGPr)
            {
                mEvents->Fill("Selected Protons", 1);
                mHProtonSpectrum->Fill(p.q[ipdg].Pt());
            }else if(p.pdgOptions[ipdg] == mPDGNe){
                mEvents->Fill("Selected Neutrons", 1);
                mHNeutronSpectrum->Fill(p.q[ipdg].Pt());
            }else if(p.pdgOptions[ipdg] == mPDGDe){
                mEvents->Fill("Selected Deuterons", 1);
                mHDeuteronSpectrum->Fill(p.q[ipdg].Pt());
            }
        }
    }
}
//_________________________________________________________________________
void runSpectrum::writeHistos()
{
    mHProtonSpectrum->Write();
    mHNeutronSpectrum->Write();
    mHDeuteronSpectrum->Write();
}
//_________________________________________________________________________
int runSpectrum::selectP(const particleCand &p)
{
    for (int i = 0; i < p.pdgOptions.size(); ++i)
    {
        if (p.pdgOptions[i] == mPDGPr || p.pdgOptions[i] == mPDGNe || p.pdgOptions[i] == mPDGDe)
        {
            return i;
        }
    }
    return -1;
}