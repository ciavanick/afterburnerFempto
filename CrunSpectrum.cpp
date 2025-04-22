#include "CrunSpectrum.h"
#include "TFile.h"

void runSpectrum::initHistos()
{
    mHProtonSpectrum = new TH1D("mHProtonSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 100, 0., 10.);
    mHNeutronSpectrum = new TH1D("mHNeutronSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 100, 0., 10.);
    mHDeuteronSpectrum = new TH1D("mHDeuteronSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 100, 0., 10.);
    mHTritiumSpectrum = new TH1D("mHTritiumSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 100, 0., 10.);
    mHHelium3Spectrum = new TH1D("mHHelium3Spectrum" + mName, "Spectrum;p_{T} (GeV/c)", 100, 0., 10.);
    mHHelium4Spectrum = new TH1D("mHHelium4Spectrum" + mName, "Spectrum;p_{T} (GeV/c)", 100, 0., 10.);
    mHProtonSpectrumA = new TH1D("mHProtonSpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 100, 0., 10.);
    mHNeutronSpectrumA = new TH1D("mHNeutronSpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 100, 0., 10.);
    mHDeuteronSpectrumA = new TH1D("mHDeuteronSpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 100, 0., 10.);
    mHTritiumSpectrumA = new TH1D("mHTritiumSpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 100, 0., 10.);
    mHHelium3SpectrumA = new TH1D("mHHelium3SpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 100, 0., 10.);
    mHHelium4SpectrumA = new TH1D("mHHelium4SpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 100, 0., 10.);
}
//_________________________________________________________________________
void runSpectrum::initEventsHisto()
{
    mEvents = new TH1D("mEvents" + mName, "Number of events", 6, 0, 6);
    mEvents->Fill("Number of Events", 0);
    mEvents->Fill("Number of accepted Events", 0);
    mEvents->Fill("Selected Particles", 0);
    mEvents->Fill("Selected Protons", 0);
    mEvents->Fill("Selected Neutrons", 0);
    mEvents->Fill("Selected Deuterons", 0);
    mEvents->Fill("Selected Tritons", 0);
    mEvents->Fill("Selected Helium3", 0);
    mEvents->Fill("Selected Helium4", 0);
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
                mHProtonSpectrumA->Fill(p.q[ipdg].Pt());
            }else if(p.pdgOptions[ipdg] == mPDGNe){
                mEvents->Fill("Selected Neutrons", 1);
                mHNeutronSpectrum->Fill(p.q[ipdg].Pt());
                mHNeutronSpectrumA->Fill(p.q[ipdg].Pt());
            }else if(p.pdgOptions[ipdg] == mPDGDe){
                mEvents->Fill("Selected Deuterons", 1);
                mHDeuteronSpectrum->Fill(p.q[ipdg].Pt());
                mHDeuteronSpectrumA->Fill(p.q[ipdg].Pt()/2);
            }else if(p.pdgOptions[ipdg] == mPDGT){
                mEvents->Fill("Selected Tritons", 1);
                mHTritiumSpectrum->Fill(p.q[ipdg].Pt());
                mHTritiumSpectrumA->Fill(p.q[ipdg].Pt()/3);
            }else if(p.pdgOptions[ipdg] == mPDGHe3){
                mEvents->Fill("Selected Helium3", 1);
                mHHelium3Spectrum->Fill(p.q[ipdg].Pt());
                mHHelium3SpectrumA->Fill(p.q[ipdg].Pt()/3);
            }else if(p.pdgOptions[ipdg] == mPDGHe4){
                mEvents->Fill("Selected Helium4", 1);
                mHHelium4Spectrum->Fill(p.q[ipdg].Pt());
                mHHelium4SpectrumA->Fill(p.q[ipdg].Pt()/4);
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
    mHTritiumSpectrum->Write();
    mHHelium3Spectrum->Write();
    mHHelium4Spectrum->Write();
    mHProtonSpectrumA->Write();
    mHNeutronSpectrumA->Write();
    mHDeuteronSpectrumA->Write();
    mHTritiumSpectrumA->Write();
    mHHelium3SpectrumA->Write();
    mHHelium4SpectrumA->Write();
}
//_________________________________________________________________________
int runSpectrum::selectP(const particleCand &p)
{
    for (int i = 0; i < p.pdgOptions.size(); ++i)
    {
        if (p.pdgOptions[i] == mPDGPr || p.pdgOptions[i] == mPDGNe || p.pdgOptions[i] == mPDGDe || p.pdgOptions[i] == mPDGT || p.pdgOptions[i] == mPDGHe3 || p.pdgOptions[i] == mPDGHe4)
        {
            return i;
        }
    }
    return -1;
}