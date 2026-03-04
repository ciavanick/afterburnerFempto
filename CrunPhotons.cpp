#include "CrunPhotons.h"
#include "TFile.h"
#include "Cutils.h"

void runPhotons::initHistos()
{
    mHPhotonsE = new TH1D("mHPhotonsE" + mName, "N;E(#gamma) [GeV]", 2000, 0., 10.);
    mHPhotonsP = new TH1D("mHPhotonsP" + mName, "N;P(#gamma) [GeV/c]", 2000, 0., 10.);
    mHPhotonsDeuteronsKStar = new TH1D("mHPhotonsDeuteronsKStar" + mName, "N;k* (#gamma-d)  [GeV/c]", 2000, 0., 10.);
    mHPhotonsDeuteronsKStark100 = new TH1D("mHPhotonsDeuteronsKStark100" + mName, "N;k* (#gamma-d)  [GeV/c]", 2000, 0., 10.);
    mHPhotonsDeuteronsKStark50 = new TH1D("mHPhotonsDeuteronsKStark50" + mName, "N;k* (#gamma-d)  [GeV/c]", 2000, 0., 10.);
    mHProtonNeutronKStar = new TH1D("mHProtonNeutronKStar" + mName, "N;k* (p-n) [GeV/c]", 2000, 0., 10.);
    mHProtonNeutronKStark100 = new TH1D("mHProtonNeutronKStak100" + mName, "N;k* (p-n) [GeV/c]", 2000, 0., 10.);
    mHProtonNeutronKStark50 = new TH1D("mHProtonNeutronKStak50" + mName, "N;k* (p-n) [GeV/c]", 2000, 0., 10.);
    mHPhotonsDeuteronsEP = new TH2D("mHPhotonsDeuteronsEP" + mName, "N;E (#gamma) [GeV]; p (d) [GeV/c]", 1000, 0., 10., 1000, 0., 10.);
}
//_________________________________________________________________________
void runPhotons::initEventsHisto()
{
    mEvents = new TH1D("mEvents" + mName, "Number of events", 7, 0, 7);
    mEvents->Fill("Number of unweighted Events", 0);
    mEvents->Fill("Number of Events", 0);
    mEvents->Fill("Number of accepted Events", 0);
    mEvents->Fill("Selected Particles", 0);
    mEvents->Fill("Number of Photons", 0);
    mEvents->Fill("Number of Deuteron per Photons", 0);
}
//_________________________________________________________________________
void runPhotons::process()
{
    for (int i1 = 0; i1 < mVect.size(); ++i1)
    {
        const particleCand &p1 = mVect[i1];
        int ipdgPhotons = selectPhotons(p1);
        if (ipdgPhotons == -1)
            continue;

        mEvents->Fill("Number of Photons", 1);

        float y = p1.q[ipdgPhotons].Rapidity();
        if (y >= mMinRapidity && y <= mMaxRapidity)
        {
            mEvents->Fill("Selected Particles", 1);
            mHPhotonsE->Fill(p1.q[ipdgPhotons].E());
            mHPhotonsP->Fill(p1.q[ipdgPhotons].P());
            for (int i2 = 0; i2 < mVect.size(); ++i2)
            {
                const particleCand &p2 = mVect[i2];
                if (i1 == i2)
                {
                    continue;
                }

                int ipdgDeuterons = selectDeuterons(p2);
                if (ipdgDeuterons == -1)
                    continue;

                mHPhotonsDeuteronsEP->Fill(p1.q[ipdgPhotons].E(), p2.q[ipdgDeuterons].P());
                auto kStar = utils::getKstar(p1, p2, ipdgPhotons, ipdgDeuterons);
                mHPhotonsDeuteronsKStar->Fill(kStar);
                if(p1.q[ipdgPhotons].E() > 0.1) mHPhotonsDeuteronsKStark100->Fill(kStar);
                if(p1.q[ipdgPhotons].E() > 0.05) mHPhotonsDeuteronsKStark50->Fill(kStar);
                mEvents->Fill("Number of Deuteron per Photons", 1);
                TLorentzVector sum = p2.q[ipdgDeuterons] + p1.q[ipdgPhotons];
                float DeltaE = sum.M() - mProtonMass - mNeutronMass;
                float E2 = sum.M() * sum.M();
                float kstarPN = 0.5 / sum.M() * sqrt(E2 * E2 + (mNeutronMass * mNeutronMass - mProtonMass * mProtonMass) * (mNeutronMass * mNeutronMass - mProtonMass * mProtonMass) - 2 * (mNeutronMass * mNeutronMass + mProtonMass * mProtonMass) * E2);  
                if (p1.q[ipdgPhotons].E() > 0.1) mHProtonNeutronKStark100->Fill(kstarPN);                
                if (p1.q[ipdgPhotons].E() > 0.05) mHProtonNeutronKStark50->Fill(kstarPN);                
                if(kstarPN > mEnergyCut) mHProtonNeutronKStar->Fill(kstarPN); // 1. / kstarPN / kstarPN
            }
        }
    }
}
//_________________________________________________________________________
void runPhotons::writeHistos()
{
    mHPhotonsE->Write();
    mHPhotonsP->Write();
    mHPhotonsDeuteronsKStar->Write();
    mHPhotonsDeuteronsKStark100->Write();
    mHPhotonsDeuteronsKStark50->Write();
    mHPhotonsDeuteronsEP->Write();
    mHProtonNeutronKStar->Write();
    mHProtonNeutronKStark100->Write();
    mHProtonNeutronKStark50->Write();
}
//_________________________________________________________________________
int runPhotons::selectPhotons(const particleCand &p)
{
    int ipdgPhotons = -1;
    for (int j = 0; j < p.pdgOptions.size(); ++j)
    {
        if (p.pdgOptions[j] == mPDGPhotons)
        {
            ipdgPhotons = j;
            break;
        }
    }
    return ipdgPhotons;
}
//_________________________________________________________________________
int runPhotons::selectDeuterons(const particleCand &p)
{
    int ipdgDeuterons = -1;
    for (int j = 0; j < p.pdgOptions.size(); ++j)
    {
        if (p.pdgOptions[j] == mPDGDeuterons)
        {
            ipdgDeuterons = j;
            break;
        }
    }
    return ipdgDeuterons;
}
//_________________________________________________________________________