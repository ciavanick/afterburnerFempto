#include "CrunPhotons.h"
#include "TFile.h"
#include "Cutils.h"

void runPhotons::initHistos()
{
    mHPhotonsE = new TH1D("mHPhotonsE" + mName, "N;E (GeV)", 200, 0., 10.);
    mHPhotonsDeuteronsKStar = new TH1D("mHPhotonsDeuteronsKStar" + mName, "N;k* (GeV/c)", 200, 0., 10.);
    mHPhotonsDeuteronsEP = new TH2D("mHPhotonsDeuteronsEP" + mName, "N;E (GeV); p (GeV/c)", 200, 0., 10., 200, 0., 10.);
}

void runPhotons::initEventsHisto()
{
    mEvents = new TH1D("mEvents" + mName, "Number of events", 2, 0, 2);
    mEvents->Fill("Number of Photons", 0);
    mEvents->Fill("Number of Photons per Deuteron", 0);
}

void runPhotons::process()
{
    for (int i = 0; i < mVect.size(); ++i)
    {
        const particleCand &p = mVect[i];
        int ipdgPhotons = -1;
        int ipdgDeuterons = -1;
        for (int j = 0; j < p.pdgOptions.size(); ++j)
        {
            if (p.pdgOptions[i] == mPDGPhotons)
            {
                ipdgPhotons = j;
                break;
            }
        }
        for (int j = 0; j < p.pdgOptions.size(); ++j)
        {
            if (p.pdgOptions[i] == mPDGDeuterons)
            {
                ipdgDeuterons = j;
                break;
            }
        }
        if (ipdgPhotons == -1)
            continue;

        mEvents->Fill("Number of Photons", 1);

        float y = p.q[ipdgPhotons].Rapidity();
        if (y >= mMinRapidity && y <= mMaxRapidity)
        {
            mHPhotonsE->Fill(p.q[ipdgPhotons].E());
            if (ipdgDeuterons == -1) continue;
            mEvents->Fill("Number of Photons per Deuteron", 1);
            mHPhotonsDeuteronsEP->Fill(p.q[ipdgPhotons].E(), p.q[ipdgDeuterons].P());
            auto kStar = utils::getKstar(p, p, mPDGPhotons, mPDGDeuterons);
            mHPhotonsDeuteronsKStar->Fill(kStar);
        }
    }
}

void runPhotons::writeHistos()
{
    mHPhotonsE->Write();
    mHPhotonsDeuteronsKStar->Write();
    mHPhotonsDeuteronsEP->Write();
}