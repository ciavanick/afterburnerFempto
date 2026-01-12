#include "CrunPhotons.h"
#include "TFile.h"
#include "Cutils.h"

void runPhotons::initHistos()
{
    mHPhotonsE = new TH1D("mHPhotonsE" + mName, "N;E (GeV)", 200, 0., 10.);
    mHPhotonsP = new TH1D("mHPhotonsE" + mName, "N;P (GeV/c)", 200, 0., 10.);
    mHPhotonsDeuteronsKStar = new TH1D("mHPhotonsDeuteronsKStar" + mName, "N;k* (GeV/c)", 200, 0., 10.);
    mHPhotonsDeuteronsEP = new TH2D("mHPhotonsDeuteronsEP" + mName, "N;E (GeV); p (GeV/c)", 200, 0., 10., 200, 0., 10.);
}

void runPhotons::initEventsHisto()
{
    mEvents = new TH1D("mEvents" + mName, "Number of events", 5, 0, 5);
    mEvents->Fill("Number of unweighted Events", 0);
    mEvents->Fill("Number of Events", 0);
    mEvents->Fill("Number of accepted Events", 0);
    mEvents->Fill("Number of Photons", 0);
    mEvents->Fill("Number of Photons per Deuteron", 0);
}

void runPhotons::process()
{
    for (int i1 = 0; i1 < mVect.size(); ++i1)
    {
        const particleCand &p1 = mVect[i1];
        int ipdgPhotons = -1;
        for (int j = 0; j < p1.pdgOptions.size(); ++j)
        {
            if (p1.pdgOptions[i1] == mPDGPhotons)
            {
                ipdgPhotons = j;
                break;
            }
        }
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

                int ipdgDeuterons = -1;
                for (int j = 0; j < p2.pdgOptions.size(); ++j)
                {
                    if (p2.pdgOptions[i2] == mPDGDeuterons)
                    {
                        ipdgDeuterons = j;
                        break;
                    }
                }
                if (ipdgDeuterons == -1)
                    continue;

                mHPhotonsDeuteronsEP->Fill(p1.q[ipdgPhotons].E(), p2.q[ipdgDeuterons].P());
                auto kStar = utils::getKstar(p1, p2, mPDGPhotons, mPDGDeuterons);
                mHPhotonsDeuteronsKStar->Fill(kStar);
            }
        }
    }
}

void runPhotons::writeHistos()
{
    mHPhotonsE->Write();
    mHPhotonsP->Write();
    mHPhotonsDeuteronsKStar->Write();
    mHPhotonsDeuteronsEP->Write();
}