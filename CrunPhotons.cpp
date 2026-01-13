#include "CrunPhotons.h"
#include "TFile.h"
#include "Cutils.h"

void runPhotons::initHistos()
{
    mHPhotonsE = new TH1D("mHPhotonsE" + mName, "N;E (GeV)", 1000, 0., 10.);
    mHPhotonsP = new TH1D("mHPhotonsP" + mName, "N;P (GeV/c)", 1000, 0., 10.);
    mHPhotonsDeuteronsKStar = new TH1D("mHPhotonsDeuteronsKStar" + mName, "N;k* (GeV/c)", 1000, 0., 10.);
    mHPhotonsDeuteronsEP = new TH2D("mHPhotonsDeuteronsEP" + mName, "N;E (GeV); p (GeV/c)", 1000, 0., 10., 1000, 0., 10.);
}

void runPhotons::initEventsHisto()
{
    mEvents = new TH1D("mEvents" + mName, "Number of events", 6, 0, 6);
    mEvents->Fill("Number of unweighted Events", 0);
    mEvents->Fill("Number of Events", 0);
    mEvents->Fill("Number of accepted Events", 0);
    mEvents->Fill("Selected Particles", 0);
    mEvents->Fill("Number of Photons", 0);
    mEvents->Fill("Number of Deuteron per Photons", 0);
    mEvents->Fill("Number of Deuteron", 0);
}

void runPhotons::process()
{
    //----------------------------------------- this part has to vanish
    for (int i = 0; i < mVect.size(); ++i)
    {
        const particleCand &p = mVect[i];
        int ipdgDeuterons = selectDeuterons(p);
        if (ipdgDeuterons == -1)
            continue;
        float y = p.q[ipdgDeuterons].Rapidity();
        if (y >= mMinRapidity && y <= mMaxRapidity){
            mEvents->Fill("Number of Deuteron", 1);
        }  
    }
    //------------------------------------------------- 
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
                mEvents->Fill("Number of Deuteron per Photons", 1);
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