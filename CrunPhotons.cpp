#include "CrunPhotons.h"
#include "TFile.h"

void runPhotons::initHistos()
{
    mHPhotonsE = new TH1D("mHPhotonsE" + mName, "N;E (GeV)", 200, 0., 10.);
}

void runPhotons::initEventsHisto()
{
    mEvents = new TH1D("mEvents" + mName, "Number of events", 1, 0, 1);
    mEvents->Fill("Number of Photons", 0);
}

void runPhotons::process()
{
    for (int i = 0; i < mVect.size(); ++i)
    {
        const particleCand &p = mVect[i];
        int ipdg = -1;
        for (int j = 0; j < p.pdgOptions.size(); ++j)
        {
            if (p.pdgOptions[i] == mPDGPhotons)
            {
                ipdg = j;
                break;
            }
        }
        if (ipdg == -1)
            continue;

        mEvents->Fill("Number of Photons", 1);

        float y = p.q[ipdg].Rapidity();
        if (y >= mMinRapidity && y <= mMaxRapidity)
        {
        }
    }
}

void runPhotons::writeHistos()
{
}