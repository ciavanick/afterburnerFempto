#ifndef CRUNPHOTONS
#define CRUNPHOTONS

#include "Cvrun.h"
#include "TH1D.h"

class runPhotons : public vrun
{
public:
    runPhotons(int nmix = 10, TString dirName = "runPhotons") : vrun(nmix, dirName) {}
    void setRapidityRange(float minRapidity, float maxRapidity) { mMinRapidity = minRapidity, mMaxRapidity = maxRapidity; }

private:
    TH1D *mHPhotonsE = nullptr;
    TH1D *mHPhotonsDeuteronsKStar = nullptr;
    
    
    TH2D *mHPhotonsDeuteronsEP = nullptr;

    void initHistos() override;
    void initEventsHisto() override;
    void process() override;
    void writeHistos() override;

    int mPDGPhotons = 22;
    int mPDGDeuterons = 4324;

    float mMinRapidity = -0.5;
    float mMaxRapidity = 0.5;
};

#endif