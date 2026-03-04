#ifndef CRUNPHOTONS
#define CRUNPHOTONS

#include "Cvrun.h"
#include "TH1D.h"

class runPhotons : public vrun
{
public:
    runPhotons(int nmix = 10, TString dirName = "runPhotons") : vrun(nmix, dirName) {}
    void setRapidityRange(float minRapidity, float maxRapidity) { mMinRapidity = minRapidity, mMaxRapidity = maxRapidity; }
    void SetEnergyCut(float energyCut = 0.){ mEnergyCut = energyCut; }
private:
    TH1D *mHPhotonsE = nullptr;
    TH1D *mHPhotonsP = nullptr;
    TH1D *mHPhotonsDeuteronsKStar = nullptr;
    TH1D *mHPhotonsDeuteronsKStark100 = nullptr;
    TH1D *mHPhotonsDeuteronsKStark50 = nullptr;
    TH1D *mHProtonNeutronKStar = nullptr;
    TH1D *mHProtonNeutronKStark100 = nullptr;
    TH1D *mHProtonNeutronKStark50 = nullptr;

    TH2D *mHPhotonsDeuteronsEP = nullptr;

    void initHistos() override;
    void initEventsHisto() override;
    void process() override;
    void writeHistos() override;

    int selectPhotons(const particleCand &p);
    int selectDeuterons(const particleCand &p);

    int mPDGPhotons = 22;
    int mPDGDeuterons = 4324;

    float mMinRapidity = -0.5;
    float mMaxRapidity = 0.5;
    float mEnergyCut = 0.;

    float mProtonMass = utils::getMass(2212);
    float mNeutronMass = utils::getMass(2112);


};

#endif