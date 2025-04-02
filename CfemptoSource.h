#ifndef CFEMPTOSOURCE
#define CFEMPTOSOURCE

#include "Cvfempto.h"
#include "CwaveUtils.h"
#include "TF1.h"
#include "TF2.h"

class femptoSource : public vfempto
{
  public:
    femptoSource() { waveUtils::setParams(17.4, 3.2, 1.44, -1.5); } // here we work with MeV (Rsource=-1 means -> wave function defined by setting Rsource to match K.E.)
    void setParams(float strong=2.2E-3, float strongR=2.4, float coulomb=1.44E-3, float sourceRadius=0, float spinFact=3./8) { waveUtils::setParams(strong, strongR, coulomb, sourceRadius, spinFact); }
    double doInteract(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong, float sumRadii = 0, float *pos=nullptr, float *posLab=nullptr); // perform interaction and return momentum exchanged
    virtual float getCoalProb(const particleMC& p1, const particleMC& p2) { return waveUtils::getCoalProb(p1, p2); }
    float calcProb() { return waveUtils::calcProb(); }
    virtual void init();

    void setSourceRadius(float radius) { waveUtils::setSourceRadius(radius); }
    void setKstar(float kstar, float kt=1.0) { waveUtils::setKstar(kstar, kt); }
    void setCharges(float cS, float cC) { waveUtils::setCharges(cS, cC); }

    float kineticSource() const { return waveUtils::kineticSource(); }
    float potentialSource() const { return waveUtils::potentialSource(); }
    float getKstarFinal(float coalProb = 0., float massRed=waveUtils::MRED_NN, float boundE=waveUtils::EBOUND_D) const { return waveUtils::getKstarFinal(coalProb, massRed, boundE); }
    float kineticDeuteron() const { return waveUtils::kineticDeuteron(); }
    float potentialDeuteron() const { return waveUtils::potentialDeuteron(); }

  private:
};

#endif

