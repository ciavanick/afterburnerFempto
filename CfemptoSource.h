#ifndef CFEMPTOSOURCE
#define CFEMPTOSOURCE

#include "Cvfempto.h"
#include "CwignerUtils.h"
#include "CwaveUtils.h"
#include "TF1.h"
#include "TF2.h"

using method = waveUtils;

class femptoSource : public vfempto
{
  public:
    femptoSource() {
      if constexpr (std::is_same_v<method, waveUtils>) {
        method::setParams(17.4, 3.2, 1.44, -3); // here we work with MeV (Rsource=-1 means -> wave function defined by setting Rsource to match K.E.)
      }
      if constexpr (std::is_same_v<method, wignerUtils>) {
        method::setParams(17.4E-3, 3.2, 1.44E-3, -3); // units in GeV (e.g., 17.4 MeV = 0.0174 GeV)
      }
    }

    void setParams(float strong=2.2E-3, float strongR=2.4, float coulomb=1.44E-3, float sourceRadius=0, float spinFact=3./8) { method::setParams(strong, strongR, coulomb, sourceRadius, spinFact);}
    double doInteract(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong,
                      float sumRadii = 0, float* pos = nullptr, float* posLab = nullptr);
    virtual float getCoalProb(const particleMC& p1, const particleMC& p2) { return method::getCoalProb(p1, p2); }
    float calcProb() { return method::calcProb(); }

    bool set(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong);

    virtual void init();

    void setSourceRadius(float radius) { method::setSourceRadius(radius); }
    void setKstar(float kstar, float kt=1.0, utils::type system=utils::nn) { method::setKstar(kstar, kt, system); }

    void setCharges(float cS, float cC) {
      if constexpr (std::is_same_v<method, waveUtils>) {
        method::setCharges(cS, cC);
      }
      if constexpr (std::is_same_v<method, wignerUtils>) {
        method::setParams(cS, 3.2, cC, wignerUtils::getRadius()); // 3.2 = strongR default
      }
    }

    float kineticSource() const { return method::kineticSource(); }
    float potentialSource() const { return method::potentialSource(); }

    float getKstarFinal(float coalProb = 0., float massRed=waveUtils::MRED_NN, float boundE=waveUtils::EBOUND_D) const { return method::getKstarFinal(coalProb, massRed, boundE); }

    float kineticDeuteron() const { return method::kineticDeuteron(); }
    float potentialDeuteron() const { return method::potentialDeuteron(); }

};

#endif
