#ifndef CFEMPTOSOURCE
#define CFEMPTOSOURCE

#include "Cvfempto.h"
#include "CWignerUtils.h"
#include "TF1.h"
#include "TF2.h"

class femptoSource : public vfempto
{
  public:
    femptoSource() {
        wignerUtils::setParams(17.4E-3, 3.2, 1.44E-3, -3); // units in GeV (e.g., 17.4 MeV = 0.0174 GeV)
    }

    void setParams(float strong = 17.4E-3, float strongR = 3.2, float coulomb = 1.44E-3,
                   float sourceRadius = 0, float spinFact = 3. / 8) {
        wignerUtils::setParams(strong, strongR, coulomb, sourceRadius, spinFact);
    }

    bool set(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong);
    double doInteract(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong,
                      float sumRadii = 0, float* pos = nullptr, float* posLab = nullptr);

    virtual float getCoalProb(const particleMC& p1, const particleMC& p2);
    float calcProb() { return wignerUtils::getcoal(); }

    virtual void init();

    void setSourceRadius(float radius) { wignerUtils::setSourceRadius(radius); }
    void setKstar(float kstar, float kt = 1.0) {
        wignerUtils::setKstar(kstar, kt);
    }

    void setCharges(float cS, float cC) {
        wignerUtils::setParams(cS, 3.2, cC, wignerUtils::getRadius()); // 3.2 = strongR default
    }

    float kineticSource() const { return wignerUtils::kineticSource(); }
    float potentialSource() const { return wignerUtils::potentialSource(); }

    float getKstarFinal(float coalProb = 0., float massRed = 938. / 2., float boundE = 2.22) const;

    float kineticDeuteron() const { return 0; }     // Placeholder
    float potentialDeuteron() const { return 0; }   // Placeholder

};

#endif
