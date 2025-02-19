#ifndef CFEMPTO
#define CFEMPTO

#include "Cutils.h"

class fempto
{
  public:
    double doInteract(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong, float sumRadii = 0, float *pos=nullptr, float *posLab=nullptr); // perform interaction and return momentum exchanged
    void doInteractAll(std::vector<particleMC>& part);

    void setParams(float strong=2.2E-3, float strongR=2.4, float coloumb=1.44E-3, float sourceRadius=0) { mStrong = strong, mStrongR = strongR, mCoulomb = coloumb, mSourceRadius = sourceRadius; }
    float getStrong() const { return mStrong; }
    float getStrongRadius() const { return mStrongR; }
    float getColoumb() const { return mCoulomb; }
    float getSourceRadius() const { return mSourceRadius; }
  private:
    float mStrong = 2.2E-3;     // attractive potential
    float mStrongR = 2.4;        // radius of box potential
    float mCoulomb = 1.44E-3;
    float mSourceRadius = 0;
};

#endif

