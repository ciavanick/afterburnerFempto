#ifndef CVFEMPTO
#define CVFEMPTO

#include "Cutils.h"
#include "TH1F.h"

class vfempto
{
  public:
    virtual double doInteract(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong, float sumRadii = 0, float *pos=nullptr, float *posLab=nullptr); // perform interaction and return momentum exchanged
    bool isCoalescence(const particleMC& p1, const particleMC& p2) const;
    virtual float getCoalProb(const particleMC& p1, const particleMC& p2) const;
    void doInteractAll(std::vector<particleMC>& part);

    void setParams(float strong=2.2E-3, float strongR=2.4, float coloumb=1.44E-3, float sourceRadius=0) { mStrong = strong, mStrongR = strongR, mCoulomb = coloumb, mSourceRadius = sourceRadius; }
    float getStrong() const { return mStrong; }
    float getStrongRadius() const { return mStrongR; }
    float getColoumb() const { return mCoulomb; }
    float getSourceRadius() const { return mSourceRadius; }
    TH1F* getHistoGroup() { return mHsizeGroup; }
    TH1F* getHistoMerge() { return mHsizeMerge; }
  protected:
    float mStrong = 2.2E-3;     // attractive potential
    float mStrongR = 2.4;        // radius of box potential
    float mCoulomb = 1.44E-3;
    float mSourceRadius = 0;

    particleMC merge(const particleMC& p1, const particleMC& p2);

  private:
    TH1F *mHsizeGroup = nullptr;
    TH1F *mHsizeMerge = nullptr;
};

#endif

