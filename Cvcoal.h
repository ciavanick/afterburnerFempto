#ifndef CVCOAL
#define CVCOAL

#include "Cutils.h"

class vcoal // this is a virtual class
{
  public:
    virtual bool docoal(const particleMC& p1, const particleMC& p2) = 0; // perform coalescence
    virtual void doCoalAll(std::vector<particleMC>& part);
    virtual particleMC merge(const particleMC& p1, const particleMC& p2);

    virtual void setParams(float *params) = 0;
    virtual float getParams(int ipar) const = 0;

    void setMassMax(float val) { mMassMax = val; }
    float getMassMax() const { return mMassMax; }
  private:
    float mMassMax = 999;  // maximum mass of nuclei which can be formed via coalescence
};

#endif

