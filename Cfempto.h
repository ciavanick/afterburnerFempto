#ifndef CFEMPTO
#define CFEMPTO

#include "Cvfempto.h"

class fempto : public vfempto
{
  public:
    fempto() { setParams(20E-3, 1.6, 1.44E-3, 0); }
    double doInteract(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong, float sumRadii = 0, float *pos=nullptr, float *posLab=nullptr); // perform interaction and return momentum exchanged
  private:
};

#endif

