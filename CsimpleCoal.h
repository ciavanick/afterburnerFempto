#ifndef CSIMPLECOAL
#define CSIMPLECOAL

#include "Cvcoal.h"

class simpleCoal : public vcoal
{
  public:
    bool docoal(const particleMC& p1, const particleMC& p2); // perform coalescence

    void setParams(float *params);
    float getParams(int ipar) const;
  private:
};

#endif

