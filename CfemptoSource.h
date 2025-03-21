#ifndef CFEMPTOSOURCE
#define CFEMPTOSOURCE

#include "Cvfempto.h"
#include "TF1.h"
#include "TF2.h"

class femptoSource : public vfempto
{
  public:
    femptoSource() { setParams(33.7, 2.1, 1.44, -2); } // here we work with MeV (Rsource=-1 means -> wave function defined by setting Rsource to match K.E.)
    double doInteract(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong, float sumRadii = 0, float *pos=nullptr, float *posLab=nullptr); // perform interaction and return momentum exchanged
    virtual float getCoalProb(const particleMC& p1, const particleMC& p2);
    float calcProb();
    virtual void init();

    void setSourceRadius(float radius);
    void setKstar(float kstar);
    void setCharges(float cS, float cC);

    TF1 *getCoalRe() { return mCoalescenceRe; }
    TF1 *getCoalIm() { return mCoalescenceIm; }

    TF1 *getUDeuteron() { return mUDeuteron; }
    TF1 *getDeuteron() { return mDeuteron; }
    TF1 *getDeuteron2int() { return mDeuteron2int; }
    TF1 *getDeuteronKin() { return mDeuteronKin; }
    TF1 *getDeuteronV() { return mDeuteronV; }

    TF1 *getSource () { return mSource; }
    TF1 *getUSource () { return mUSource; }
    TF1 *getSource2int () { return mSource; }
    TF1 *getSourceKin () { return mSourceKin; }
    TF1 *getSourceV () { return mSourceV; }
    float kineticSource() const { double err; return mSourceKin->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err); }
    float potentialSource() const { double err; return mSourceV->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err); }
    float getKstarFinal(float coalProb = 0.) const;
    float kineticDeuteron() const { double err; return mDeuteronKin->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err); }
    float potentialDeuteron() const { double err; return mDeuteronV->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err); }

  private:
    static constexpr double HCUT = 197.3; // [fm*Mev]
    static constexpr double EBOUND = -2.22; // [MeV] deuteron binding energy
    static constexpr float MRED = 938/2.; // [MeV] reducted mass
    // detueron wave-function params
    double mK1 = 0; // sqrt(2*MRED*(EBOUND+V0))/HCUT;
    double mK2 = 0; // sqrt(-2*MRED*EBOUND)/HCUT;
    double mNormRight = 0; // sin(k1*mStrongR) * TMath::Exp(k2*mStrongR);

    // coalescence
    static double coalescenceRe(double *x,double *pm);
    static double coalescenceIm(double *x,double *pm);
    TF1 *mCoalescenceRe;
    TF1 *mCoalescenceIm;

    // deuteron functions
    static double uDeuteron(double *x,double *pm);
    static double waveDeuteron(double *x,double *pm);
    static double intDeuteron(double *x,double *pm);
    static double deuteronKin(double *x,double *pm);
    static double deuteronV(double *x,double *pm);
    TF1 *mDeuteron;
    TF1 *mUDeuteron;
    TF1 *mDeuteron2int;
    TF1 *mDeuteronKin;
    TF1 *mDeuteronV;

    // source functions
    float mMaxIntRange = 20;
    static double source(double *x,double *pm);
    static double usource(double *x,double *pm);
    static double source2int(double *x,double *pm);
    static double sourceKin(double *x,double *pm);
    static double sourceV(double *x,double *pm);
    static double sourceCos(double *x,double *pm);
    static double sourceSin(double *x,double *pm);
    TF1 *mSource;
    TF1 *mUSource;
    TF1 *mSource2int;
    TF1 *mSourceKin;
    TF1 *mSourceV;
};

#endif

