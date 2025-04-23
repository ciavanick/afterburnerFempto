#ifndef CWIGNERUTILS
#define CWIGNERUTILS

#include "TF2.h"
#include "TFile.h"
#include "TH2.h"

class wignerUtils{
    public:
        static void init();
        static void setFunctionsParameters();
        static void setRanges(double xmin, double ymin, double xmax, double ymax);
        static void setMu(double mu);
        static void setRWidth(double rWidth);
        static void setV0(double v0);
        static double wignerSource(double *x, double *pm);
        static double wignerSource2(double *x, double *pm);
        static double jacobianFun(double *x, double *pm);
        static double jacobianW2(double *x, double *pm);
        static double kineticEnergy(double *x, double *pm);
        static double potentialEnergy(double *x, double *pm);
        static double hamiltonian(double *x, double *pm);
        static double wK(double *x, double *pm);
        static double wH(double *x, double *pm);
        static double wV(double *x, double *pm);
        static double wignerDeuteron(double *x, double *pm);
        static double wignerDeuteronIntegral(double *x, double *pm);
        static double coalescenceProbability(double *x, double *pm);
        static double radius(double k, double r0);
        static double kStarEff(double k, double radius);
        static void setRadiusK(double k);

        static double integral(TF2 *function, double minX = mRMin, double maxX = mRMax, double minP = mPMin, double maxP = mPMax);

        static double getMinX();
        static double getMaxX();
        static double getMinP();
        static double getMaxP();
        static double getHCut();

        static void setMinX(double minX);
        static void setMaxX(double maxX);
        static void setMinP(double minP);
        static void setMaxP(double maxP);
        static void setIntegrationRanges(double minX, double maxX, double minP, double maxP);

        static void setParams(float strong=17.4, float strongR=3.2, float coloumb=1.44, float sourceRadius=0, float spinFact=3./8);
        static void setSourceRadius(float radius);
        static void setKstar(float kstar, float kt);
        static void setRadiusK(float kstar);

        static float kineticSource();
        static float potentialSource();
        static float kineticDeuteron();
        static float potentialDeuteron();

        static TF2* getWignerFunction();
        static TF2* getWignerFunctionForItself();
        static TF2* getWignerFunctionForJacobian();
        static TF2* getWignerFunction2ForJacobian();
        static TF2* getKineticEnergyFunction();
        static TF2* getPotentialEnergyFunction();
        static TF2* getHamiltonianFunction();
        static TF2* getWignerKinetic();
        static TF2* getWignerPotential();
        static TF2* getWignerHamiltonan();
        static TF2* getWignerDeuteron();
        static TF2* getWignerDeuteronIntegral();
        static TF2* getCoalescenceProbability();

        static double getNorm();
        static double getWH();
        static double getRadius();
        static double getKStar();
        static double getRMin();
        static double getRMax();
        static double getPMin();
        static double getPMax();
        static double getMu();
        static double getRWidth();
        static double getV0();
        static double getcoal();
        static double getDeuteronInt();
        
        static double checkWxW();

    private:
        static double mHCut;
        static double mRMin;
        static double mRMax;
        static double mPMin;
        static double mPMax;
        static double mDx;
        static double mDp;
        static double mFactor;

        static TF2 *mW;
        static TF2 *mWxW;
        static TF2 *mWxJ;
        static TF2 *mWxJforItself;
        static TF2 *mK;
        static TF2 *mV;
        static TF2 *mH;
        static TF2 *mWK;
        static TF2 *mWV;
        static TF2 *mWH;
        static TF2 *mD;
        static TF2 *mDInt;
        static TF2 *mC;

        static TH2D* mDeuteronH;
        static TFile *mFileDeuteron;

        static void setThreeParam(TF2 *function);
        static void setFourParam(TF2 *function);
        static void setSixParam(TF2 *function);

        static double mR0;
        static double mRadius;
        static double mKStar;
        static double mKin;
        static double mNorm;
        static double mMu;
        static double mRWidth; //2.1
        static double mV0; //-0.0337
        static double mBoundE;

        static void normalization();
        static void reSetKStar();
        static void reSetRadius();
        static void reSetNorm();
        static void reSetMu();
        static void reSetRWidth();
        static void reSetV0();
};

#endif