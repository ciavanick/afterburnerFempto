#ifndef CWIGNERSOURCE
#define CWIGNERSOURCE

#include "TF2.h"
#include "Cvfempto.h"
#include "Cvread.h"

class wignerSource{
    public:
        wignerSource(TString name = "") {}
        void initFunctions(); //changed in init and removed
        void setFunctionsParameters(); //removed 
        void setRadius(double radius);
        void setR0(double r0);
        void setRadiusK(double r0);
        void setKstar(double k);
        void setKIn(double k);
        void setRanges(double xmin, double ymin, double xmax, double ymax);
        void setMu(double mu);
        void setRWidth(double rWidth);
        void setV0(double v0);

        double doInteract(particleMC& p1, particleMC& p2);
        float getCoalProb(const particleMC& p1, const particleMC& p2);

        TF2* getWignerFunction();
        TF2* getWignerFunctionForItself();
        TF2* getWignerFunctionForJacobian();
        TF2* getWignerFunction2ForJacobian();
        TF2* getKineticEnergyFunction();
        TF2* getPotentialEnergyFunction();
        TF2* getHamiltonianFunction();
        TF2* getWignerKinetic();
        TF2* getWignerPotential();
        TF2* getWignerHamiltonan();
        TF2* getWignerDeuteron();
        TF2* getWignerDeuteronIntegral();
        TF2* getCoalescenceProbability();

        double getNorm();
        double getwK();
        double getwV();
        double getwH();
        double getRadius();
        double getKStar();
        double getRMin();
        double getRMax();
        double getPMin();
        double getPMax();
        double getMu();
        double getRWidth();
        double getV0();
        double getcoal();
        double getDeuteronInt();
        
        double checkWxW();

    private:

        void normalization(); //removed
        void reSetNorm();
        void reSetRadius();
        void reSetKStar();
        void reSetMu();
        void reSetRWidth();
        void reSetV0();

};


#endif