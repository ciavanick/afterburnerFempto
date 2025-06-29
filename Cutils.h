#ifndef CUTILS
#define CUTILS

#include "TLorentzVector.h"
#include "TDatabasePDG.h"

class particle
{
public:
  float pt;
  float p;
  float ptpc;
  float dedx;
  float eta;
  float phi;
  float chi2;
  float trdchi2;
  float l;
  float time;
  float exppi;
  float expka;
  float exppr;
  float expde;
  float exphe;
  float sigmapi;
  float sigmaka;
  float sigmapr;
  float sigmade;
  float sigmahe;
  bool ispileup;
};

class particleMC
{
public:
  TLorentzVector q;
  int pdg = 0;
  int mother = -1;
  int StrongC = false;
  int ColoumbC = false;
  std::vector<int> daughters;

  particleMC() {}
  particleMC(double px, double py, double pz, double e) : q(px, py, pz, e) {}
  void print() const;
};

class particleCand
{
public:
  std::vector<TLorentzVector> q;
  std::vector<int> pdgOptions;

  particleCand() {}
  void addOption(double px, double py, double pz, int pdg);
};

class utils
{
  public:
    enum type{
      nn=0,
      pn=1,
      pp=2,
      Dn=3,
      Dp=4,
      DD=5,
      Tn=6,
      Tp=7,
      Hen=8,
      Hep=9,
      none=100
    };

    static double getKstar(const particleMC& p1, const particleMC& p2);                         // return kstar
    static double getKstarAsPr(const particleMC& p1, const particleMC& p2);                     // return kstar
    static double getKt(const particleMC& p1, const particleMC& p2);                            // return kT
    static double getKtAsPr(const particleMC& p1, const particleMC& p2);                        // return kT
    static double getKstar(const particleCand& p1, const particleCand& p2, int iPdg1, int iPdg2);           // return kstar
    static double getKstarAsPr(const particleCand& p1, const particleCand& p2, int iPdg1, int iPdg2);       // return kstar
    static double getKt(const particleCand& p1, const particleCand& p2, int iPdg1, int iPdg2);              // return kT
    static double getKtAsPr(const particleCand& p1, const particleCand& p2, int iPdg1, int iPdg2);          // return kT
    static double getMass(int pdg);
    static double getDPhi(const particleCand& p1, const particleCand& p2, int iPdg1, int iPdg2, double rangeMin = mRangeMin, double rangeMax = mRangeMax); //return phi
    static double getDEta(const particleCand& p1, const particleCand& p2, int iPdg1, int iPdg2); //return eta
    static void setDPhiRange(double rangeMin, double rangeMax);
  private:
    static double mRangeMin;
    static double mRangeMax;
};

#endif
