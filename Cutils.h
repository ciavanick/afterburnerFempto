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
};

class utils
{
  public:
    static double getKstar(const particleMC& p1, const particleMC& p2);           // return kstar
    static double getKt(const particleMC& p1, const particleMC& p2);              // return kT
};

#endif
