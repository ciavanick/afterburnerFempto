#include "Cutils.h"

double utils::getKstar(const particleMC& p1, const particleMC& p2){
  TLorentzVector pSum = p1.q + p2.q;
  double Minv = pSum.M();
  double mass1square = p1.q.M()*p1.q.M();
  double mass2square = p2.q.M()*p2.q.M();
  double A = Minv*Minv - mass1square - mass2square;

  return sqrt(A*A - 4 * mass1square * mass2square) * 0.5 / Minv;
}
//_________________________________________________________________________
double utils::getKt(const particleMC& p1, const particleMC& p2){
  TLorentzVector pSum = 0.5*(p1.q + p2.q);
  return pSum.Pt();
}
//_________________________________________________________________________
double utils::getKstarAsPr(const particleMC& p1, const particleMC& p2){
  float m1scaling = 0.938/p1.q.M();
  float m2scaling = 0.938/p2.q.M();
  TLorentzVector pSum = p1.q*m1scaling + p2.q*m2scaling;
  double Minv = pSum.M();
  double mass1square = p1.q.M()*p1.q.M();
  double mass2square = p2.q.M()*p2.q.M();
  double A = Minv*Minv - mass1square - mass2square;

  return sqrt(A*A - 4 * mass1square * mass2square) * 0.5 / Minv;
}
//_________________________________________________________________________
double utils::getKtAsPr(const particleMC& p1, const particleMC& p2){
  float m1scaling = 0.938/p1.q.M();
  float m2scaling = 0.938/p2.q.M();
  TLorentzVector pSum = 0.5*(p1.q*m1scaling + p2.q*m2scaling);
  return pSum.Pt();
}

