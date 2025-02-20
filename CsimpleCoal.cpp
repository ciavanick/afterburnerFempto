#include "CsimpleCoal.h"

bool simpleCoal::docoal(const particleMC& p1, const particleMC& p2){ // perform coalescence
  float kstar = utils::getKstar(p1,p2);

  if(p1.daughters.size() || p1.daughters.size()){
    return false;
  }

  if(std::abs(p1.pdg) < 999 || std::abs(p2.pdg) < 999){
    return false;
  }

  if(p1.pdg * p2.pdg < 0){
    return false;
  }

  int charge = std::abs(p1.ColoumbC + p2.ColoumbC);
  float mass = p1.q.M() + p2.q.M();

  if(charge < 1){ // no nuclei with Z=0
    return false;
  }

  if(charge == 1 && mass > 3){ // only D,T allowed
    return false;
  }

  if(charge == 2 && (mass < 2 || mass > 4)){ // He: A=3,4 allowed
    return false;
  }

  if(charge == 3 && (mass < 5 || mass > 7)){ // Li: A=6,7 allowed
    return false;
  }

  if(charge > 3){
    return false;
  }

  if(p1.ColoumbC != 0 || p2.ColoumbC != 0){

  }

  if(kstar < 0.05){
//    printf("%d %d \n",p1.pdg,p2.pdg);
    return true;
  }

  return false;
}
//_________________________________________________________________________
void simpleCoal::setParams(float *params){

}
//_________________________________________________________________________
float simpleCoal::getParams(int ipar) const {
  return 0.;
}
