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
