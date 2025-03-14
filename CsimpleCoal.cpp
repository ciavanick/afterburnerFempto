#include "CsimpleCoal.h"
#include "TRandom.h"

bool simpleCoal::docoal(const particleMC& p1, const particleMC& p2){ // perform coalescence
  float kstar = utils::getKstar(p1,p2);

  if(p1.daughters.size() || p1.daughters.size()){
    return false;
  }

  int mergedPdg = p1.pdg+p2.pdg;

  if(std::abs(mergedPdg) != 4324 && std::abs(mergedPdg) != 6436 && std::abs(mergedPdg) != 6536 && std::abs(mergedPdg) != 8648 && std::abs(mergedPdg) != 12972 && std::abs(mergedPdg) != 15084){
    return false;
  }

  if(p1.pdg * p2.pdg < 0){
    return false;
  }

  int charge = std::abs(p1.ColoumbC + p2.ColoumbC);
  float mass = p1.q.M() + p2.q.M();

  float prob = TMath::Exp(-kstar/0.05);

  if(gRandom->Rndm() < prob){
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
