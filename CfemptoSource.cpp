#include "CfemptoSource.h"
#include "TRandom.h"

void femptoSource::init(){
  waveUtils::init();
  vfempto::init();
}
//_________________________________________________________________________
bool femptoSource::set(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong){
 if(utils::getKstar(p1,p2) > mThreshold){
    return false;
  }
  int pdgM = std::abs(p1.pdg);
  int pdgL = std::abs(p2.pdg);

  if(pdgM < pdgL){
    int dummy = pdgL;
    pdgL = pdgM;
    pdgM = dummy;
  }
  
  waveUtils::type system=waveUtils::none;
  if(pdgM == 2112){ // nn
    system = waveUtils::nn;
  } else if (pdgM == 2212) { // pn or pp
    if(pdgL == 2112){
      system = waveUtils::pn;
    } else {
      system = waveUtils::pp;
    }
  } else if(pdgM == 4324) { // Dn or Dp or DD
    if(pdgL == 2112){
      system = waveUtils::Dn;
    } else if(pdgL == 2212){
      system = waveUtils::Dp;
    } else {
      system = waveUtils::DD;
    }
  } else if(pdgM == 6436) { // Tn or Tp
    if(pdgL == 2112){
      system = waveUtils::Tn;
    } else if(pdgL == 2212){
      system = waveUtils::Tp;
    }
  } else if(pdgM == 6536) { // 3Hen o 3Hep
    if(pdgL == 2112){
      system = waveUtils::Hen;
    } else if(pdgL == 2212){
      system = waveUtils::Hep;
    }
  } 
  
  if(system == waveUtils::none){
    return false;
  }
  
  setCharges(chargeStrong,chargeColoumb);
  double kt = utils::getKt(p1,p2);
  double kstar = utils::getKstar(p1,p2) * 1E3; // to MeV

  // compute scaling at threshold
  setKstar(kstar, kt, system);

  return true;
}
//_________________________________________________________________________
double femptoSource::doInteract(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong, float sumRadii, float *pos, float *posLab){
  int pdgM = std::abs(p1.pdg);
  int pdgL = std::abs(p2.pdg);

  if(pdgM < pdgL){
    int dummy = pdgL;
    pdgL = pdgM;
    pdgM = dummy;
  }

  waveUtils::type system=waveUtils::none;
  if(pdgM == 2112){ // nn
    system = waveUtils::nn;
  } else if (pdgM == 2212) { // pn or pp
    if(pdgL == 2112){
      system = waveUtils::pn;
    } else {
      system = waveUtils::pp;
    }
  } else if(pdgM == 4324) { // Dn or Dp or DD
    if(pdgL == 2112){
      system = waveUtils::Dn;
    } else if(pdgL == 2212){
      system = waveUtils::Dp;
    } else {
      system = waveUtils::DD;
    }
  } else if(pdgM == 6436) { // Tn or Tp
    if(pdgL == 2112){
      system = waveUtils::Tn;
    } else if(pdgL == 2212){
      system = waveUtils::Tp;
    }
  } else if(pdgM == 6536) { // 3Hen o 3Hep
    if(pdgL == 2112){
      system = waveUtils::Hen;
    } else if(pdgL == 2212){
      system = waveUtils::Hep;
    }
  }

  if(system == waveUtils::none){
    return false;
  }

  // compute scaling at threshold
  double kt = utils::getKt(p1,p2);
  setKstar(mThreshold * 1E3, kt, system);

  float coalProbAtTh = 0;
  if(std::abs(chargeStrong) > 0.9 && std::abs(chargeColoumb) < 1E-3) { // coalescence allowed
    coalProbAtTh = calcProb();
  }
  double momFinalAtTh = getKstarFinal(coalProbAtTh) * 1E-3;
  float scalingAtTh = momFinalAtTh / mThreshold;

  if(scalingAtTh > 1){
    scalingAtTh = 1;
  }

  double kstar = utils::getKstar(p1,p2) * 1E3; // to MeV
  setKstar(kstar,kt,system);

  TLorentzVector pSum = p1.q + p2.q;
  TVector3 b = pSum.BoostVector();
  TVector3 bInv = -b;

  p1.q.Boost(bInv);
  p2.q.Boost(bInv);

  float coalProb = 0;
  if(std::abs(chargeStrong) > 0.9 && std::abs(chargeColoumb) < 1E-3) { // coalescence allowed
    coalProb = calcProb();
  }

  double momFinal = getKstarFinal(coalProb) * 1E-3 / scalingAtTh;
  float scaling = momFinal / p1.q.P();

//  if(scaling < 0.8){
//    printf("kt=%f - k*=%f -> scaling=%f - K = %f\n",kt,p1.q.P(),scaling,p1.q.P()*p1.q.P()/938.);
//    scaling = 0.8;
//  }

  double m1 = p1.q.M();
  double px1= p1.q.Px()*scaling;
  double py1= p1.q.Py()*scaling;
  double pz1= p1.q.Pz()*scaling;
  double e1 = sqrt(px1*px1 + py1*py1 + pz1*pz1 + m1*m1);
  double ptEx = std::abs(momFinal) - std::abs(p1.q.P());
  p1.q.SetPxPyPzE(px1,py1,pz1,e1);
  double m2 = p2.q.M();
  double px2= p2.q.Px()*scaling;
  double py2= p2.q.Py()*scaling;
  double pz2= p2.q.Pz()*scaling;
  double e2 = sqrt(px2*px2 + py2*py2 + pz2*pz2 + m2*m2);
  p2.q.SetPxPyPzE(px2,py2,pz2,e2);

  p1.q.Boost(b);
  p2.q.Boost(b);

  return ptEx;
}
