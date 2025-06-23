#include "CfemptoSource.h"
#include "CwaveUtils.h"
#include "TRandom.h"

void femptoSource::init() {
  method::init();
  vfempto::init();
}

//_________________________________________________________________________
bool femptoSource::set(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong) {
  if (utils::getKstar(p1, p2) > mThreshold) {
    return false;
  }

  int pdgM = std::abs(p1.pdg);
  int pdgL = std::abs(p2.pdg);
  if(pdgM < pdgL){
    int dummy = pdgL;
    pdgL = pdgM;
    pdgM = dummy;
  }

  utils::type system = utils::none;
  if (pdgM == 2112) {
    system = utils::nn;
  } else if (pdgM == 2212) {
    system = (pdgL == 2112) ? utils::pn : utils::pp;
  } else if (pdgM == 4324) {
    if (pdgL == 2112) system = utils::Dn;
    else if (pdgL == 2212) system = utils::Dp;
    else system = utils::DD;
  } else if (pdgM == 6436) {
    system = (pdgL == 2112) ? utils::Tn : utils::Tp;
  } else if (pdgM == 6536) {
    system = (pdgL == 2112) ? utils::Hen : utils::Hep;
  }

  if (system == utils::none) {
    return false;
  }

  setCharges(chargeStrong, chargeColoumb);

  double kt = utils::getKt(p1, p2);
  double kstar = utils::getKstar(p1, p2) * 1E3; // to MeV
  setKstar(kstar, kt, system);

  return true;
}

//_________________________________________________________________________
double femptoSource::doInteract(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong,
                                 float sumRadii, float* pos, float* posLab) {
  int pdgM = std::abs(p1.pdg);
  int pdgL = std::abs(p2.pdg);
  if(pdgM < pdgL){
    int dummy = pdgL;
    pdgL = pdgM;
    pdgM = dummy;
  }

  bool coalAllow = false;
  utils::type system = utils::none;

  if (pdgM == 2112) {
    system = utils::nn;
  } else if (pdgM == 2212) {
    if (pdgL == 2112) {
      system = utils::pn;
      coalAllow = true;
    } else {
      system = utils::pp;
    }
  } else if (pdgM == 4324) {
    if (pdgL == 2112) {
      system = utils::Dn;
      coalAllow = true;
    } else if (pdgL == 2212) {
      system = utils::Dp;
      coalAllow = true;
    } else {
      system = utils::DD;
      coalAllow = true;
    }
  } else if (pdgM == 6436) {
    if (pdgL == 2112) {
      system = utils::Tn;
    } else if (pdgL == 2212) {
      system = utils::Tp;
      coalAllow = true;
    }
  } else if (pdgM == 6536) {
    if (pdgL == 2112) {
      system = utils::Hen;
      coalAllow = true;
    } else if (pdgL == 2212) {
      system = utils::Hep;
    }
  }

  if (system == utils::none) {
    return false;
  }

  double kt = utils::getKt(p1, p2);
  setKstar(mThreshold * 1E3, kt, system);

  float coalProbAtTh = 0;
  if (std::abs(chargeStrong) > 0.9 && std::abs(chargeColoumb) < 1E-3) {
    coalProbAtTh = calcProb();
  }
  double momFinalAtTh = getKstarFinal(coalProbAtTh) * 1E-3;
  float scalingAtTh = momFinalAtTh / mThreshold;
  if (scalingAtTh > 1) scalingAtTh = 1;

  double kstar = utils::getKstar(p1, p2) * 1E3;
  setKstar(kstar, kt, system);

  TLorentzVector pSum = p1.q + p2.q;
  TVector3 b = pSum.BoostVector();
  TVector3 bInv = -b;

  p1.q.Boost(bInv);
  p2.q.Boost(bInv);

  float coalProb = 0;
  if (coalAllow) {
    coalProb = calcProb();
  }

  double momFinal = getKstarFinal(coalProb) * 1E-3 / scalingAtTh;
  float scaling = momFinal / p1.q.P();

  double m1 = p1.q.M();
  double px1 = p1.q.Px() * scaling;
  double py1 = p1.q.Py() * scaling;
  double pz1 = p1.q.Pz() * scaling;
  double e1 = sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + m1 * m1);
  double ptEx = std::abs(momFinal) - std::abs(p1.q.P());
  p1.q.SetPxPyPzE(px1, py1, pz1, e1);

  double m2 = p2.q.M();
  double px2 = p2.q.Px() * scaling;
  double py2 = p2.q.Py() * scaling;
  double pz2 = p2.q.Pz() * scaling;
  double e2 = sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + m2 * m2);
  p2.q.SetPxPyPzE(px2, py2, pz2, e2);

  p1.q.Boost(b);
  p2.q.Boost(b);

  return ptEx;
}
