#include "Cfempto.h"
#include "TRandom.h"

double fempto::doInteract(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong, float sumRadii, float *pos, float *posLab){
  double kstar = utils::getKstar(p1,p2);

  float *lpos,*lposLab; // 3+3 position vectors

  if(pos){ // return also position in the output (working on it)
    lpos = pos;
  } else { // working on a local variable
    lpos = new float[8];
  }
  if(posLab){ // return also lab position in the output (working on it)
    lposLab = posLab;
  } else { // working on a local variable
    lposLab = new float[8];
  }

  TVector3 impactparam(0,0,0);

  float sumRadiiOr = TMath::Power(p1.q.M(),1./3) + TMath::Power(p2.q.M(),1./3);//sumRadii;

  sumRadii = 0.1973 * 0.5 / kstar;

  if(sumRadii < sumRadiiOr){
    sumRadii = sumRadiiOr;
  }

  float extraSmear = gRandom->Gaus(0,mSourceRadius);

  lposLab[0] = sumRadii * 0.5;
  lposLab[4] = -lposLab[0];
  lposLab[1] = extraSmear*0.5;
  lposLab[5] = -lposLab[1];
  lposLab[2] = lposLab[3] = lposLab[6] = lposLab[7] = 0;

  double scaling;

  TLorentzVector pSum = p1.q + p2.q;
  TVector3 b = pSum.BoostVector();
  TVector3 bInv = -b;

  p1.q.Boost(bInv);
  p2.q.Boost(bInv);

  double beta1 = p1.q.Beta();
  double beta2 = p2.q.Beta();

  TLorentzVector v1(lposLab[0],lposLab[1],lposLab[2],lposLab[3]);
  v1.Boost(bInv);
  TLorentzVector v2(lposLab[4],lposLab[5],lposLab[6],lposLab[7]);
  v2.Boost(bInv);
  lpos[0]=lposLab[0];
  lpos[1]=lposLab[1];
  lpos[2]=lposLab[2];
  lpos[3]=lposLab[3];
  lpos[4]=lposLab[4];
  lpos[5]=lposLab[5];
  lpos[6]=lposLab[6];
  lpos[7]=lposLab[7];
  impactparam.SetXYZ(lpos[0]-lpos[4], lpos[1]-lpos[5], lpos[2]-lpos[6]);

  impactparam.SetXYZ(lpos[0]-lpos[4], lpos[1]-lpos[5], lpos[2]-lpos[6]);

  float dist = sqrt(impactparam.Mag2());
  float deltaE = 0;
  deltaE = chargeColoumb * mCoulomb / dist;         // Columb contribution (repulsive/attractive)
//  deltaE -= chargeStrong*mStrong*(dist < mStrongR); // Strong potential attractive

  float rad6 = mStrongR*sqrt(sumRadiiOr*0.5)/dist;
  rad6 = rad6*rad6; // ^2
  rad6 = rad6*rad6*rad6; // ^6
  float rad12 = rad6*rad6; // ^12
  deltaE += chargeStrong*TMath::Min(0.,mStrong*double(-rad6 + 2*rad12)); // Lenard-Jones limited to +200 MeV

  double E0 = p1.q.Energy()+p2.q.Energy();
  double M0 = p1.q.M() + p2.q.M();
  double EK0 = E0 - M0;
  double EKF = EK0 + deltaE;

  if(EKF < 0){
    EKF = 0;
  }
  double EF = M0 + EKF;

  double momFinal = sqrt(EF*EF - M0*M0)*0.5;
  scaling = momFinal / p1.q.P();

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

  if(!pos){ // clean up local
    delete lpos;
  }
  if(!posLab){ // clean up local
    delete lposLab;
  }

  return ptEx;
}
