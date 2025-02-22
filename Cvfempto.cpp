#include "Cvfempto.h"
#include "TRandom.h"

void vfempto::doInteractAll(std::vector<particleMC>& part){
  // build pairs only if both charged
  std::vector<std::pair<int,int>> intPairs;
  for(int i1 = 0; i1 < part.size(); i1++){
    int sC = part[i1].StrongC;
    int cC = part[i1].ColoumbC;

    if(part[i1].mother > -1){ // only mother inteacts
      continue;
    }

    if(!(sC || cC)){ // particle can interact otherwise skip
      continue;
    }

    for(int i2 = i1+1; i2 < part.size(); i2++){
      if(part[i2].mother > -1){ // only mother inteacts
        continue;
      }

      if(!(part[i2].StrongC && sC || part[i2].ColoumbC && cC)){ // two particles can interact otherwise skip
        continue;
      }

      double kstar = utils::getKstar(part[i1],part[i2]);
      if(kstar > 0.4){
        continue;
      }

      intPairs.push_back(std::make_pair(i1,i2));
    }
  }

  // sorting by kstar
  std::sort(intPairs.begin(), intPairs.end(), [part](const auto &a, const auto &b)
  {
     double ks1 = utils::getKstar(part[a.first],part[a.second]);
     double ks2 = utils::getKstar(part[b.first],part[b.second]);
     return ks1 < ks2;
  });

  for(const auto& o : intPairs){ // do interactions starting from lower kstar
     particleMC original1 = part[o.first];
     particleMC original2 = part[o.second];

     doInteract(part[o.first],part[o.second],part[o.first].ColoumbC*part[o.second].ColoumbC,part[o.first].StrongC*part[o.second].StrongC,0);

     // propagate effect of the inteaction to daugher particles
     for(const auto& id : original1.daughters){
       particleMC& dau = part[id];
       TVector3 bF = original1.q.BoostVector();       // before interaction
       TVector3 bFinv = -bF;
       TVector3 bA = part[o.first].q.BoostVector();   // after interaction

       dau.q.Boost(bFinv);
       dau.q.Boost(bA);
     }
     for(const auto& id : original2.daughters){
       particleMC& dau = part[id];
       TVector3 bF = original2.q.BoostVector();       // before interaction
       TVector3 bFinv = -bF;
       TVector3 bA = part[o.second].q.BoostVector();  // after interaction

       dau.q.Boost(bFinv);
       dau.q.Boost(bA);
     }
  }
}
//_________________________________________________________________________
double vfempto::doInteract(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong, float sumRadii, float *pos, float *posLab){
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

  float sumRadiiOr = sumRadii;

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
  deltaE -= chargeStrong*mStrong*(dist < mStrongR); // Strong potential attractive

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
