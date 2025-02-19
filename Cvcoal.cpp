#include "Cvcoal.h"

void vcoal::doCoalAll(std::vector<particleMC>& part){
  // build pairs only if one charged and one neutral (both particles or antiparticles)
  std::vector<std::pair<int,int>> intPairs;
  int initSize = 0;
  while(part.size() > initSize){
    int alreadydone = initSize;
    initSize = part.size();
    for(int i1 = 0; i1 < initSize; i1++){
      int cC = part[i1].ColoumbC;

      if(part[i1].daughters.size()){
        continue;
      }

      for(int i2 = TMath::Max(i1+1,alreadydone); i2 < part.size(); i2++){
        if(part[i1].daughters.size()){
          continue;
        }

        int cC2 = part[i2].ColoumbC;

        if(cC * cC2 != 0){ // one neutral particle expected
          continue;
        }

        if(part[i1].pdg * part[i2].pdg < 0){ // one particle and one antiparticle -> skip
          continue;
        }

        double kstar = utils::getKstar(part[i1],part[i2]);
        if(kstar > 0.2){
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
      if(docoal(part[o.first],part[o.second])){
        particleMC merged = merge(part[o.first],part[o.second]);
        if(merged.q.M() > mMassMax){ // skip it
          continue;
        }

        int pos = part.size();

        part[o.first].daughters.push_back(pos);
        part[o.second].daughters.push_back(pos);
        part.push_back(merged);
      }
    }
  }
}
//_________________________________________________________________________
particleMC vcoal::merge(const particleMC& p1, const particleMC& p2){
  particleMC pSum;
  pSum.q = p1.q + p2.q;
  int pdgSum = p1.pdg + p2.pdg;
  double mass = p1.q.M() + p2.q.M();
  double p = pSum.q.P();
  double e = sqrt(p*p + mass*mass);
  pSum.q.SetE(e);

  pSum.pdg = p1.pdg + p2.pdg;
  pSum.mother = -1;              // not tracing mothers
  pSum.StrongC = p1.StrongC + p2.StrongC;
  pSum.ColoumbC = p1.ColoumbC + p2.ColoumbC;

//  printf("merged = %d - %f + %f = %f\n",pSum.pdg,p1.q.P(),p2.q.P(),pSum.q.P());

  return pSum;
}
