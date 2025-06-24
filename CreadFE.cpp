#include "CreadFE.h"

bool readerFE::LoadEvent(int iev) {
  if(vreader::LoadEvent(iev)){
    return 1;
  }
  int ntrack = mArr->GetEntriesFast();

  mParticles.clear();

  for(int i=0; i < ntrack; i++){
    TParticle *p = (TParticle*) mArr->At(i);

    int pdg = p->GetPdgCode();
    if(pdg == 2212000){
      pdg = 2212;
    }
    if(pdg == 2112000){
      pdg = 2112;
    }
    if(pdg == -2212000){
      pdg = -2212;
    }
    if(pdg == -2112000){
      pdg = -2112;
    }
    if(pdg == 1000010020 || pdg == -1000010020){ //deuteron
      continue;
    }
    if(! TDatabasePDG::Instance()->GetParticle(pdg)){
      printf("%d doesn't exist\n",pdg);
      continue;
    }

    // all stable particles in the tree
    int charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();

    double px = p->Px();
    double py = p->Py();
    double pz = p->Pz();

    double energy = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
    energy = sqrt(energy*energy + px*px + py*py + pz*pz);

    particleMC part(px,py,pz,energy);

    //applying mMinEta, mMaxEta cuts
    if(part.q.Eta() < mMinEta || part.q.Eta() > mMaxEta){
      continue;
    }

    // add particle
    part.pdg = pdg;
    part.mother = -1;
    if(part.pdg > 999){
      part.StrongC = 1;
    } else if(part.pdg < -999){
      part.StrongC = -1;
    } else {
      part.StrongC = 0;
    }
    part.ColoumbC = TDatabasePDG::Instance()->GetParticle(part.pdg)->Charge()/3;
//  part.daughters
    if(std::abs(part.pdg) > 999 && std::abs(part.pdg) < 9999){//== 2212 || std::abs(part.pdg) == 2112){
      mParticles.push_back(part);
    }
  }

  return 0;
}
//_________________________________________________________________________
void readerFE::initTree(){
  if(! mIsOpened){
    printf("init tree cannot be run since file was not open\n");
    return;
  }

  mArr = new TClonesArray("TParticle");

  mTree->SetBranchAddress("Particles",&mArr);
}
