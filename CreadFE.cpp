#include "CreadFE.h"

bool readerFE::LoadEvent(int iev) {
  if(vreader::LoadEvent(iev)){
    return 1;
  }

  for(int i=0; i < mNtrack; i++){

    // all stable particles in the tree
    int charge = TDatabasePDG::Instance()->GetParticle(mPID[i])->Charge();

    if(charge % 3){ // no fractional charges
      continue;
    }

    double energy = TDatabasePDG::Instance()->GetParticle(mPID[i])->Mass();
    energy = sqrt(energy*energy + mPx[i]*mPx[i] + mPy[i]*mPy[i] + mPz[i]*mPz[i]);

    particleMC part(mPx[i],mPy[i],mPz[i],energy);

    //applying mMinEta, mMaxEta cuts
    if(part.q.Eta() < mMinEta || part.q.Eta() > mMaxEta){
      continue;
    }

    // add particle
    part.pdg = mPID[i];
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

  mTree->SetBranchAddress("ntrack",&mNtrack);
  mTree->SetBranchAddress("PID",mPID);
  mTree->SetBranchAddress("px",mPx);
  mTree->SetBranchAddress("py",mPy);
  mTree->SetBranchAddress("pz",mPz);
}
