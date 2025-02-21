#include "CreadMC.h"

bool readerMC::LoadEvent(int iev) {
  if(vreader::LoadEvent(iev)){
    return 1;
  }

  return 0;
}
//_________________________________________________________________________
void readerMC::initTree(){
  if(! mIsOpened){
    printf("init tree cannot be run since file was not open\n");
    return;
  }

  static std::vector<particleMC>* point = &mParticles;
  mTree->SetBranchAddress("tracks",&point);
}
