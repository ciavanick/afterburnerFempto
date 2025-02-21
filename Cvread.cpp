#include "Cvread.h"

void vreader::openFile(const char *fname, bool fromCollection){
  if(mIsOpened){
    printf("file reader already set, clear it if you want to reuse it\n");
  }

  if(! fname) {
    fname = mFileName.Data();
  }

  mTree = new TChain(getTreeName());
  if(fromCollection){
    FILE *fin = fopen(fname,"r");
    char namefile[1000];
    while(fscanf(fin,"%s",namefile) == 1){
      mTree->AddFile(namefile);
    }
    fclose(fin);
  } else {
    mTree->AddFile(fname);
  }

  if(mTree->GetEntries() < 1){
    printf("no events found\n");
    return;
  }

  mIsOpened = true;

  initTree();
}
//_________________________________________________________________________
void vreader::clear(){
  mIsOpened = false;
  mIEvent = -1;
}
//_________________________________________________________________________
bool vreader::LoadEvent(int iev){
  if(! mTree){
    printf("tree not set");
    return 1;
  }

  if(iev >= mTree->GetEntries()){
    printf("event %d >= %llu\n",iev,mTree->GetEntries());
    return 1;
  }

  mTree->GetEvent(iev);
  mParticles.clear();

  return 0;
}
//_________________________________________________________________________
int vreader::getNevents() const {
  if(mTree){
    return mTree->GetEntries();
  }

  return 0;
}
