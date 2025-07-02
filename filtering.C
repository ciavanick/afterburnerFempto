std::vector<particleMC>& readEvent();
vreader *reader;

float sourceSize = 1.5;
float spinFactor=1;//3./8;

void filtering(){
  vfempto *interactor;
  interactor = new femptoSource;      // Lenard-Jones  strong potential
  interactor->setParams(17.4, 3.2, 1.44, -2*sourceSize, spinFactor);
  interactor->setThreshold(0.4);

  reader = new readerFE();
  reader->setEtaRange(-100,100);
  reader->openFile("lista",true); // true to read a collection
  const int nev = reader->getNevents();

  std::vector<particleMC>& eventMC = reader->getParticles();
  std::vector<particleMC>& event = eventMC;
  std::vector<particleMC>* pevent = &event;

  TFile *fout = new TFile("output.root","RECREATE");
  TTree *tout = new TTree("tree","tree");
  tout->Branch("tracks",&pevent);

  for(int i=0; i < nev; i++){
    event.clear();
    readEvent();
    interactor->doInteractAll(eventMC);

    bool write=false;

    for(auto& o : eventMC){
      if(std::abs(o.pdg) > 3000){
        write = true;
      }
    }

    if(write){
      tout->Fill();
    }
  }

  fout->cd();
  tout->Write();
  fout->Close();
}

std::vector<particleMC>& readEvent(){
  reader->NextEvent();
  return reader->getParticles();
}
