std::vector<particleMC>& readEvent();
vreader *reader;

void filtering(){
  reader = new readerFE();
  reader->openFile("lista",true); // true to read a collection
  const int nev = reader->getNevents();

  std::vector<particleMC>& event = reader->getParticles();
  std::vector<particleMC>* pevent = &event;

  TFile *fout = new TFile("output.root","RECREATE");
  TTree *tout = new TTree("tree","tree");
  tout->Branch("tracks",&pevent);

  for(int i=0; i < nev; i++){
    readEvent();
    tout->Fill();
  }

  fout->cd();
  tout->Write();
  fout->Close();
}

std::vector<particleMC>& readEvent(){
  reader->NextEvent();
  return reader->getParticles();
}
