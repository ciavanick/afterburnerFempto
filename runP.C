int pdgPr = 2212; // proton
int pdgNe = 2112; // neutron
int pdgDe = 4324; // deuteron
int pdgAntiPr = -2212; // antiproton
int pdgAntiNe = -2112; // antineutron
int pdgAntiDe = -4324; // antideuteron

bool isWigner = false;

#define MAXD(...) std::max(std::initializer_list<double>{__VA_ARGS__})

void runP(TString listname, TString outfolder, TString outname, double R_D = 3.2, double R_0 = 1.125)
{
  // set file reader
  vreader *reader = new readerFE();
  //reader->setEtaRange(-100.,100);
  reader->openFile(listname, true); // true to read a collection
  int nev = reader->getNevents();
  printf("N events found = %d\n", nev);

  // set analisys
  int nmixedEvents = 10;

  runPhotons *photons = new runPhotons(); //spectrum analysis
  //photons->setRapidityRange(-100.,100);
  photons->init();

  //waveUtils::setParams(10, R);

  //waveUtils::setParams(waveUtils::calculateV0(), R);

  //-------- interactor ------------------------
  double vStrong = 0;
  vfempto *interactor;
  auto R_0_eff = R_0*2*(1); 
  if(! isWigner){
   interactor = new femptoSource;      // Lenard-Jones  strong potential
   interactor->setParams(17.4, R_D, 1.44, R_0_eff, 3./8);
   vStrong = waveUtils::calculateV0();
   interactor->setParams(vStrong, R_D, 1.44, R_0_eff, 3./8);
  } else {
   interactor = new femptoSource<wignerUtils>;      // Lenard-Jones trong potential
   interactor->setParams(17.4E-3, 3.2, 1.44E-3, 1.5, 3./8);
  }
  interactor->setThreshold(0.4);

  //-------- interactor ------------------------

  std::cout<<"Strong Radius: "<< R_D << "\t Strong Potential: "<< vStrong <<"\n";

  for (int i = 0; i < nev; i++)
  {
     // read event
    reader->NextEvent();
    std::vector<particleMC> &event = reader->getParticles();

    if (i % 100000 == 0)
      std::cout << "Event: " << i << "/" << nev << " - " << ((double)i / (double)nev) * 100 << "\n";

    // for MC add it the effect of afterburner
    interactor->doInteractAll(event, false, true);
    //---------------------------------------

    // then analyse
    photons->setEvent(event);
    photons->doAnalysis();
  }

  // finalization
  TFile *foutPhotons = new TFile(outfolder + "Photons" + outname , "RECREATE");
  photons->write();
  foutPhotons->Close();

}
