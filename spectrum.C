//Macro to run the pT spectrum analysis

int pdgPr = 2212; // proton
int pdgNe = 2112; // neutron
int pdgDe = 4324; // deuteron

void spectrum(TString listname,  TString outfolder, TString outname)
{
  // set file reader
  vreader *reader = new readerFE();
  reader->openFile(listname, true); // true to read a collection
  int nev = reader->getNevents();
  printf("N events found = %d\n", nev);


  runSpectrum *spectra = new runSpectrum();
  spectra->init();

  vfempto *interactor = new femptoSource();
  // interactor->setParams(33.7, 2.1, 1.44, -1.5, 0.);
  interactor->setThreshold(0.4);

  for (int i = 0; i < nev; i++)
  {
    // read event
    reader->NextEvent();
    std::vector<particleMC> &event = reader->getParticles();

    if (i % 100000 == 0)
      std::cout << "Event: " << i << "/" << nev << " - " << ((double)i / (double)nev) * 100 << "\n";

    // for MC add it the effect of afterburner
    interactor->doInteractAll(event, true, true);
    //---------------------------------------

    spectra->setEvent(event);
    spectra->doAnalysis();
  }

  TFile *foutSpec = new TFile(outfolder + "Spectrum" + outname , "RECREATE");
  spectra->write();
  foutSpec->Close();
}
