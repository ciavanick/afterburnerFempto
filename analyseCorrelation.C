int pdgPr = 2212; // proton
int pdgNe = 2112; // neutron
int pdgDe = 4324; // deuteron

void analyseCorrelation(TString listname,  TString outfolder, TString outname)
{
  // set file reader
  vreader *reader = new readerFE();
  reader->openFile(listname, true); // true to read a collection
  int nev = reader->getNevents();
  printf("N events found = %d\n", nev);

  // set analisys
  int nmixedEvents = 10;
  vrun *runPrPr = new vrun(nmixedEvents);
  runPrPr->setIDName("PrPr");
  runPrPr->selectPDG(pdgPr, pdgPr); // pdg of the two particles to be correlated
  runPrPr->init();                  // init analysis

  vrun *runPrNe = new vrun(nmixedEvents);
  runPrNe->setIDName("PrNe");
  runPrNe->selectPDG(pdgPr, pdgNe); // pdg of the two particles to be correlated
  runPrNe->init();                  // init analysis

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

    // then analyze
    runPrPr->setEvent(event); // pass the event to the analyzer
    runPrPr->doAnalysis();       // process the event

    runPrNe->setEvent(event); // pass the event to the analyzer
    runPrNe->doAnalysis();       // process the event
  }

  // finalization
  runPrPr->finalize(); // finalize the analysis
  TFile *foutPrPr = new TFile(outfolder + "pp_" + outname, "RECREATE");
  runPrPr->write(); // write outputs of the analysis
  //hCorrelationPrPr->Write();
  foutPrPr->Close();

  // finalization
  runPrNe->finalize(); // finalize the analysis
  TFile *foutPrNe = new TFile(outfolder + "pn_" + outname, "RECREATE");
  runPrNe->write(); // write outputs of the analysis
  //hCorrelationPrNe->Write();
  foutPrNe->Close();
}
