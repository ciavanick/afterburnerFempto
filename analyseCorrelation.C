int pdgPr = 2212; // proton
int pdgNe = 2112; // neutron
int pdgDe = 4324; // deuteron
int pdgAntiPr = -2212; // antiproton
int pdgAntiNe = -2112; // antineutron
int pdgAntiDe = -4324; // antideuteron

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

  vrun *runPrDe = new vrun(nmixedEvents);
  runPrDe->setIDName("PrDe");
  runPrDe->selectPDG(pdgPr, pdgDe); // pdg of the two particles to be correlated
  runPrDe->init();                  // init analysis

  vrun *runAntiPrPr = new vrun(nmixedEvents);
  runAntiPrPr->setIDName("AntiPrPr");
  runAntiPrPr->selectPDG(pdgAntiPr, pdgAntiPr); // pdg of the two particles to be correlated
  runAntiPrPr->init();                  // init analysis

  vrun *runAntiPrNe = new vrun(nmixedEvents);
  runAntiPrNe->setIDName("AntiPrNe");
  runAntiPrNe->selectPDG(pdgAntiPr, pdgAntiNe); // pdg of the two particles to be correlated
  runAntiPrNe->init();                  // init analysis

  vrun *runAntiPrDe = new vrun(nmixedEvents);
  runAntiPrDe->setIDName("AntiPrDe");
  runAntiPrDe->selectPDG(pdgAntiPr, pdgAntiDe); // pdg of the two particles to be correlated
  runAntiPrDe->init();                  // init analysis

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

    runPrDe->setEvent(event); // pass the event to the analyzer
    runPrDe->doAnalysis();       // process the event

    runAntiPrPr->setEvent(event); // pass the event to the analyzer
    runAntiPrPr->doAnalysis();       // process the event

    runAntiPrNe->setEvent(event); // pass the event to the analyzer
    runAntiPrNe->doAnalysis();       // process the event

    runAntiPrDe->setEvent(event); // pass the event to the analyzer
    runAntiPrDe->doAnalysis();       // process the event
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

  runPrDe->finalize(); // finalize the analysis
  TFile *foutPrDe = new TFile(outfolder + "pd_" + outname, "RECREATE");
  runPrDe->write(); // write outputs of the analysis
  //hCorrelationPrNe->Write();
  foutPrDe->Close();

  runAntiPrPr->finalize(); // finalize the analysis
  TFile *foutAntiPrPr = new TFile(outfolder + "Antipp_" + outname, "RECREATE");
  runAntiPrPr->write(); // write outputs of the analysis
  //hCorrelationPrPr->Write();
  foutAntiPrPr->Close();

  // finalization
  runAntiPrNe->finalize(); // finalize the analysis
  TFile *foutAntiPrNe = new TFile(outfolder + "Antipn_" + outname, "RECREATE");
  runAntiPrNe->write(); // write outputs of the analysis
  //hCorrelationPrNe->Write();
  foutAntiPrNe->Close();

  runAntiPrDe->finalize(); // finalize the analysis
  TFile *foutAntiPrDe = new TFile(outfolder + "Antipd_" + outname, "RECREATE");
  runAntiPrDe->write(); // write outputs of the analysis
  //hCorrelationPrNe->Write();
  foutAntiPrDe->Close();
}
