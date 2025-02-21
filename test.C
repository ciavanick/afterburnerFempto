std::vector<particleMC>&doToyEvent();
std::vector<particleMC>& readEvent();
bool sel1(const particleMC& p);
bool sel2(const particleMC& p);
bool isproton(const particleMC& p);

float meson_downscale = 0.1;
bool fromFile = false;

vreader *reader;

void test(bool interact=true, bool coalescence=true, bool noHe=false){
  if(fromFile){
    reader = new readerMC();
    reader->openFile("listaMC",true); // true to read a collection
  }

  bool samePart = true;

  const int nbins = 20;
  const int nev = fromFile ? reader->getNevents() : 1E5;
  const int nmix = 5;

  int chargeComb = 1; // 1=same charge, -1=opposite charge
  if(chargeComb < 0){
    samePart = false;
  }

  float sumRadii = 0;

  float sourceSize = 0;

  fempto interactor;
  vcoal *coalescer = new simpleCoal();

  if(noHe){
    coalescer->setMassMax(2);
  }

  interactor.setParams(2.2E-3, 2.4, 1.44E-3, sourceSize);

  TH1D *h = new TH1D("h",";k* (GeV/c)",nbins,0,0.2);
  TH1D *h2 = new TH1D("h2",";k* (GeV/c)",nbins,0,0.2);
  TH1D *hb = new TH1D("hb",";k* (GeV/c)",nbins,0,0.2);

  TH1D *hPr = new TH1D("hPr",";p/A (GeV/c)",100,0,5);
  TH1D *hNe = new TH1D("hNe",";p/A (GeV/c)",100,0,5);
  TH1D *hDe = new TH1D("hDe",";p/A (GeV/c)",100,0,5);
  TH1D *hHe = new TH1D("hHe",";p/A (GeV/c)",100,0,5);

  TH1D *hB2 = new TH1D("hB2",";p/A (GeV/c);N(p/A)_{D} / N(p/A)_{p} / N(p/A)_{n} (GeV/c)",100,0,5);
  TH1D *hB3 = new TH1D("hB3",";p/A (GeV/c);N(p/A)_{D} / N(p/A)_{p} / N(p/A)_{n} (GeV/c)",100,0,5);

  TH1D *hMass = new TH1D("hMass",";mass (GeV/c^{2})",1000,0,10);

  h->SetLineColor(2);
  h2->SetLineColor(1);
  hb->SetLineColor(4);

  std::vector<particleMC> evPrev[nmix];

  int nDeu=0;
  int nAntiDeu=0;

  for(int i=1; i < nmix; i++){
    if(fromFile){
      evPrev[i] = readEvent();
    } else {
      evPrev[i] = doToyEvent();
    }
    if(coalescence){
      coalescer->doCoalAll(evPrev[i]);
    }

    if(interact){
      interactor.doInteractAll(evPrev[i]);
    }
  }

  printf("start running events \n");
  for(int i=0; i < nev-nmix; i++){
    if(! (i%10000)) printf("%d/%d\n",i,nev);
    int lev = i % nmix;
    std::vector<particleMC>& event = fromFile ? readEvent() : doToyEvent();

    if(coalescence){
      coalescer->doCoalAll(event);
    }

    if(interact){
      interactor.doInteractAll(event);
    }

    for(int i1=0; i1 < event.size(); i1++){ // same event
      const particleMC& p1 = event[i1];

      if(p1.daughters.size()){
        continue;
      }

      hMass->Fill(p1.q.M());

      if(p1.pdg == 6536){
        hHe->Fill(p1.q.P()/3);
      } else if(p1.pdg == 4324){
        nDeu++;
        hDe->Fill(p1.q.P()*0.5);
      }
      if(p1.pdg == -6536){
        hHe->Fill(p1.q.P()/3);
      } else if(p1.pdg == -4324){
        nAntiDeu++;
        hDe->Fill(p1.q.P()*0.5);
      }

      if(std::abs(p1.pdg) == 2212){
        hPr->Fill(p1.q.P());
      }
      if(std::abs(p1.pdg) == 2112){
        hNe->Fill(p1.q.P());
      }

      if(! isproton(p1)){
        continue;
      }
      for(int i2=(i1+1)*samePart; i2 < event.size(); i2++){
        const particleMC& p2 = event[i2];
        if(p2.daughters.size()){
          continue;
        }
        if((!isproton(p2)) || chargeComb*p1.ColoumbC*p2.ColoumbC < 0){
          continue;
        }

        double kt = utils::getKt(p1,p2);
        if(fabs(kt-1) < 0.2) {
          double kstar = utils::getKstar(p1,p2);
          h->Fill(kstar);
          h2->Fill(kstar);
        }
      }
    }

    for(int j=0; j < nmix;j++){ // mixed event
      if(j==lev){
        evPrev[j] = event;
        continue;
      }
      for(int i1=0; i1 < event.size(); i1++){
        const particleMC& p1 = event[i1];
        if(p1.daughters.size()){
          continue;
        }

        if(! isproton(p1)){
          continue;
        }

        for(int i2=0; i2 < evPrev[j] .size(); i2++){
          const particleMC& p2 = evPrev[j][i2];
          if(p2.daughters.size()){
            continue;
          }

          if((!isproton(p2)) || chargeComb*p1.ColoumbC*p2.ColoumbC < 0){
            continue;
          }

          double kt = utils::getKt(p1,p2);
          if(fabs(kt-1) < 0.2) {
            double kstar = utils::getKstar(p1,p2);
            hb->Fill(kstar);
          }
        }
      }
    }
  }

  printf("nDeu = %d - nAntiDeu = %d\n",nDeu, nAntiDeu);

  float adjustScaling = h->Integral(nbins/2+1,nbins) / hb->Integral(nbins/2+1,nbins);
  hb->Scale(adjustScaling);

  h->Sumw2();
  h2->Sumw2();
  h2->Divide(hb);

  hb->Scale(1./adjustScaling); //remove scaling before to store

  TCanvas *c = new TCanvas;
  h2->Draw();

  new TCanvas;
  hPr->Draw();
  hNe->Draw("SAME");
  hDe->Draw("SAME");
  hHe->Draw("SAME");
  hPr->SetLineColor(1);
  hNe->SetLineColor(4);
  hDe->SetLineColor(2);
  hHe->SetLineColor(6);
  hPr->SetMinimum(1);

  hDe->Sumw2();
  hPr->Sumw2();
  hNe->Sumw2();
  hHe->Sumw2();

  hB2->Divide(hDe,hPr);
  hB2->Divide(hB2,hNe);
  hB2->Scale(2.*nev);
  hB3->Divide(hHe,hPr);
  hB3->Divide(hB3,hPr);
  hB3->Divide(hB3,hNe);
  hB3->Scale(4.*nev*nev);

  new TCanvas;
  hB2->Draw();
  hB3->Draw("SAME");

  c = new TCanvas;
  c->SetLogy();
  hMass->Draw();

  TFile *fout = new TFile("resNew.root","RECREATE");
  h->Write();
  h2->Write();
  hb->Write();
  hPr->Write();
  hNe->Write();
  hDe->Write();
  hHe->Write();
  hMass->Write();
  fout->Close();
}

std::vector<particleMC>& readEvent(){
  reader->NextEvent();
  return reader->getParticles();
}

std::vector<particleMC>& doToyEvent(){
  static std::vector<particleMC> part;

  part.clear();

  int np = gRandom->Gaus(15,3);

  for(int i=0; i < np; i++){
    float probPart = gRandom->Rndm();
    float mass;
    int pdg;
    if(probPart < 0.6*meson_downscale){
      pdg = 211;
      mass = 0.139;
    } else if(probPart < 0.8*meson_downscale){
      pdg = 321;
      mass = 0.493;
    } else {
      if(gRandom->Rndm() < 0.5){
        pdg = 2212;
        mass = 0.938;
      } else {
        pdg = 2112;
        mass = 0.940;
      }
    }

    double phi1 = gRandom->Rndm()*2*TMath::Pi();
    double eta1 = (gRandom->Rndm()-0.5);
    double p1 = - TMath::Log(gRandom->Rndm());
    double pt1 = p1/TMath::CosH(eta1);
    double mt1 = sqrt(mass*mass + pt1*pt1);
    double pz1 = mt1 * TMath::SinH(eta1);

    particleMC p(pt1*cos(phi1), pt1*sin(phi1),pz1,sqrt(mt1*mt1 + pz1*pz1));
    p.mother = -1;
    p.pdg = pdg;
    if(pdg > 999){
      p.StrongC = 1;
    } else {
      p.StrongC = 0;
    }

    p.ColoumbC = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/3;

    if(gRandom->Rndm() > 0.5){ // anti-particle
      p.pdg *= -1;
      p.StrongC *= -1;
      p.ColoumbC *= -1;
    }
    part.push_back(p);
  }
  return part;
}

bool sel1(const particleMC& p){
  return isproton(p);
}

bool sel2(const particleMC& p){
  return isproton(p);
}

bool isproton(const particleMC& p){
  if(std::abs(p.pdg) == 2212){
    return true;
  }
  return false;
}

