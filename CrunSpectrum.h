#ifndef CRUNSPECTRUM
#define CRUNSPECTRUM


#include "Cvrun.h"
#include "TH1D.h"

class runSpectrum : public vrun
{
    public:
      runSpectrum(int nmix=10, TString dirName = "runSpectrum") : vrun(nmix, dirName){}
    private:
      TH1D *mHProtonSpectrum = nullptr;
      TH1D *mHNeutronSpectrum = nullptr;
      TH1D *mHDeuteronSpectrum = nullptr;
      void initHistos() override;
      void initEventsHisto() override;
      void process() override;
      void writeHistos() override;

      int selectP(const particleCand& p);

      int mPDGPr = 2212;
      int mPDGNe = 2112;
      int mPDGDe = 4324;
};

#endif