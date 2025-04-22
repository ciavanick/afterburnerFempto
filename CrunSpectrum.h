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
      TH1D *mHTritiumSpectrum = nullptr;
      TH1D *mHHelium3Spectrum = nullptr;
      TH1D *mHHelium4Spectrum = nullptr;
      TH1D *mHProtonSpectrumA = nullptr;
      TH1D *mHNeutronSpectrumA = nullptr;
      TH1D *mHDeuteronSpectrumA = nullptr;
      TH1D *mHTritiumSpectrumA = nullptr;
      TH1D *mHHelium3SpectrumA = nullptr;
      TH1D *mHHelium4SpectrumA = nullptr;

      void initHistos() override;
      void initEventsHisto() override;
      void process() override;
      void writeHistos() override;

      int selectP(const particleCand& p);

      int mPDGPr = 2212;
      int mPDGNe = 2112;
      int mPDGDe = 4324;
      int mPDGT = 2212 + 2*2112;
      int mPDGHe3 = 2*2212 + 2112;
      int mPDGHe4 = 2*2212 + 2*2112;

};

#endif