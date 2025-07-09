#ifndef CRUNSPECTRUM
#define CRUNSPECTRUM


#include "Cvrun.h"
#include "TH1D.h"

class runSpectrum : public vrun
{
    public:
      runSpectrum(int nmix=10, TString dirName = "runSpectrum") : vrun(nmix, dirName){}
      void setRapidityRange(float minRapidity, float maxRapidity) {mMinRapidity = minRapidity, mMaxRapidity = maxRapidity;}
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

      TH1D *mHAntiProtonSpectrum = nullptr;
      TH1D *mHAntiNeutronSpectrum = nullptr;
      TH1D *mHAntiDeuteronSpectrum = nullptr;
      TH1D *mHAntiTritiumSpectrum = nullptr;
      TH1D *mHAntiHelium3Spectrum = nullptr;
      TH1D *mHAntiHelium4Spectrum = nullptr;
      TH1D *mHAntiProtonSpectrumA = nullptr;
      TH1D *mHAntiNeutronSpectrumA = nullptr;
      TH1D *mHAntiDeuteronSpectrumA = nullptr;
      TH1D *mHAntiTritiumSpectrumA = nullptr;
      TH1D *mHAntiHelium3SpectrumA = nullptr;
      TH1D *mHAntiHelium4SpectrumA = nullptr;

      TH2D *mHProtonPEta = nullptr;
      TH2D *mHNeutronPEta = nullptr;
      TH2D *mHDeuteronPEta = nullptr;
      TH2D *mHTritiumPEta = nullptr;
      TH2D *mHHelium3PEta = nullptr;
      TH2D *mHHelium4PEta = nullptr;

      TH2D *mHAntiProtonPEta = nullptr;
      TH2D *mHAntiNeutronPEta = nullptr;
      TH2D *mHAntiDeuteronPEta = nullptr;
      TH2D *mHAntiTritiumPEta = nullptr;
      TH2D *mHAntiHelium3PEta = nullptr;
      TH2D *mHAntiHelium4PEta = nullptr;

      TH2D *mHProtonPtEta = nullptr;
      TH2D *mHNeutronPtEta = nullptr;
      TH2D *mHDeuteronPtEta = nullptr;
      TH2D *mHTritiumPtEta = nullptr;
      TH2D *mHHelium3PtEta = nullptr;
      TH2D *mHHelium4PtEta = nullptr;

      TH2D *mHAntiProtonPtEta = nullptr;
      TH2D *mHAntiNeutronPtEta = nullptr;
      TH2D *mHAntiDeuteronPtEta = nullptr;
      TH2D *mHAntiTritiumPtEta = nullptr;
      TH2D *mHAntiHelium3PtEta = nullptr;
      TH2D *mHAntiHelium4PtEta = nullptr;

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

      int mPDGAntiPr = -2212;
      int mPDGAntiNe = -2112;
      int mPDGAntiDe = -4324;
      int mPDGAntiT = -2212 - 2*2112;
      int mPDGAntiHe3 = -2*2212 - 2112;
      int mPDGAntiHe4 = -2*2212 - 2*2112;

      float mMinRapidity = -0.5;
      float mMaxRapidity = 0.5;

};

#endif