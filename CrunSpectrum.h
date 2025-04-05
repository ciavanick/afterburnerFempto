#ifndef CRUNSPECTRUM
#define CRUNSPECTRUM


#include "Cvrun.h"
#include "TH1D.h"

class runSpectrum : public vrun
{
    public:
      runSpectrum(int nmix=10) : vrun(nmix){}
      void initSpectrum();
      void process() override;
      void write() override;
      void setPDG(int pdg);
    private:
      TH1D *mHSpectrum = nullptr;
      TH1D *mNEvents = nullptr;

      int mPDG = 2212;
};

#endif