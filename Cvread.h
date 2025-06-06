#ifndef CVREAD
#define CVREAD

#include "Cutils.h"
#include "TString.h"
#include "TChain.h"

class vreader
{
  public:
    vreader() {}
    vreader(const char* fname) : mFileName(fname) {}
    void setFileName(const char* fname) { mFileName = fname; }
    virtual bool LoadEvent(int iev);
    virtual int getNevents() const;
    void NextEvent() { mIEvent++; LoadEvent(mIEvent); };
    int getCurrentEventID() const { return mIEvent; }
    virtual const char* getTreeName() const = 0;
    virtual void openFile(const char *fname = nullptr, bool fromCollection=false);
    virtual void clear();
    void setEtaRange(float minE, float maxE) { mMinEta = minE, mMaxEta = maxE; }
    float getMinEta() const { return mMinEta; }
    float getMaxEta() const { return mMaxEta; }
    std::vector<particleMC>& getParticles() { return mParticles; }

  protected:
    virtual void initTree() = 0;

    TChain *mTree = nullptr;
    int mIEvent=-1;
    bool mIsOpened = false;
    float mMinEta = -3;
    float mMaxEta = 3;
    std::vector<particleMC> mParticles;

  private:
    TString mFileName;
};

#endif
