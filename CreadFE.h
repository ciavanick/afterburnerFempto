#ifndef CREADFE
#define CREADFE

#include "Cvread.h"

class readerFE : public vreader
{
  public:
    readerFE() : vreader() {}
    readerFE(const char* fname) : vreader(fname) {}
    virtual bool LoadEvent(int iev);
    virtual const char* getTreeName() const { return "h76"; }
    virtual void clear() { vreader::clear(); }

  private:
    virtual void initTree();
    // variables in the tree
    int mNtrack;
    int mPID[10000];              // pdg
//    int mPID_mother1[10000];
//    int mPID_mother2[10000];
//    double mTau0[10000];           // ctau in fm
//    double mTau0_mother1[10000];
    double mPx[10000],mPy[10000],mPz[10000];

};

#endif
