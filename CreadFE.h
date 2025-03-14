#ifndef CREADFE
#define CREADFE

#include "Cvread.h"

#include "TClonesArray.h"
#include "TParticle.h"

class readerFE : public vreader
{
  public:
    readerFE() : vreader() { mArr = new TClonesArray("TParticle"); }
    readerFE(const char* fname) : vreader(fname) { mArr = new TClonesArray("TParticle"); }
    virtual bool LoadEvent(int iev);
    virtual const char* getTreeName() const { return "T"; }
    virtual void clear() { vreader::clear(); }

  private:
    virtual void initTree();
    // variables in the tree
    TClonesArray *mArr = nullptr;
};

#endif
