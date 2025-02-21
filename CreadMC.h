#ifndef CREADMC
#define CREADMC

#include "Cvread.h"

class readerMC : public vreader
{
  public:
    readerMC() : vreader() {}
    readerMC(const char* fname) : vreader(fname) {}
    virtual bool LoadEvent(int iev);
    virtual const char* getTreeName() const { return "tree"; }
    virtual void clear() { vreader::clear(); }

  private:
    virtual void initTree();

};

#endif
