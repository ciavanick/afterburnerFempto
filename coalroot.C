#include <stdio.h>
#include "TApplication.h"
#include "TRint.h"
#include "version/version.h"

int main(int argc, char **argv){
  TRint *theApp = new TRint("coalroot", &argc, argv);

  theApp->Run();

  return 0;
}
