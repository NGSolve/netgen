//File for handling warnings, errors, messages
#include <meshing.hpp>

namespace netgen
{
  int printwarnings = 1;
  int printerrors = 1;
  int printdots = 1;
  int printfnstart = 0;

  //the dots for progression of program
  void PrintDot(char ch)
  {
    if (printmessage_importance >= 4)
      {
        char st[2] = { ch, 0 };
        Ng_PrintDest(st);
      }
  }
}
