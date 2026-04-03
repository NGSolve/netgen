//File for handling warnings, errors, messages
#include <meshing.hpp>

namespace netgen
{
  // int printmessage_importance = 3;
int printwarnings = 1;
int printerrors = 1;
int printdots = 1;
int printfnstart = 0;

// extern void Ng_PrintDest(const MyStr& s);
extern void Ng_PrintDest(const char * s);

//the dots for progression of program
void PrintDot(char ch)
{
  // if (printdots)
  if (printmessage_importance >= 4)
    {
      char st[2];
      st[0] = ch;
      st[1] = 0;
      Ng_PrintDest(st);
    }
}

void PrintMessage(int importance, 
		  const MyStr& s1, const MyStr& s2)
{
  if (importance <= printmessage_importance)
    {
      Ng_PrintDest(MyStr(" ")+s1+s2+MyStr("\n"));
    }
}

void PrintMessage(int importance, 
		  const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4)
{
  if (importance <= printmessage_importance)
    {
      Ng_PrintDest(MyStr(" ")+s1+s2+s3+s4+MyStr("\n"));
    }
}

void PrintMessage(int importance, 
		  const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		  const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (importance <= printmessage_importance)
    {
      Ng_PrintDest(MyStr(" ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
    }
}

void PrintMessageCR(int importance, 
		    const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		    const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (importance <= printmessage_importance)
    {
      Ng_PrintDest(MyStr(" ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\r"));
    }
}

void PrintFnStart(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		  const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (printfnstart)
    Ng_PrintDest(MyStr(" Start Function: ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}

void PrintWarning(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		  const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (printwarnings)
    Ng_PrintDest(MyStr(" WARNING: ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}

void PrintError(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (printerrors)
    Ng_PrintDest(MyStr(" ERROR: ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}

void PrintFileError(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		    const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (printerrors)
    Ng_PrintDest(MyStr(" FILE ERROR: ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}

void PrintUserError(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  Ng_PrintDest(MyStr(" USER ERROR: ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}

void PrintSysError(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (printerrors)
    Ng_PrintDest(MyStr(" SYSTEM ERROR: ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}

void PrintTime(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
	       const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (printmessage_importance >= 3)
    Ng_PrintDest(MyStr(" Time = ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}


/*
#ifdef SMALLLIB
#define SMALLLIBORNOTCL
#endif
#ifdef NOTCL
#define SMALLLIBORNOTCL
#endif

#ifdef SMALLLIBORNOTCL
void Ng_PrintDest(const char * s){cout << s <<flush;}
double GetTime(){return 0;}
void MyError(const char * ch)
{
  cerr << ch << endl;
}
#endif
*/


}
