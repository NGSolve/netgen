#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#define maxlen 1000

int main (int argc, char ** argv)
{
  if (argc != 4)
    {
      cout << "use:  makerlsfile  infile outfile rulename" << endl;
      exit(1);
    }
  

  char line[maxlen], infile[maxlen], outfile[maxlen];\
  char ch;
  int i, j;

  /*
  cout << "infile: ";
  cin >> infile;
  cout << "outfile: ";
  cin >> outfile;
  */

  ifstream inf (argv[1]);
  ofstream outf (argv[2]);
  string rulename = argv[3];

  outf << "namespace netgen" << endl << '{' << endl;
  outf << "const char * " << rulename << "[] = {" << endl;
  while (inf.good())
    {
      i = 0;

      inf.get(ch);
      while (ch != '\n' && inf.good() && i < maxlen)
	{
	  if (ch == '\"')
	    {
	      line[i] = '\\';
	      line[i+1] = '\"';
	      i+= 2;
	    }
	  else if (ch == '\\')
	    {
	      line[i] = '\\';
	      line[i+1] = '\\';
	      i+= 2;
	    }
	  else
	    {
	      line[i] = ch;
	      i++;
	    }
	  inf.get(ch);
	}
      line[i] = 0;
      // cout << line << endl;
      outf << "\"" << line << "\\n\",\\" << endl;
    }
  outf << "0};" << endl;
  outf << '}' << endl;
  return 0;
}
