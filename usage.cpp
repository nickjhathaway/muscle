#include "muscle.h"
#include <stdio.h>

void Usage()
	{
	fprintf(stderr,
"Basic usage\n"
"\n"
"    muscle -in <inputfile> -out <outputfile>\n"
"\n"
"Common options\n"
"\n"
"    -in <inputfile>    Input file in FASTA format (default stdin)\n"
"    -out <outputfile>  Output alignment in FASTA format (default stdout)\n"
"    -maxiters <n>      Maximum number of iterations (integer, default 16)\n"
"    -maxhours <h>      Maximum time to iterate in hours (default no limit)\n"
"    -clw               Write output in CLUSTALW format (default FASTA)\n"
"    -html              Write output in HTML format (default FASTA)\n"
"    -msf               Write output in MSF format (default FASTA)\n"
"    -log[a] <logfile>  Log to file (append if -loga, overwrite if -log)\n"
"    -quiet             Do not write progress messages to stderr\n"
"    -verbose           Write parameter information to log file\n"
"    -nocore            Trap exceptions and exit gracefully\n"
"    -core              Allow default exception handling, e.g. core dump\n"
"\n"
"Progressive only (very fast, avg accuracy ~= T-Coffee): -maxiters 2\n"
"Fastest options: -sv -maxiters 1 -diags1 -distance1 kbit20_3\n"
"\n"
"MUSCLE home page http://www.drive5.com/muscle\n");
	}
