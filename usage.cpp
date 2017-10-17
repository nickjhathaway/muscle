#include "muscle.h"
#include <stdio.h>

void Credits()
	{
	static bool Displayed = false;
	if (Displayed)
		return;

	fprintf(stderr, "\n" MUSCLE_LONG_VERSION "\n\n");
	fprintf(stderr, "http://www.drive5.com/muscle\n");
	fprintf(stderr, "This software is donated to the public domain.\n");
	fprintf(stderr, "Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.\n\n");
	Displayed = true;
	}

void Usage()
	{
	Credits();
	fprintf(stderr,
"\n"
"Basic usage\n"
"\n"
"    muscle -in <inputfile> -out <outputfile>\n"
"\n"
"Most important options (for a complete list please see the User Guide):\n"
"\n"
"    -in <inputfile>    Input file in FASTA format (default stdin)\n"
"    -out <outputfile>  Output alignment in FASTA format (default stdout)\n"
"    -maxiters <n>      Maximum number of iterations (integer, default 16)\n"
"    -maxhours <h>      Maximum time to iterate in hours (default no limit)\n"
"    -clw               Write output in CLUSTALW format (default FASTA)\n"
"    -clwstrict         As -clw, with 'CLUSTAL W (1.81)' header\n"
"    -html              Write output in HTML format (default FASTA)\n"
"    -msf               Write output in MSF format (default FASTA)\n"
"    -log[a] <logfile>  Log to file (append if -loga, overwrite if -log)\n"
"    -quiet             Do not write progress messages to stderr\n"
"    -stable            Output sequences in input order (default is -group)\n"
"    -group             Group sequences by similarity (this is the default)\n"
"    -version           Display version information and exit\n"
"\n"
"Progressive only (very fast, avg accuracy ~= T-Coffee): -maxiters 2\n"
"Fastest options: -sv -maxiters 1 -diags1 -distance1 kbit20_3\n");
	}
