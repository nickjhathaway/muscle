#include "muscle.h"
#include <stdio.h>

int g_argc;
char **g_argv;

int main(int argc, char **argv)
	{
	//extern void SPTest();
	//SPTest();
	//exit(0);

	g_argc = argc;
	g_argv = argv;

	SetNewHandler();
	SetStartTime();
	ProcessArgVect(argc - 1, argv + 1);
	SetParams();
	SetLogFile();

	if (!g_bQuiet)
		{
		fprintf(stderr, "\n" MUSCLE_LONG_VERSION "\n\n");
		fprintf(stderr, "http://www.drive5.com/muscle\n");
		fprintf(stderr, "This software is donated to the public domain.\n");
		fprintf(stderr, "Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.\n\n");
		}

	if (g_bCatchExceptions)
		{
		try
			{
			Run();
			}
		catch (...)
			{
			OnException();
			exit(EXIT_Except);
			}
		}
	else
		Run();

	exit(EXIT_Success);
	}
