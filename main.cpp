#include "muscle.h"
#include <stdio.h>
#ifdef	WIN32
#include <io.h>	// for isatty()
#else
#include <unistd.h>	// for isatty()
#endif

int g_argc;
char **g_argv;

int main(int argc, char **argv)
	{
	g_argc = argc;
	g_argv = argv;

	SetNewHandler();
	SetStartTime();
	ProcessArgVect(argc - 1, argv + 1);
	SetParams();
	SetLogFile();

	if (g_bVersion)
		{
		printf(MUSCLE_LONG_VERSION "\n");
		exit(EXIT_SUCCESS);
		}

	if (!g_bQuiet)
		Credits();

	if (MissingCommand() && isatty(0))
		{
		Usage();
		exit(EXIT_SUCCESS);
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
