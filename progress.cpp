#include "muscle.h"
#include <stdio.h>
#include <time.h>

// Functions that provide visible feedback to the user
// that progress is being made.

static unsigned g_uIter = 0;		// Main MUSCLE iteration 1, 2..
static unsigned g_uLocalMaxIters = 0;	// Max iters
static FILE *g_fProgress = stderr;	// Default to standard error
static char g_strFileName[32];		// File name
static time_t g_tLocalStart;				// Start time
static char g_strDesc[32];			// Description
static bool g_bWipeDesc = false;
static int g_nPrevDescLength;
static unsigned g_uTotalSteps;

const char *ElapsedTimeAsStr()
	{
	time_t Now = time(0);
	unsigned long ElapsedSecs = (unsigned long) (Now - g_tLocalStart);
	return SecsToStr(ElapsedSecs);
	}

const char *MemUseAsStr()
	{
	static char Str[9];
	static double MaxMB = 0;

	double MB = GetMemUseMB();
	if (MB > MaxMB)
		MaxMB = MB;
	sprintf(Str, "%.0f", MaxMB);
	return Str;
	}

void SetInputFileName(const char *pstrFileName)
	{
	NameFromPath(pstrFileName, g_strFileName, sizeof(g_strFileName));
	}

void SetSeqStats(unsigned uSeqCount, unsigned uMaxL, unsigned uAvgL)
	{
	if (g_bQuiet)
		return;

	fprintf(g_fProgress, "%s %u seqs, max length %u, avg  length %u\n",
	  g_strFileName, uSeqCount, uMaxL, uAvgL);
	if (g_bVerbose)
		Log("%u seqs, max length %u, avg  length %u\n",
		  uSeqCount, uMaxL, uAvgL);
	}

void SetStartTime()
	{
	time(&g_tLocalStart);
	}

unsigned long GetStartTime()
	{
	return (unsigned long) g_tLocalStart;
	}

void SetIter(unsigned uIter)
	{
	g_uIter = uIter;
	}

void IncIter()
	{
	++g_uIter;
	}

void SetMaxIters(unsigned uMaxIters)
	{
	g_uLocalMaxIters = uMaxIters;
	}

void SetProgressDesc(const char szDesc[])
	{
	strncpy(g_strDesc, szDesc, sizeof(g_strDesc));
	g_strDesc[sizeof(g_strDesc) - 1] = 0;
	}

static void Wipe(int n)
	{
	for (int i = 0; i < n; ++i)
		fprintf(g_fProgress, " ");
	}

void Progress(unsigned uStep, unsigned uTotalSteps)
	{
	CheckMaxTime();

	if (g_bQuiet)
		return;

	double dPct = ((uStep + 1)*100.0)/uTotalSteps;
	fprintf(g_fProgress, "%8.8s  %4.4s Mb  Iter %3u  %6.2f%%  %s",
	  ElapsedTimeAsStr(),
	  MemUseAsStr(),
	  g_uIter,
	  dPct,
	  g_strDesc);

	if (g_bWipeDesc)
		{
		int n = g_nPrevDescLength - (int) strlen(g_strDesc);
		Wipe(n);
		g_bWipeDesc = false;
		}

	fprintf(g_fProgress, "\r");

	g_uTotalSteps = uTotalSteps;
	}

void ProgressStepsDone()
	{
	CheckMaxTime();

	if (g_bVerbose)
		Log("Elaspsed time %8.8s  Peak memory use %4.4s Mb  Iteration %3u %s\n",
		 ElapsedTimeAsStr(),
		 MemUseAsStr(),
		 g_uIter,
		 g_strDesc);

	if (g_bQuiet)
		return;

	Progress(g_uTotalSteps - 1, g_uTotalSteps);
	fprintf(g_fProgress, "\n");
	g_bWipeDesc = true;
	g_nPrevDescLength = (int) strlen(g_strDesc);
	}
