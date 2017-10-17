#include "muscle.h"
#include "objscore.h"
#include "profile.h"
#include "enumopts.h"
#ifdef	WIN32
#include <io.h>	// for isatty()
#else
#include <unistd.h>	// for isatty()
#endif

SCORE g_scoreCenter = 0;
SCORE g_scoreGapExtend = 0;
SCORE g_scoreGapAmbig = 0;
SCORE g_scoreAmbigFactor = 0;

const char *g_pstrInFileName = "-";
const char *g_pstrOutFileName = "-";

const char *g_pstrFileName1 = 0;
const char *g_pstrFileName2 = 0;

const char *g_pstrSPFileName = 0;

const char *g_pstrProf1FileName = 0;
const char *g_pstrProf2FileName = 0;

unsigned g_uSmoothWindowLength = 7;
unsigned g_uAnchorSpacing = 32;
unsigned g_uMaxTreeRefineIters = 1;

unsigned g_uHydrophobicRunLength = 5;
float g_dHydroFactor = (float) 1.2;

unsigned g_uMinDiagLength = 24;
unsigned g_uMaxDiagBreak = 1;
unsigned g_uDiagMargin = 5;

float g_dSUEFF = (float) 0.1;

bool g_bPrecompiledCenter = true;
bool g_bTermGapsHalf = true;
bool g_bTermGapsHalfLonger = false;
bool g_bNormalizeCounts = false;
bool g_bDiags1 = false;
bool g_bDiags2 = false;
bool g_bAnchors = true;
bool g_bQuiet = false;
bool g_bVerbose = false;
bool g_bRefine = false;
bool g_bLow = false;
bool g_bSW = false;
bool g_bTermGaps4 = true;
bool g_bCluster = false;
bool g_bProfile = false;
bool g_bBrenner = false;
bool g_bDimer = false;

#if	DEBUG
bool g_bCatchExceptions = false;
#else
bool g_bCatchExceptions = true;
#endif

bool g_bMSF = false;
bool g_bAln = false;
bool g_bHTML = false;

unsigned g_uMaxIters = 16;
unsigned long g_ulMaxSecs = 0;

PPSCORE g_PPScore = PPSCORE_LE;
OBJSCORE g_ObjScore = OBJSCORE_SPM;

SEQWEIGHT g_SeqWeight1 = SEQWEIGHT_ClustalW;
SEQWEIGHT g_SeqWeight2 = SEQWEIGHT_ClustalW;

DISTANCE g_Distance1 = DISTANCE_Kmer6_6;
DISTANCE g_Distance2 = DISTANCE_PctIdKimura;

CLUSTER g_Cluster1 = CLUSTER_UPGMB;
CLUSTER g_Cluster2 = CLUSTER_UPGMB;

ROOT g_Root1 = ROOT_Pseudo;
ROOT g_Root2 = ROOT_Pseudo;

bool g_bDiags;

//------------------------------------------------------
// These parameters depending on the chosen prof-prof
// score (g_PPScore), initialized to "Undefined".
float g_dSmoothScoreCeil = fInsane;
float g_dMinBestColScore = fInsane;
float g_dMinSmoothScore = fInsane;
SCORE g_scoreGapOpen = fInsane;
//------------------------------------------------------

static unsigned atou(const char *s)
	{
	return (unsigned) atoi(s);
	}

const char *MaxSecsToStr()
	{
	if (0 == g_ulMaxSecs)
		return "(No limit)";
	return SecsToStr(g_ulMaxSecs);
	}

void ListParams()
	{
	Log("\n");
	Log("%s\n", MUSCLE_LONG_VERSION);
	Log("http://www.drive5.com/muscle\n");
	Log("\n");
	Log("Profile-profile score    %s\n", PPSCOREToStr(g_PPScore));
	Log("Max iterations           %u\n", g_uMaxIters);
	Log("Max trees                %u\n", g_uMaxTreeRefineIters);
	Log("Max time                 %s\n", MaxSecsToStr());
	Log("Gap open                 %g\n", g_scoreGapOpen);
	Log("Gap extend (dimer)       %g\n", g_scoreGapExtend);
	Log("Gap ambig factor         %g\n", g_scoreAmbigFactor);
	Log("Gap ambig penalty        %g\n", g_scoreGapAmbig);
	Log("Center (LE)              %g\n", g_scoreCenter);

	Log("Smooth window length     %u\n", g_uSmoothWindowLength);
	Log("Min anchor spacing       %u\n", g_uAnchorSpacing);
	Log("Min diag length (lambda) %u\n", g_uMinDiagLength);
	Log("Diag margin (mu)         %u\n", g_uDiagMargin);
	Log("Min diag break           %u\n", g_uMaxDiagBreak);
	Log("Hydrophobic window       %u\n", g_uHydrophobicRunLength);

	Log("Hydrophobic gap factor   %g\n", g_dHydroFactor);
	Log("Smooth score ceiling     %g\n", g_dSmoothScoreCeil);
	Log("Min best col score       %g\n", g_dMinBestColScore);
	Log("Min anchor score         %g\n", g_dMinSmoothScore);
	Log("SUEFF                    %g\n", g_dSUEFF);

	Log("Term gaps half           %s\n", BoolToStr(g_bTermGapsHalf));
	Log("Term gaps half longer    %s\n", BoolToStr(g_bTermGapsHalfLonger));
	Log("Term gaps four ways      %s\n", BoolToStr(g_bTermGaps4));
	Log("Brenner root MSA         %s\n", BoolToStr(g_bBrenner));
	Log("Normalize counts         %s\n", BoolToStr(g_bNormalizeCounts));
	Log("Diagonals (1)            %s\n", BoolToStr(g_bDiags1));
	Log("Diagonals (2)            %s\n", BoolToStr(g_bDiags2));
	Log("Anchors                  %s\n", BoolToStr(g_bAnchors));
	Log("MSF output format        %s\n", BoolToStr(g_bMSF));
	Log("ClustalW output format   %s\n", BoolToStr(g_bAln));
	Log("Catch exceptions         %s\n", BoolToStr(g_bCatchExceptions));
	Log("Quiet                    %s\n", BoolToStr(g_bQuiet));
	Log("Refine                   %s\n", BoolToStr(g_bRefine));
	Log("Low complexity profiles  %s\n", BoolToStr(g_bLow));

	Log("Objective score          %s\n", OBJSCOREToStr(g_ObjScore));

	Log("Distance method (1)      %s\n", DISTANCEToStr(g_Distance1));
	Log("Clustering method (1)    %s\n", CLUSTERToStr(g_Cluster1));
	Log("Root method (1)          %s\n", ROOTToStr(g_Root1));
	Log("Sequence weighting (1)   %s\n", SEQWEIGHTToStr(g_SeqWeight1));

	Log("Distance method (2)      %s\n", DISTANCEToStr(g_Distance2));
	Log("Clustering method (2)    %s\n", CLUSTERToStr(g_Cluster2));
	Log("Root method (2)          %s\n", ROOTToStr(g_Root2));
	Log("Sequence weighting (2)   %s\n", SEQWEIGHTToStr(g_SeqWeight2));

	Log("\n");
	}

static void SetDefaultsLE()
	{
	g_scoreGapOpen = (SCORE) -3.00;
	g_scoreCenter = (SCORE) -0.55;

	g_bNormalizeCounts = true;

	//g_dSmoothScoreCeil = 5.0;
	//g_dMinBestColScore = 4.0;
	//g_dMinSmoothScore = 2.0;
	g_dSmoothScoreCeil = 3.0;
	g_dMinBestColScore = 2.0;
	g_dMinSmoothScore = 1.0;
	}

void SetDefaultsSP()
	{
	g_scoreGapOpen = -1439;
	g_scoreCenter = 0.0;	// center pre-added into score mx

	g_bNormalizeCounts = false;

	g_dSmoothScoreCeil = 200.0;
	g_dMinBestColScore = 300.0;
	g_dMinSmoothScore = 125.0;
	}

static void SetDefaultsSV()
	{
	g_scoreGapOpen = -300;
	g_scoreCenter = 0.0;	// center pre-added into score mx

	g_bNormalizeCounts = false;

	g_dSmoothScoreCeil = 90.0;
	g_dMinBestColScore = 130.0;
	g_dMinSmoothScore = 40.0;
	}

static void FlagParam(const char *OptName, bool *ptrParam, bool bValueIfFlagSet)
	{
	bool bIsSet = FlagOpt(OptName);
	if (bIsSet)
		*ptrParam = bValueIfFlagSet;
	}

static void StrParam(const char *OptName, const char **ptrptrParam)
	{
	const char *opt = ValueOpt(OptName);
	if (0 != opt)
		*ptrptrParam = opt;
	}

static void FloatParam(const char *OptName, float *ptrParam)
	{
	const char *opt = ValueOpt(OptName);
	if (0 != opt)
		*ptrParam = (float) atof(opt);
	}

static void UintParam(const char *OptName, unsigned *ptrParam)
	{
	const char *opt = ValueOpt(OptName);
	if (0 != opt)
		*ptrParam = atou(opt);
	}

static void EnumParam(const char *OptName, EnumOpt *Opts, int *Param)
	{
	const char *Value = ValueOpt(OptName);
	if (0 == Value)
		return;

	for (;;)
		{
		if (0 == Opts->pstrOpt)
			Quit("Invalid parameter -%s %s", OptName, Value);
		if (0 == stricmp(Value, Opts->pstrOpt))
			{
			*Param = Opts->iValue;
			return;
			}
		++Opts;
		}
	}

static void SetPPScore()
	{
	if (FlagOpt("SP"))
		g_PPScore = PPSCORE_SP;
	else if (FlagOpt("LE"))
		g_PPScore = PPSCORE_LE;
	else if (FlagOpt("SV"))
		g_PPScore = PPSCORE_SV;

	switch (g_PPScore)
		{
	case PPSCORE_SP:
		SetDefaultsSP();
		break;

	case PPSCORE_LE:
		SetDefaultsLE();
		break;

	case PPSCORE_SV:
		SetDefaultsSV();
		break;

	default:
		Quit("Invalid g_PPScore");
		}

	g_scoreGapExtend = (SCORE) (g_scoreCenter/2.0);
	g_scoreGapAmbig = (SCORE) (g_scoreGapOpen*g_scoreAmbigFactor);

	SetScoreMatrix();
	}

void SetPPScore(PPSCORE p)
	{
	g_PPScore = p;
	SetPPScore();
	}

static void SetMaxSecs()
	{
	float fMaxHours = 0.0;
	FloatParam("MaxHours", &fMaxHours);
	if (0.0 == fMaxHours)
		return;
	g_ulMaxSecs = (unsigned long) (fMaxHours*60*60);
	}

static bool CanDoLowComplexity()
	{
	if (g_SeqWeight1 != SEQWEIGHT_ClustalW)
		return false;
	if (1 == g_uMaxIters)
		return true;
	return g_SeqWeight2 == SEQWEIGHT_ClustalW;
	}

static bool MissingCommand()
	{
	if (strcmp(g_pstrInFileName, "-"))
		return false;
	if (0 != g_pstrFileName1)
		return false;
	if (0 != g_pstrSPFileName)
		return false;
	return true;
	}

void SetParams()
	{
	SetPPScore();
	SetMaxSecs();

	StrParam("in", &g_pstrInFileName);
	StrParam("out", &g_pstrOutFileName);

	StrParam("in1", &g_pstrFileName1);
	StrParam("in2", &g_pstrFileName2);

	StrParam("SPScore", &g_pstrSPFileName);

	FlagParam("TermGapsHalf", &g_bTermGapsHalf, true);
	FlagParam("TermGapsHalfLonger", &g_bTermGapsHalfLonger, true);

	FlagParam("TermGapsFull", &g_bTermGapsHalf, false);
	FlagParam("TermGapsFull", &g_bTermGapsHalfLonger, false);

	FlagParam("Core", &g_bCatchExceptions, false);
	FlagParam("NoCore", &g_bCatchExceptions, true);

	FlagParam("Diags1", &g_bDiags1, true);
	FlagParam("Diags2", &g_bDiags2, true);

	FlagParam("Anchors", &g_bAnchors, true);

	FlagParam("Quiet", &g_bQuiet, true);
	FlagParam("Verbose", &g_bVerbose, true);
	FlagParam("Refine", &g_bRefine, true);
	FlagParam("SW", &g_bSW, true);
	FlagParam("Cluster", &g_bCluster, true);
	FlagParam("Profile", &g_bProfile, true);
	FlagParam("TermGaps4", &g_bTermGaps4, true);
	FlagParam("Brenner", &g_bBrenner, true);
	FlagParam("Dimer", &g_bDimer, true);

	FlagParam("MSF", &g_bMSF, true);
	FlagParam("clw", &g_bAln, true);
	FlagParam("HTML", &g_bHTML, true);

	UintParam("MaxIters", &g_uMaxIters);
	UintParam("MaxTrees", &g_uMaxTreeRefineIters);
	UintParam("SmoothWindow", &g_uSmoothWindowLength);
	UintParam("AnchorSpacing", &g_uAnchorSpacing);
	UintParam("DiagLength", &g_uMinDiagLength);
	UintParam("DiagMargin", &g_uDiagMargin);
	UintParam("DiagBreak", &g_uMaxDiagBreak);
	UintParam("Hydro", &g_uHydrophobicRunLength);

	FloatParam("GapOpen", &g_scoreGapOpen);
	FloatParam("GapExtend", &g_scoreGapExtend);
	FloatParam("GapAmbig", &g_scoreAmbigFactor);
	FloatParam("Center", &g_scoreCenter);
	FloatParam("SmoothScoreCeil", &g_dSmoothScoreCeil);
	FloatParam("MinBestColScore", &g_dMinBestColScore);
	FloatParam("MinSmoothScore", &g_dMinSmoothScore);
	FloatParam("SUEFF", &g_dSUEFF);
	FloatParam("HydroFactor", &g_dHydroFactor);

	EnumParam("ObjScore", OBJSCORE_Opts, (int *) &g_ObjScore);

	EnumParam("Distance1", DISTANCE_Opts, (int *) &g_Distance1);
	EnumParam("Distance2", DISTANCE_Opts, (int *) &g_Distance2);

	EnumParam("Weight1", SEQWEIGHT_Opts, (int *) &g_SeqWeight1);
	EnumParam("Weight2", SEQWEIGHT_Opts, (int *) &g_SeqWeight2);

	EnumParam("Cluster1", CLUSTER_Opts, (int *) &g_Cluster1);
	EnumParam("Cluster2", CLUSTER_Opts, (int *) &g_Cluster2);

	EnumParam("Root1", ROOT_Opts, (int *) &g_Root1);
	EnumParam("Root2", ROOT_Opts, (int *) &g_Root2);

	g_scoreGapAmbig = g_scoreGapOpen*g_scoreAmbigFactor;
	g_bLow = CanDoLowComplexity();

	if (g_bDimer)
		g_bPrecompiledCenter = false;

	if (MissingCommand() && isatty(0))
		{
		Usage();
		exit(0);
		}
	}
