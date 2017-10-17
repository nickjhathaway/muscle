#include "muscle.h"
#include "textfile.h"
#include "seqvect.h"
#include "distfunc.h"
#include "msa.h"
#include "tree.h"
#include "profile.h"

extern void MHackStart(SeqVect &v);
extern void MHackEnd(MSA &msa);

static void Muscle()
	{
	SetOutputFileName(g_pstrInFileName);
	SetInputFileName(g_pstrInFileName);
	SetStartTime();

	SetMaxIters(g_uMaxIters);
	SetSeqWeightMethod(g_SeqWeight1);

	TextFile fileIn(g_pstrInFileName);
	SeqVect v;
	v.FromFASTAFile(fileIn);
	const unsigned uSeqCount = v.Length();

	if (0 == uSeqCount)
		Quit("No sequences in input file");

	ALPHA Alpha = ALPHA_Undefined;
	switch (g_SeqType)
		{
	case SEQTYPE_Auto:
		Alpha = v.GuessAlpha();
		break;

	case SEQTYPE_Protein:
		Alpha = ALPHA_Amino;
		break;

	case SEQTYPE_Nucleo:
		Alpha = ALPHA_Nucleo;
		break;

	default:
		Quit("Invalid SeqType");
		}
	SetAlpha(Alpha);
	v.FixAlpha();

	if (ALPHA_Nucleo == Alpha)
		{
		SetPPScore(PPSCORE_SPN);
		g_Distance1 = DISTANCE_Kmer4_6;
		}

	if (g_bVerbose)
		ListParams();

	unsigned uMaxL = 0;
	unsigned uTotL = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned L = v.GetSeq(uSeqIndex).Length();
		uTotL += L;
		if (L > uMaxL)
			uMaxL = L;
		}

	SetIter(1);
	g_bDiags = g_bDiags1;
	SetSeqStats(uSeqCount, uMaxL, uTotL/uSeqCount);

	SetMuscleSeqVect(v);

	MSA::SetIdCount(uSeqCount);

// Initialize sequence ids.
// From this point on, ids must somehow propogate from here.
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		v.SetSeqId(uSeqIndex, uSeqIndex);

	if (uSeqCount > 1)
		MHackStart(v);

	if (0 == uSeqCount)
		Quit("Input file '%s' has no sequences", g_pstrInFileName);
	if (1 == uSeqCount)
		{
		TextFile fileOut(g_pstrOutFileName, true);
		v.ToFile(fileOut);
		return;
		}
/////
//	{
//#include "pwpath.h"
//	const Seq &seq1 = v.GetSeq(0);
//	const Seq &seq2 = v.GetSeq(1);
//	PWPath Path;
//	SCORE s = GlobalAlignSS(seq1, seq2, Path);
//	Log("SS = %.3g\n", s);
//	Log("Path=\n");
//	Path.LogMe();
//	}
/////

// First iteration
	Tree GuideTree;
	TreeFromSeqVect(v, GuideTree, g_Cluster1, g_Distance1, g_Root1);

	const char *Tree1 = ValueOpt("Tree1");
	if (0 != Tree1)
		{
		TextFile f(Tree1, true);
		GuideTree.ToFile(f);
		if (g_bCluster)
			return;
		}

	SetMuscleTree(GuideTree);
	ValidateMuscleIds(GuideTree);

	MSA msa;
	ProgNode *ProgNodes = 0;
	if (g_bLow)
		ProgNodes = ProgressiveAlignE(v, GuideTree, msa);
	else
		ProgressiveAlign(v, GuideTree, msa);
	SetCurrentAlignment(msa);

	ValidateMuscleIds(msa);

	if (1 == g_uMaxIters || 2 == uSeqCount)
		{
		TextFile fileOut(g_pstrOutFileName, true);
		MHackEnd(msa);
		msa.ToFile(fileOut);
		return;
		}

	g_bDiags = g_bDiags2;
	SetIter(2);

	if (g_bLow)
		{
		if (0 != g_uMaxTreeRefineIters)
			RefineTreeE(msa, v, GuideTree, ProgNodes);
		}
	else
		RefineTree(msa, GuideTree);

	const char *Tree2 = ValueOpt("Tree2");
	if (0 != Tree2)
		{
		TextFile f(Tree2, true);
		GuideTree.ToFile(f);
		}

	SetSeqWeightMethod(g_SeqWeight2);
	SetMuscleTree(GuideTree);

	if (g_bAnchors)
		RefineVert(msa, GuideTree, g_uMaxIters - 2);
	else
		RefineHoriz(msa, GuideTree, g_uMaxIters - 2, false, false);


#if	0
// Refining by subfamilies is disabled as it didn't give better
// results. I tried doing this before and after RefineHoriz.
// Should get back to this as it seems like this should work.
	RefineSubfams(msa, GuideTree, g_uMaxIters - 2);
#endif

	ValidateMuscleIds(msa);
	ValidateMuscleIds(GuideTree);

	TextFile fileOut(g_pstrOutFileName, true);
	MHackEnd(msa);
	msa.ToFile(fileOut);
	}

void Run()
	{
	Log("Started %s\n", GetTimeAsStr());
	for (int i = 0; i < g_argc; ++i)
		Log("%s ", g_argv[i]);
	Log("\n");

	if (g_bRefine)
		Refine();
	else if (g_bSW)
		Local();
	else if (0 != g_pstrSPFileName)
		DoSP();
	else if (g_bProfile)
		Profile();
	else
		Muscle();

	ListDiagSavings();
	Log("Finished %s\n", GetTimeAsStr());
	}
