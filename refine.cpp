#include "muscle.h"
#include "textfile.h"
#include "seqvect.h"
#include "distfunc.h"
#include "msa.h"
#include "tree.h"
#include "clust.h"
#include "profile.h"
#include "clustsetmsa.h"

void Refine()
	{
	SetOutputFileName(g_pstrInFileName);
	SetInputFileName(g_pstrInFileName);
	SetStartTime();

	SetMaxIters(g_uMaxIters);
	SetSeqWeightMethod(g_SeqWeight1);

	TextFile fileIn(g_pstrInFileName);
	MSA msa;
	msa.FromFile(fileIn);

	const unsigned uSeqCount = msa.GetSeqCount();
	MSA::SetIdCount(uSeqCount);

// Initialize sequence ids.
// From this point on, ids must somehow propogate from here.
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		msa.SetSeqId(uSeqIndex, uSeqIndex);
	SetMuscleInputMSA(msa);

	Tree GuideTree;
	TreeFromMSA(msa, GuideTree, g_Cluster2, g_Distance2, g_Root2);
	SetMuscleTree(GuideTree);

	if (g_bAnchors)
		RefineVert(msa, GuideTree, g_uMaxIters);
	else
		RefineHoriz(msa, GuideTree, g_uMaxIters, false, false);

	ValidateMuscleIds(msa);
	ValidateMuscleIds(GuideTree);

	TextFile fileOut(g_pstrOutFileName, true);
	msa.ToFile(fileOut);
	}
