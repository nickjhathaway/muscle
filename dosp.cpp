#include "muscle.h"
#include "textfile.h"
#include "msa.h"
#include "objscore.h"
#include "tree.h"
#include "profile.h"

void DoSP()
	{
	TextFile f(g_pstrSPFileName);

	if (PPSCORE_LE == g_PPScore)
		{
		g_PPScore = PPSCORE_SV;
		SetScoreMatrix();
		}

	MSA a;
	a.FromFile(f);

	const unsigned uSeqCount = a.GetSeqCount();
	if (0 == uSeqCount)
		Quit("No sequences in input file %s", g_pstrSPFileName);

	MSA::SetIdCount(uSeqCount);
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		a.SetSeqId(uSeqIndex, uSeqIndex);

	SetSeqWeightMethod(g_SeqWeight1);
	Tree tree;
	TreeFromMSA(a, tree, g_Cluster2, g_Distance2, g_Root2);
	SetMuscleTree(tree);
	SetMSAWeightsMuscle((MSA &) a);

	SCORE SP = ObjScoreSP(a);

	Log("File=%s;SP=%.4g\n", g_pstrSPFileName, SP);
	fprintf(stderr, "File=%s;SP=%.4g\n", g_pstrSPFileName, SP);
	}
