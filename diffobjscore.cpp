#include "muscle.h"
#include "msa.h"
#include "objscore.h"
#include "profile.h"

SCORE DiffObjScore(
  const MSA &msa1, const PWPath &Path1, const unsigned Edges1[], unsigned uEdgeCount1, 
  const MSA &msa2, const PWPath &Path2, const unsigned Edges2[], unsigned uEdgeCount2)
	{
	SCORE score1 = ObjScoreEdges(msa1, Path1, Edges1, uEdgeCount1);
	SCORE score2 = ObjScoreEdges(msa1, Path1, Edges1, uEdgeCount1);
	return score1 - score2;
	}

SCORE ObjScoreEdges(const MSA &msa, const PWPath &path, const unsigned Edges[],
  unsigned uEdgeCount)
	{
	if (g_ObjScore != OBJSCORE_SP)
		Quit("ObjScoreEdges: requires SP");

	const unsigned uSeqCount = msa.GetSeqCount();
	const unsigned uColCount = msa.GetColCount();

// Letters
	SCORE s = 0;
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const unsigned uColIndex = Edges[uEdgeIndex];
		assert(uColIndex < uColCount);
		for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount; ++uSeqIndex1)
			{
			WEIGHT w1 = msa.GetSeqWeight(uSeqIndex1);
			unsigned uLetter1 = msa.GetLetterEx(uSeqIndex1, uColIndex);
			if (uLetter1 >= 20)
				continue;
			SCORE *Row = (*g_ptrScoreMatrix)[uLetter1];
			for (unsigned uSeqIndex2 = 0; uSeqIndex2 < uSeqIndex1; ++uSeqIndex2)
				{
				WEIGHT w2 = msa.GetSeqWeight(uSeqIndex2);
				unsigned uLetter2 = msa.GetLetterEx(uSeqIndex2, uColIndex);
				if (uLetter2 >= 20)
					continue;
				s += w1*w2*Row[uLetter2];
				}
			}
		}

// Gaps
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const unsigned uColIndex = Edges[uEdgeIndex];
		assert(uColIndex < uColCount);
		for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount; ++uSeqIndex1)
			{
			WEIGHT w1 = msa.GetSeqWeight(uSeqIndex1);
			bool bGap1 = msa.IsGap(uSeqIndex1, uColIndex);
			bool bGapPrev1 =
			  (uColIndex == 0 ? false : msa.IsGap(uSeqIndex1, uColIndex - 1));
			for (unsigned uSeqIndex2 = 0; uSeqIndex2 < uSeqIndex1; ++uSeqIndex2)
				{
				WEIGHT w2 = msa.GetSeqWeight(uSeqIndex2);
				bool bGapPrev1 =
				  (uColIndex == 0 ? false : msa.IsGap(uSeqIndex1, uColIndex - 1));
				}
			}
		}
	return s;
	}
