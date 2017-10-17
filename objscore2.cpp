#include "muscle.h"
#include "msa.h"
#include "profile.h"
#include "objscore.h"

#define TRACE			0
#define TEST_SPFAST		0

static SCORE ScoreSeqPair(const MSA &msa1, unsigned uSeqIndex1,
  const MSA &msa2, unsigned uSeqIndex2, WEIGHT w, SCORE MatchScore[] = 0)
	{
	const unsigned uColCount = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();
	if (uColCount != uColCount2)
		Quit("ScoreSeqPair, different lengths");

	SCORE scoreTotal = 0;
	bool bGapping1 = false;
	bool bGapping2 = false;
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		bool bGap1 = msa1.IsGap(uSeqIndex1, uColIndex);
		bool bGap2 = msa2.IsGap(uSeqIndex2, uColIndex);

		if (bGap1 && bGap2)
			continue;

		if (bGap1)
			{
			if (!bGapping1)
				{
				scoreTotal += g_scoreGapOpen;
				bGapping1 = true;
				}
			//else
			//	scoreTotal += g_scoreGapExtend;
			continue;
			}

		else if (bGap2)
			{
			if (!bGapping2)
				{
				scoreTotal += g_scoreGapOpen;
				bGapping2 = true;
				}
			//else
			//	scoreTotal += g_scoreGapExtend;
			continue;
			}

		bGapping1 = false;
		bGapping2 = false;

		if (msa1.IsWildcard(uSeqIndex1, uColIndex) ||
		  msa2.IsWildcard(uSeqIndex2, uColIndex))
			continue;

		unsigned uLetter1 = msa1.GetLetter(uSeqIndex1, uColIndex);
		unsigned uLetter2 = msa2.GetLetter(uSeqIndex2, uColIndex);

		SCORE scoreMatch = (*g_ptrScoreMatrix)[uLetter1][uLetter2];
		if (0 != MatchScore)
			MatchScore[uColIndex] += w*scoreMatch;
		scoreTotal += scoreMatch;
		}
	return scoreTotal;
	}

// The usual sum-of-pairs objective score: sum the score
// of the alignment of each pair of sequences.
SCORE ObjScoreSP(const MSA &msa, SCORE MatchScore[])
	{
	if (0 != MatchScore)
		{
		const unsigned uColCount = msa.GetColCount();
		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			MatchScore[uColIndex] = 0;
		}

	const unsigned uSeqCount = msa.GetSeqCount();
	SCORE scoreTotal = 0;
	unsigned uPairCount = 0;
#if	TRACE
	Log("     Score  Weight  Weight       Total\n");
	Log("----------  ------  ------  ----------\n");
#endif
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount; ++uSeqIndex1)
		{
		const WEIGHT w1 = msa.GetSeqWeight(uSeqIndex1);
		for (unsigned uSeqIndex2 = uSeqIndex1 + 1; uSeqIndex2 < uSeqCount; ++uSeqIndex2)
			{
			const WEIGHT w2 = msa.GetSeqWeight(uSeqIndex2);
			const WEIGHT w = w1*w2;
			SCORE scorePair = ScoreSeqPair(msa, uSeqIndex1, msa, uSeqIndex2, w,
			  MatchScore);
			scoreTotal += w1*w2*scorePair;
			++uPairCount;
#if	TRACE
			Log("%10.2f  %6.3f  %6.3f  %10.2f  >%s >%s\n",
			  scorePair,
			  w1,
			  w2,
			  scorePair*w1*w2,
			  msa.GetSeqName(uSeqIndex1),
			  msa.GetSeqName(uSeqIndex2));
#endif
			}
		}
#if	TEST_SPFAST
	{
	SCORE f = ObjScoreSPFast(msa);
	Log("Fast  = %.6g\n", f);
	Log("Brute = %.6g\n", scoreTotal);
	if (BTEq(f, scoreTotal))
		Log("Agree\n");
	else
		Log("** DISAGREE **\n");
	}
#endif
//	return scoreTotal / uPairCount;
	return scoreTotal;
	}

// Objective score defined as the dynamic programming score.
// Input is two alignments, which must be of the same length.
// Result is the same profile-profile score that is optimized
// by dynamic programming.
SCORE ObjScoreDP(const MSA &msa1, const MSA &msa2, SCORE MatchScore[])
	{
	const unsigned uColCount = msa1.GetColCount();
	if (msa2.GetColCount() != uColCount)
		Quit("ObjScoreDP, must be same length");

	const unsigned uColCount1 = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();

	bool bTermGapsFree1 = TermGapsFree(uColCount1, uColCount2);
	bool bTermGapsFree2 = TermGapsFree(uColCount2, uColCount1);

	const ProfPos *PA = ProfileFromMSA(msa1, g_scoreGapOpen, bTermGapsFree1);
	const ProfPos *PB = ProfileFromMSA(msa2, g_scoreGapOpen, bTermGapsFree2);

#if	TRACE
	Log("Profile 1:\n");
	ListProfile(PA, uColCount, &msa1);

	Log("Profile 2:\n");
	ListProfile(PB, uColCount, &msa2);
#endif

	SCORE scoreTotal = 0;

	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		const ProfPos &PPA = PA[uColIndex];
		const ProfPos &PPB = PB[uColIndex];

		SCORE scoreGap = 0;
		SCORE scoreMatch = 0;
	// If gapped column...
		if (PPA.m_bAllGaps && PPB.m_bAllGaps)
			scoreGap = 0;
		else if (PPA.m_bAllGaps)
			{
			if (uColCount - 1 == uColIndex || !PA[uColIndex+1].m_bAllGaps)
				scoreGap = PPB.m_scoreGapClose;
			if (0 == uColIndex || !PA[uColIndex-1].m_bAllGaps)
				scoreGap += PPB.m_scoreGapOpen;
			//if (0 == scoreGap)
			//	scoreGap = PPB.m_scoreGapExtend;
			}
		else if (PPB.m_bAllGaps)
			{
			if (uColCount - 1 == uColIndex || !PB[uColIndex+1].m_bAllGaps)
				scoreGap = PPA.m_scoreGapClose;
			if (0 == uColIndex || !PB[uColIndex-1].m_bAllGaps)
				scoreGap += PPA.m_scoreGapOpen;
			//if (0 == scoreGap)
			//	scoreGap = PPA.m_scoreGapExtend;
			}
		else
			scoreMatch = ScoreProfPos2(PPA, PPB);

		if (0 != MatchScore)
			MatchScore[uColIndex] = scoreMatch;

		scoreTotal += scoreMatch + scoreGap;

#if	TRACE
		{
		const unsigned uSeqCount1 = msa1.GetSeqCount();
		const unsigned uSeqCount2 = msa2.GetSeqCount();

		for (unsigned n = 0; n < uSeqCount1; ++n)
			Log("%c", msa1.GetChar(n, uColIndex));
		Log("  ");
		for (unsigned n = 0; n < uSeqCount2; ++n)
			Log("%c", msa2.GetChar(n, uColIndex));
		Log("  %10.3f", scoreMatch);
		if (scoreGap != 0)
			Log("  %10.3f", scoreGap);
		Log("\n");
		}
#endif
		}

	delete[] PA;
	delete[] PB;

	return scoreTotal;
	}

// Objective score defined as the sum of profile-sequence
// scores for each sequence in the alignment. The profile
// is computed from the entire alignment, so this includes
// the score of each sequence against itself. This is to
// avoid recomputing the profile each time, so we reduce
// complexity but introduce a questionable approximation.
// The goal is to see if we can exploit the apparent
// improvement in performance of log-expectation score
// over the usual sum-of-pairs by optimizing this
// objective score in the iterative refinement stage.
SCORE ObjScorePS(const MSA &msa, SCORE MatchScore[])
	{
	if (g_PPScore != PPSCORE_LE)
		Quit("FastScoreMSA_LASimple: LA");

	const unsigned uSeqCount = msa.GetSeqCount();
	const unsigned uColCount = msa.GetColCount();

	const ProfPos *Prof = ProfileFromMSA(msa, g_scoreGapOpen, false);

	if (0 != MatchScore)
		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			MatchScore[uColIndex] = 0;

	SCORE scoreTotal = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const WEIGHT weightSeq = msa.GetSeqWeight(uSeqIndex);
		SCORE scoreSeq = 0;
		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			const ProfPos &PP = Prof[uColIndex];
			if (msa.IsGap(uSeqIndex, uColIndex))
				{
				bool bOpen = (0 == uColIndex ||
				  !msa.IsGap(uSeqIndex, uColIndex - 1));
				bool bClose = (uColCount - 1 == uColIndex ||
				  !msa.IsGap(uSeqIndex, uColIndex + 1));

				if (bOpen)
					scoreSeq += PP.m_scoreGapOpen;
				if (bClose)
					scoreSeq += PP.m_scoreGapClose;
				//if (!bOpen && !bClose)
				//	scoreSeq += PP.m_scoreGapExtend;
				}
			else if (msa.IsWildcard(uSeqIndex, uColIndex))
				continue;
			else
				{
				unsigned uLetter = msa.GetLetter(uSeqIndex, uColIndex);
				const SCORE scoreMatch = PP.m_AAScores[uLetter];
				if (0 != MatchScore)
					MatchScore[uColIndex] += weightSeq*scoreMatch;
				scoreSeq += scoreMatch;
				}
			}
		scoreTotal += weightSeq*scoreSeq;
		}

	delete[] Prof;
	return scoreTotal;
	}

// The XP score is the sum of the score of each pair of
// sequences between two profiles which are aligned to
// each other. Notice that for two given profiles aligned
// in different ways, the difference in XP score must be
// the same as the difference in SP score because the
// score of a pair of sequences in one profile doesn't
// depend on the alignment.
SCORE ObjScoreXP(const MSA &msa1, const MSA &msa2)
	{
	const unsigned uColCount1 = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();
	if (uColCount1 != uColCount2)
		Quit("ObjScoreXP, alignment lengths differ %u %u", uColCount1, uColCount2);

	const unsigned uSeqCount1 = msa1.GetSeqCount();
	const unsigned uSeqCount2 = msa2.GetSeqCount();

#if	TRACE
	Log("     Score  Weight  Weight       Total\n");
	Log("----------  ------  ------  ----------\n");
#endif

	SCORE scoreTotal = 0;
	unsigned uPairCount = 0;
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount1; ++uSeqIndex1)
		{
		const WEIGHT w1 = msa1.GetSeqWeight(uSeqIndex1);
		for (unsigned uSeqIndex2 = 0; uSeqIndex2 < uSeqCount2; ++uSeqIndex2)
			{
			const WEIGHT w2 = msa2.GetSeqWeight(uSeqIndex2);
			const WEIGHT w = w1*w2;
			SCORE scorePair = ScoreSeqPair(msa1, uSeqIndex1, msa2, uSeqIndex2, w);
			scoreTotal += w1*w2*scorePair;
			++uPairCount;
#if	TRACE
			Log("%10.2f  %6.3f  %6.3f  %10.2f  >%s >%s\n",
			  scorePair,
			  w1,
			  w2,
			  scorePair*w1*w2,
			  msa1.GetSeqName(uSeqIndex1),
			  msa2.GetSeqName(uSeqIndex2));
#endif
			}
		}
	if (0 == uPairCount)
		Quit("0 == uPairCount");

#if	TRACE
	Log("msa1=\n");
	msa1.LogMe();
	Log("msa2=\n");
	msa2.LogMe();
	Log("XP=%g\n", scoreTotal);
#endif
//	return scoreTotal / uPairCount;
	return scoreTotal;
	}
