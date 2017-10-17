#include "muscle.h"
#include <math.h>
#include "pwpath.h"
#include "profile.h"
#include <stdio.h>

#define	TRACE	0

SCORE ScoreProfPos2LA(const ProfPos &PPA, const ProfPos &PPB)
	{
	SCORE Score = 0;
	for (unsigned n = 0; n < 20; ++n)
		{
		const unsigned uLetter = PPA.m_uSortOrder[n];
		const FCOUNT fcLetter = PPA.m_fcCounts[uLetter];
		if (0 == fcLetter)
			break;
		Score += fcLetter*PPB.m_AAScores[uLetter];
		}
	if (0 == Score)
		return -2.5;
	SCORE logScore = logf(Score);
	return (SCORE) ((logScore - g_scoreCenter)*(PPA.m_fOcc * PPB.m_fOcc));
	}

SCORE ScoreProfPos2NS(const ProfPos &PPA, const ProfPos &PPB)
	{
	SCORE Score = 0;
	for (unsigned n = 0; n < 20; ++n)
		{
		const unsigned uLetter = PPA.m_uSortOrder[n];
		const FCOUNT fcLetter = PPA.m_fcCounts[uLetter];
		if (0 == fcLetter)
			break;
		Score += fcLetter*PPB.m_AAScores[uLetter];
		}
	return Score - g_scoreCenter;
	}

SCORE ScoreProfPos2SP(const ProfPos &PPA, const ProfPos &PPB)
	{
	SCORE Score = 0;
	for (unsigned n = 0; n < 20; ++n)
		{
		const unsigned uLetter = PPA.m_uSortOrder[n];
		const FCOUNT fcLetter = PPA.m_fcCounts[uLetter];
		if (0 == fcLetter)
			break;
		Score += fcLetter*PPB.m_AAScores[uLetter];
		}
	return Score - g_scoreCenter;
	}

static const char *LocalScoreToStr(SCORE s)
	{
	static char str[16];
	if (MINUS_INFINITY == s)
		return "     *";
	sprintf(str, "%6.0f", s);
	return str;
	}

static char ConsensusChar(const ProfPos &PP)
	{
	unsigned uMostCommonLetter = 0;
	FCOUNT fcMostCommon = PP.m_fcCounts[0];
	bool bMoreThanOneLetter = false;
	bool bAnyLetter = false;
	for (unsigned uLetter = 0; uLetter < 20; ++uLetter)
		{
		const FCOUNT fc = PP.m_fcCounts[uLetter];
		if (fc > 0)
			{
			if (bAnyLetter)
				bMoreThanOneLetter = true;
			bAnyLetter = true;
			}
		if (fc > fcMostCommon)
			{
			uMostCommonLetter = uLetter;
			fcMostCommon = fc;
			}
		}
	if (!bAnyLetter)
		return '-';
	char c = LetterToCharAmino(uMostCommonLetter);
	if (bMoreThanOneLetter)
		return tolower(c);
	return c;
	}

static void ListDP(const SCORE *DPM_, const ProfPos *PA, const ProfPos *PB,
  unsigned uPrefixCountA, unsigned uPrefixCountB)
	{
	Log("        ");
	for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		char c = ' ';
		if (uPrefixLengthB > 0)
			c = ConsensusChar(PB[uPrefixLengthB - 1]);
		Log(" %4u:%c", uPrefixLengthB, c);
		}
	Log("\n");
	for (unsigned uPrefixLengthA = 0; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
		char c = ' ';
		if (uPrefixLengthA > 0)
			c = ConsensusChar(PA[uPrefixLengthA - 1]);
		Log("%4u:%c  ", uPrefixLengthA, c);
		for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
			Log(" %s", LocalScoreToStr(DPM(uPrefixLengthA, uPrefixLengthB)));
		Log("\n");
		}
	}

SCORE ScoreProfPos2(const ProfPos &PPA, const ProfPos &PPB)
	{
	if (PPSCORE_SP == g_PPScore)
		return ScoreProfPos2NS(PPA, PPB);
	else if (PPSCORE_LE == g_PPScore)
		return ScoreProfPos2LA(PPA, PPB);
	else if (PPSCORE_SV == g_PPScore)
		return ScoreProfPos2SP(PPA, PPB);
	Quit("Invalid g_PPScore");
	return 0;
	}

SCORE GlobalAlignSimple(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	assert(uLengthB > 0 && uLengthA > 0);

	const unsigned uPrefixCountA = uLengthA + 1;
	const unsigned uPrefixCountB = uLengthB + 1;

// Allocate DP matrices
	const size_t LM = uPrefixCountA*uPrefixCountB;
	SCORE *DPM_ = new SCORE[LM];
	SCORE *DPD_ = new SCORE[LM];
	SCORE *DPI_ = new SCORE[LM];

	DPM(0, 0) = 0;
	DPD(0, 0) = MINUS_INFINITY;
	DPI(0, 0) = MINUS_INFINITY;

	DPM(1, 0) = MINUS_INFINITY;
	DPD(1, 0) = PA[0].m_scoreGapOpen;
	DPI(1, 0) = MINUS_INFINITY;

	DPM(0, 1) = MINUS_INFINITY;
	DPD(0, 1) = MINUS_INFINITY;
	DPI(0, 1) = PB[0].m_scoreGapOpen;

// Empty prefix of B is special case
	for (unsigned uPrefixLengthA = 2; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
	// M=LetterA+LetterB, impossible with empty prefix
		DPM(uPrefixLengthA, 0) = MINUS_INFINITY;

	// D=LetterA+GapB
		DPD(uPrefixLengthA, 0) = DPD(uPrefixLengthA - 1, 0);

	// I=GapA+LetterB, impossible with empty prefix
		DPI(uPrefixLengthA, 0) = MINUS_INFINITY;
		}

// Empty prefix of A is special case
	for (unsigned uPrefixLengthB = 2; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
	// M=LetterA+LetterB, impossible with empty prefix
		DPM(0, uPrefixLengthB) = MINUS_INFINITY;

	// D=LetterA+GapB, impossible with empty prefix
		DPD(0, uPrefixLengthB) = MINUS_INFINITY;

	// I=GapA+LetterB
		DPI(0, uPrefixLengthB) = DPI(0, uPrefixLengthB - 1);
		}

// ============
// Main DP loop
// ============
	SCORE scoreGapCloseB = MINUS_INFINITY;
	for (unsigned uPrefixLengthB = 1; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		const ProfPos &PPB = PB[uPrefixLengthB - 1];

		SCORE scoreGapCloseA = MINUS_INFINITY;
		for (unsigned uPrefixLengthA = 1; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
			{
			const ProfPos &PPA = PA[uPrefixLengthA - 1];

			{
		// Match M=LetterA+LetterB
			SCORE scoreLL = ScoreProfPos2(PPA, PPB);

			SCORE scoreMM = DPM(uPrefixLengthA-1, uPrefixLengthB-1);
			SCORE scoreDM = DPD(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapCloseA;
			SCORE scoreIM = DPI(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapCloseB;

			SCORE scoreBest;
			if (scoreMM >= scoreDM && scoreMM >= scoreIM)
				scoreBest = scoreMM;
			else if (scoreDM >= scoreMM && scoreDM >= scoreIM)
				scoreBest = scoreDM;
			else 
				{
				assert(scoreIM >= scoreMM && scoreIM >= scoreDM);
				scoreBest = scoreIM;
				}
			DPM(uPrefixLengthA, uPrefixLengthB) = scoreBest + scoreLL;
			}

			{
		// Delete D=LetterA+GapB
			SCORE scoreMD = DPM(uPrefixLengthA-1, uPrefixLengthB) +
			  PA[uPrefixLengthA-1].m_scoreGapOpen;
			SCORE scoreDD = DPD(uPrefixLengthA-1, uPrefixLengthB);

			SCORE scoreBest;
			if (scoreMD >= scoreDD)
				scoreBest = scoreMD;
			else
				{
				assert(scoreDD >= scoreMD);
				scoreBest = scoreDD;
				}
			DPD(uPrefixLengthA, uPrefixLengthB) = scoreBest;
			}

		// Insert I=GapA+LetterB
			{
			SCORE scoreMI = DPM(uPrefixLengthA, uPrefixLengthB-1) +
			  PB[uPrefixLengthB - 1].m_scoreGapOpen;
			SCORE scoreII = DPI(uPrefixLengthA, uPrefixLengthB-1);

			SCORE scoreBest;
			if (scoreMI >= scoreII)
				scoreBest = scoreMI;
			else 
				{
				assert(scoreII > scoreMI);
				scoreBest = scoreII;
				}
			DPI(uPrefixLengthA, uPrefixLengthB) = scoreBest;
			}

			scoreGapCloseA = PPA.m_scoreGapClose;
			}
		scoreGapCloseB = PPB.m_scoreGapClose;
		}

#if TRACE
	Log("DPM:\n");
	ListDP(DPM_, PA, PB, uPrefixLengthA, uPrefixLengthB);
	Log("DPD:\n");
	ListDP(DPD_, PA, PB, uPrefixLengthA, uPrefixLengthB);
	Log("DPI:\n");
	ListDP(DPI_, PA, PB, uPrefixLengthA, uPrefixLengthB);
#endif

	SCORE Score = TraceBack(PA, uLengthA, PB, uLengthB, DPM_, DPD_, DPI_, Path);

#if	TRACE
	SCORE scorePath = FastScorePath2(PA, uLengthA, PB, uLengthB, Path);
	Path.LogMe();
	Log("Score = %s Path = %s\n", LocalScoreToStr(Score), LocalScoreToStr(scorePath));
#endif

	delete[] DPM_;
	delete[] DPD_;
	delete[] DPI_;

	return Score;
	}
