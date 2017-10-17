#include "muscle.h"
#include "msa.h"
#include "profile.h"

#define TRACE	0
#define PAF		0

bool TermGapsFree(unsigned uLength1, unsigned uLength2)
	{
	if (g_bTermGapsHalf)
		return true;
	if (!g_bTermGapsHalfLonger)
		return false;
	return uLength1 > uLength2;
	}

void SortCounts(const FCOUNT fcCounts[], unsigned SortOrder[])
	{
	static unsigned InitialSortOrder[MAX_ALPHA] =
		{
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
		};
	memcpy(SortOrder, InitialSortOrder, sizeof(InitialSortOrder));

	bool bAny = true;
	while (bAny)
		{
		bAny = false;
		for (unsigned n = 0; n < MAX_ALPHA - 1; ++n)
			{
			unsigned i1 = SortOrder[n];
			unsigned i2 = SortOrder[n+1];
			if (fcCounts[i1] < fcCounts[i2])
				{
				SortOrder[n+1] = i1;
				SortOrder[n] = i2;
				bAny = true;
				}
			}
		}
	}

unsigned AminoGroupFromFCounts(const FCOUNT fcCounts[])
	{
	bool bAny = false;
	unsigned uConsensusAminoGroup = AAGROUP_MULTIPLE;
	for (unsigned uLetter = 0; uLetter < MAX_ALPHA; ++uLetter)
		{
		if (0 == fcCounts[uLetter])
			continue;
		const unsigned uAminoGroup = AminoGroup[uLetter];
		if (bAny)
			{
			if (uAminoGroup != uConsensusAminoGroup)
				return AAGROUP_MULTIPLE;
			}
		else
			{
			bAny = true;
			uConsensusAminoGroup = uAminoGroup;
			}
		}
	return uConsensusAminoGroup;
	}

ProfPos *ProfileFromMSA(const MSA &a, SCORE scoreGapOpen, bool bTermGapsFree)
	{
	const unsigned uSeqCount = a.GetSeqCount();
	const unsigned uColCount = a.GetColCount();

// Yuck -- cast away const (inconsistent design here).
	SetMSAWeightsMuscle((MSA &) a);

	ProfPos *Pos = new ProfPos[uColCount];

	unsigned uHydrophobicRunLength = 0;
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		ProfPos &PP = Pos[uColIndex];

		PP.m_bAllGaps = a.IsGapColumn(uColIndex);

		FCOUNT fcGapStart;
		FCOUNT fcGapEnd;
		FCOUNT fcGapExtend;
		FCOUNT fOcc;
		a.GetFractionalWeightedCounts(uColIndex, g_bNormalizeCounts, PP.m_fcCounts,
		  &fcGapStart, &fcGapEnd, &fcGapExtend, &fOcc,
		  &PP.m_LL, &PP.m_LG, &PP.m_GL, &PP.m_GG);
		PP.m_fOcc = fOcc;

		SortCounts(PP.m_fcCounts, PP.m_uSortOrder);

		PP.m_uAminoGroup = AminoGroupFromFCounts(PP.m_fcCounts);

		for (unsigned i = 0; i < 20; ++i)
			{
			SCORE scoreSum = 0;
			for (unsigned j = 0; j < 20; ++j)
				scoreSum += PP.m_fcCounts[j]*(*g_ptrScoreMatrix)[i][j];
			PP.m_AAScores[i] = scoreSum;
			}

		SCORE sStartOcc = (SCORE) (1.0 - fcGapStart);
		SCORE sEndOcc = (SCORE) (1.0 - fcGapEnd);

		PP.m_fcStartOcc = sStartOcc;
		PP.m_fcEndOcc = sEndOcc;

		PP.m_scoreGapOpen = sStartOcc*scoreGapOpen/2;
		PP.m_scoreGapClose = sEndOcc*scoreGapOpen/2;
//		PP.m_scoreGapExtend = (SCORE) ((1.0 - fcGapExtend)*scoreGapExtend);

#if	PAF
		if (sStartOcc > 0.5)
			{
			extern SCORE PAFactor(const FCOUNT fcCounts[]);
			SCORE paf = PAFactor(PP.m_fcCounts);
			PP.m_scoreGapOpen *= paf;
			PP.m_scoreGapClose *= paf;
			}
#endif

		if (bTermGapsFree)
			{
			if (0 == uColIndex)
				PP.m_scoreGapOpen = 0;
			if (uColCount - 1 == uColIndex)
				PP.m_scoreGapClose = 0;
			}
		}

	Hydro(Pos, uColCount);

#if	TRACE
	{
	Log("ProfileFromMSA\n");
	ListProfile(Pos, uColCount, &a);
	}
#endif
	return Pos;
	}

static void LogF(FCOUNT f)
	{
	if (f > -0.00001 && f < 0.00001)
		Log("       ");
	else
		Log("  %5.3f", f);
	}

void ListProfile(const ProfPos *Prof, unsigned uLength, const MSA *ptrMSA)
	{
	Log("  Pos  Occ     LL     LG     GL     GG     Open  Close\n");
	Log("  ---  ---     --     --     --     --     ----  -----\n");
	for (unsigned n = 0; n < uLength; ++n)
		{
		const ProfPos &PP = Prof[n];
		Log("%5u", n);
		LogF(PP.m_fOcc);
		LogF(PP.m_LL);
		LogF(PP.m_LG);
		LogF(PP.m_GL);
		LogF(PP.m_GG);
		Log("  %5.1f", -PP.m_scoreGapOpen);
		Log("  %5.1f", -PP.m_scoreGapClose);
		if (0 != ptrMSA)
			{
			const unsigned uSeqCount = ptrMSA->GetSeqCount();
			Log("  ");
			for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
				Log("%c", ptrMSA->GetChar(uSeqIndex, n));
			}
		Log("\n");
		}

	Log("\n");
	Log("  Pos G");
	for (unsigned n = 0; n < 20; ++n)
		Log("     %c", LetterToCharAmino(n));
	Log("\n");
	Log("  --- -");
	for (unsigned n = 0; n < 20; ++n)
		Log(" -----");
	Log("\n");

	for (unsigned n = 0; n < uLength; ++n)
		{
		const ProfPos &PP = Prof[n];
		Log("%5u", n);
		if (-1 == PP.m_uAminoGroup)
			Log(" -", PP.m_uAminoGroup);
		else
			Log(" %d", PP.m_uAminoGroup);

		for (unsigned uLetter = 0; uLetter < 20; ++uLetter)
			{
			FCOUNT f = PP.m_fcCounts[uLetter];
			if (f == 0.0)
				Log("      ");
			else
				Log(" %5.3f", f);
			}
		if (0 != ptrMSA)
			{
			const unsigned uSeqCount = ptrMSA->GetSeqCount();
			Log("  ");
			for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
				Log("%c", ptrMSA->GetChar(uSeqIndex, n));
			}
		Log("\n");
		}
	}
