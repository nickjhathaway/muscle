#include "muscle.h"
#include "msa.h"
#include "profile.h"
#include "pwpath.h"
#include "textfile.h"
#include "timing.h"

SCORE GlobalAlign4(ProfPos *PA, unsigned uLengthA, ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
static SCORE GlobalAlignXX(int NPenalty, int CPenalty, ProfPos *PA, unsigned uLengthA,
  ProfPos *PB, unsigned uLengthB, PWPath &Path, SCORE *ptrMatch);
static SCORE MatchScore(ProfPos *PA, unsigned uLengthA, ProfPos *PB,
  unsigned uLengthB, PWPath &Path);

SCORE AlignTwoMSAs(const MSA &msa1, const MSA &msa2, MSA &msaOut, PWPath &Path,
  bool bLockLeft, bool bLockRight)
	{
	const unsigned uLengthA = msa1.GetColCount();
	const unsigned uLengthB = msa2.GetColCount();

	SCORE scoreGapOpen = (SCORE) g_scoreGapOpen;

	bool bTermGapsFreeA = TermGapsFree(uLengthA, uLengthB);
	bool bTermGapsFreeB = TermGapsFree(uLengthB, uLengthA);

	ProfPos *PA = ProfileFromMSA(msa1, scoreGapOpen, bTermGapsFreeA);
	ProfPos *PB = ProfileFromMSA(msa2, scoreGapOpen, bTermGapsFreeB);

	if (bLockLeft)
		{
		PA[0].m_scoreGapOpen = MINUS_INFINITY;
		PB[0].m_scoreGapOpen = MINUS_INFINITY;
		}

	if (bLockRight)
		{
		PA[uLengthA-1].m_scoreGapClose = MINUS_INFINITY;
		PB[uLengthB-1].m_scoreGapClose = MINUS_INFINITY;
		}

	float r = (float) uLengthA/ (float) (uLengthB + 1); // +1 to prevent div 0
	if (r < 1)
		r = 1/r;

	SCORE Score;
	if (g_bTermGaps4 && !bLockLeft && !bLockRight && r > 1.2)
		Score = GlobalAlign4(PA, uLengthA, PB, uLengthB, Path);
	else
		Score = GlobalAlign(PA, uLengthA, PB, uLengthB, Path);

	AlignTwoMSAsGivenPath(Path, msa1, msa2, msaOut);

	delete[] PA;
	delete[] PB;

	return Score;
	}

SCORE GlobalAlign4(ProfPos *PA, unsigned uLengthA, ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	SCORE Score11;
	SCORE Score10;
	SCORE Score01;
	SCORE Score00;

	SCORE Match11;
	SCORE Match10;
	SCORE Match01;
	SCORE Match00;

	PWPath Path11;
	PWPath Path10;
	PWPath Path01;
	PWPath Path00;

	Score11 = GlobalAlignXX(1, 1, PA, uLengthA, PB, uLengthB, Path11, &Match11);
	Score10 = GlobalAlignXX(1, 0, PA, uLengthA, PB, uLengthB, Path10, &Match10);
	Score01 = GlobalAlignXX(0, 1, PA, uLengthA, PB, uLengthB, Path01, &Match01);
	Score00 = GlobalAlignXX(0, 0, PA, uLengthA, PB, uLengthB, Path00, &Match00);

	SCORE MaxMatch = Match11;
	int Best = 1;
	if (Match10 > MaxMatch)
		{
		Best = 2;
		MaxMatch = Match10;
		}
	if (Match01 > MaxMatch)
		{
		Best = 3;
		MaxMatch = Match01;
		}
	if (Match00 > MaxMatch)
		{
		Best = 4;
		MaxMatch = Match00;
		}

	switch (Best)
		{
	case 1:
		Path.Copy(Path11);
		return Score11;
	case 2:
		Path.Copy(Path10);
		return Score10;
	case 3:
		Path.Copy(Path01);
		return Score01;
	case 4:
		Path.Copy(Path00);
		return Score00;
		}
	Quit("GlobalAlign4 invalid Best");
	return MINUS_INFINITY;
	}

static SCORE GlobalAlignXX(int NPenalty, int CPenalty, ProfPos *PA,
  unsigned uLengthA, ProfPos *PB, unsigned uLengthB, PWPath &Path, SCORE *ptrMatch)
	{
	bool bALonger = (uLengthA > uLengthB);
	ProfPos &PPN = (bALonger ? PA[0] : PB[0]);
	ProfPos &PPC = (bALonger ? PA[uLengthA - 1] : PB[uLengthB - 1]);

	if (NPenalty)
		PPN.m_scoreGapOpen = g_scoreGapOpen;
	else
		PPN.m_scoreGapOpen = 0;

	if (CPenalty)
		PPC.m_scoreGapClose = g_scoreGapOpen;
	else
		PPC.m_scoreGapClose = 0;

	SCORE Score = GlobalAlign(PA, uLengthA, PB, uLengthB, Path);
	*ptrMatch = MatchScore(PA, uLengthA, PB, uLengthB, Path);
	return Score;
	}

static SCORE MatchScore(ProfPos *PA, unsigned uLengthA, ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	SCORE Score = 0;
	const unsigned uEdgeCount = Path.GetEdgeCount();
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = Path.GetEdge(uEdgeIndex);
		if (Edge.cType != 'M')
			continue;
		const ProfPos &PPA = PA[Edge.uPrefixLengthA - 1];
		const ProfPos &PPB = PB[Edge.uPrefixLengthB - 1];
		Score += ScoreProfPos2(PPA, PPB);
		}
	return Score;
	}
