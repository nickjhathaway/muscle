#include "muscle.h"
#include "msa.h"
#include "objscore.h"
#include "profile.h"

SCORE ObjScore(const MSA &msa, const unsigned SeqIndexes1[],
  unsigned uSeqCount1, const unsigned SeqIndexes2[], unsigned uSeqCount2)
	{
	const unsigned uSeqCount = msa.GetSeqCount();

	OBJSCORE OS = g_ObjScore;
	if (g_ObjScore == OBJSCORE_SPM)
		{
        if (uSeqCount <= 100)
			OS = OBJSCORE_XP;
		else
			OS = OBJSCORE_SPF;
		}

	MSA msa1;
	MSA msa2;

	switch (OS)
		{
	case OBJSCORE_DP:
	case OBJSCORE_XP:
		MSAFromSeqSubset(msa, SeqIndexes1, uSeqCount1, msa1);
		MSAFromSeqSubset(msa, SeqIndexes2, uSeqCount2, msa2);

		SetMSAWeightsMuscle(msa1);
		SetMSAWeightsMuscle(msa2);
		break;

	case OBJSCORE_SP:
	case OBJSCORE_SPF:
	case OBJSCORE_PS:
	// Yuck -- casting away const (design flaw)
		SetMSAWeightsMuscle((MSA &) msa);
		break;
		}

	switch (OS)
		{
	case OBJSCORE_SP:
		return ObjScoreSP(msa);

	case OBJSCORE_DP:
		return ObjScoreDP(msa1, msa2);

	case OBJSCORE_XP:
		return ObjScoreXP(msa1, msa2);

	case OBJSCORE_PS:
		return ObjScorePS(msa);

	case OBJSCORE_SPF:
		return ObjScoreSPDimer(msa);
		}

	Quit("Invalid g_ObjScore=%d", g_ObjScore);
	return 0;
	}

SCORE ObjScoreIds(const MSA &msa, const unsigned Ids1[],
  unsigned uCount1, const unsigned Ids2[], unsigned uCount2)
	{
	unsigned *SeqIndexes1 = new unsigned[uCount1];
	unsigned *SeqIndexes2 = new unsigned[uCount2];

	for (unsigned n = 0; n < uCount1; ++n)
		SeqIndexes1[n] = msa.GetSeqIndex(Ids1[n]);

	for (unsigned n = 0; n < uCount2; ++n)
		SeqIndexes2[n] = msa.GetSeqIndex(Ids2[n]);

	SCORE dObjScore = ObjScore(msa, SeqIndexes1, uCount1, SeqIndexes2, uCount2);

	delete[] SeqIndexes1;
	delete[] SeqIndexes2;

	return dObjScore;
	}
