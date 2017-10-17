#ifndef ObjScore_h
#define ObjScore_h

SCORE ObjScore(const MSA &msa, const unsigned SeqIndexes1[],
  unsigned uSeqCount1, const unsigned SeqIndexes2[], unsigned uSeqCount2);

SCORE ObjScoreIds(const MSA &msa, const unsigned Ids1[],
  unsigned uCount1, const unsigned Ids2[], unsigned uCount2);

SCORE ObjScoreDP(const MSA &msa1, const MSA &msa2, SCORE MatchScore[] = 0);
SCORE ObjScorePS(const MSA &msa, SCORE MatchScore[] = 0);
SCORE ObjScoreSP(const MSA &msa, SCORE MatchScore[] = 0);
SCORE ObjScoreXP(const MSA &msa, const MSA &msa2);
SCORE ObjScoreSPDimer(const MSA &msa);

SCORE DiffObjScore(
  const MSA &msa1, const PWPath &Path1, const unsigned Edges1[], unsigned uEdgeCount1, 
  const MSA &msa2, const PWPath &Path2, const unsigned Edges2[], unsigned uEdgeCount2);

SCORE ObjScoreEdges(const MSA &msa, const PWPath &path, const unsigned Edges[],
  unsigned uEdgeCount);

#endif // ObjScore_h
