#include "muscle.h"
#include "msa.h"
#include "profile.h"
#include "pwpath.h"

SCORE GlobalAlign4(ProfPos *PA, unsigned uLengthA, ProfPos *PB,
  unsigned uLengthB, PWPath &Path);

SCORE AlignTwoProfs(
  const ProfPos *PA, unsigned uLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uLengthB, WEIGHT wB,
  PWPath &Path, ProfPos **ptrPout, unsigned *ptruLengthOut)
	{
	assert(uLengthA < 100000);
	assert(uLengthB < 100000);

	float r = (float) uLengthA/ (float) (uLengthB + 1); // +1 to prevent div 0
	if (r < 1)
		r = 1/r;

	SCORE Score;
	if (g_bTermGaps4 && r > 1.2)
		Score = GlobalAlign4((ProfPos *) PA, uLengthA, (ProfPos *) PB, uLengthB, Path);
	else
		Score = GlobalAlign(PA, uLengthA, PB, uLengthB, Path);

	AlignTwoProfsGivenPath(Path, PA, uLengthB, wA/(wA + wB), PB, uLengthB, wB/(wA + wB),
	  ptrPout, ptruLengthOut);

	if (ALPHA_Amino == g_Alpha)
		Hydro(*ptrPout, *ptruLengthOut);

	return Score;
	}
