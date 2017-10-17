#include "muscle.h"
#include "pwpath.h"
#include "timing.h"
#include "textfile.h"
#include "msa.h"
#include "profile.h"

SCORE GlobalAlign(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	if (g_bDiags)
		return GlobalAlignDiags(PA, uLengthA, PB, uLengthB, Path);
	else
		return GlobalAlignNoDiags(PA, uLengthA, PB, uLengthB, Path);
	}

SCORE GlobalAlignNoDiags(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	if (g_bDimer)
		return GlobalAlignDimer(PA, uLengthA, PB, uLengthB, Path);

	switch (g_PPScore)
		{
	case PPSCORE_LE:
		return GlobalAlignLE(PA, uLengthA, PB, uLengthB, Path);

	case PPSCORE_SP:
	case PPSCORE_SV:
		return GlobalAlignSP(PA, uLengthA, PB, uLengthB, Path);

	case PPSCORE_SPN:
		return GlobalAlignSimple(PA, uLengthA, PB, uLengthB, Path);
		}

	Quit("Invalid PP score (GlobalAlignNoDiags)");
	return 0;
	}
