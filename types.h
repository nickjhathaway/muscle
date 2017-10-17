#ifndef types_h
#define types_h

typedef unsigned char byte;
typedef unsigned short ushort;

typedef float SCOREMATRIX[32][32];
typedef SCOREMATRIX *PTR_SCOREMATRIX;

class MSA;
class Seq;
class ClusterTree;
class DistFunc;
class TextFile;
class PWPath;
class Tree;
class SeqVect;
class DistCalc;

struct ProgNode;
struct ProfPos;

enum EXIT
	{
	EXIT_Success = 0,
	EXIT_NotStarted = 1,
	EXIT_FatalError = 2,
	EXIT_Except = 3,
	};

enum NODECMP
	{
	NODECMP_Undefined = 0,
	NODECMP_Same = 0,		// equivalent to node in old tree
	NODECMP_Diff = 1,		// equivalent & parent is changed
	NODECMP_Changed = 2		// no equivalent node in old tree
	};

// Declare enums using macro hacks (see enums.h).
#define s(t)	enum t { t##_Undefined = 0,
#define c(t, x)	t##_##x,
#define e(t)	};
#include "enums.h"

// Declare conversion function XXXToStr(XXX x)
// for each enum type XXX.
#define	s(t)	const char *t##ToStr(t x);
#define c(t, x)	/* empty */
#define e(t)	/* empty */
#include "enums.h"

// Declare conversion function StrToXXX(const char *Str)
// for each enum type XXX.
#define	s(t)	t StrTo##t(const char *Str);
#define c(t, x)	/* empty */
#define e(t)	/* empty */
#include "enums.h"

const char *BoolToStr(bool b);
const char *SecsToStr(unsigned long Secs);

#endif // types_h
