#include "muscle.h"
#include <ctype.h>

/***
Nucleotide alphabet is AGCUTRYN
	T and U are equivalent in scoring etc.
	R = purine (G or A)
	Y = pyrimidine (C or T/U)
	N = any
***/

unsigned g_CharToLetter[MAX_CHAR];
unsigned g_CharToLetterEx[MAX_CHAR];

char g_LetterToChar[MAX_ALPHA];
char g_LetterExToChar[MAX_ALPHA_EX];

char g_UnalignChar[MAX_CHAR];
char g_AlignChar[MAX_CHAR];

bool g_IsWildcardChar[MAX_CHAR];
bool g_IsResidueChar[MAX_CHAR];

ALPHA g_Alpha = ALPHA_Undefined;
unsigned g_AlphaSize = 0;

#define Res(c, Letter)												\
	{																\
	const unsigned char Upper = (unsigned char) toupper(c);			\
	const unsigned char Lower = (unsigned char) tolower(c);			\
	g_CharToLetter[Upper] = Letter;									\
	g_CharToLetter[Lower] = Letter;									\
	g_CharToLetterEx[Upper] = Letter;								\
	g_CharToLetterEx[Lower] = Letter;								\
	g_LetterToChar[Letter] = Upper;									\
	g_LetterExToChar[Letter] = Upper;								\
	g_IsResidueChar[Upper] = true;									\
	g_IsResidueChar[Lower] = true;									\
	g_AlignChar[Upper] = Upper;										\
	g_AlignChar[Lower] = Upper;										\
	g_UnalignChar[Upper] = Lower;									\
	g_UnalignChar[Lower] = Lower;									\
	}

#define Wild(c, Letter)												\
	{																\
	const unsigned char Upper = (unsigned char) toupper(c);			\
	const unsigned char Lower = (unsigned char) tolower(c);			\
	g_CharToLetterEx[Upper] = Letter;								\
	g_CharToLetterEx[Lower] = Letter;								\
	g_LetterExToChar[Letter] = Upper;								\
	g_IsResidueChar[Upper] = true;									\
	g_IsResidueChar[Lower] = true;									\
	g_AlignChar[Upper] = Upper;										\
	g_AlignChar[Lower] = Upper;										\
	g_UnalignChar[Upper] = Lower;									\
	g_UnalignChar[Lower] = Lower;									\
	g_IsWildcardChar[Lower] = true;									\
	g_IsWildcardChar[Upper] = true;									\
	}

static unsigned GetAlphaSize(ALPHA Alpha)
	{
	switch (Alpha)
		{
	case ALPHA_Amino:
		return 20;

	case ALPHA_Nucleo:
		return 4;
		}
	Quit("Invalid Alpha=%d", Alpha);
	return 0;
	}

static void InitArrays()
	{
	memset(g_CharToLetter, 0xff, sizeof(g_CharToLetter));
	memset(g_CharToLetterEx, 0xff, sizeof(g_CharToLetterEx));

	memset(g_LetterToChar, '?', sizeof(g_LetterToChar));
	memset(g_LetterExToChar, '?', sizeof(g_LetterExToChar));

	memset(g_AlignChar, '?', sizeof(g_UnalignChar));
	memset(g_UnalignChar, '?', sizeof(g_UnalignChar));

	memset(g_IsWildcardChar, 0, sizeof(g_IsWildcardChar));
	}

static void SetGapChar(char c)
	{
	unsigned char u = (unsigned char) c;

	g_CharToLetterEx[u] = AX_GAP;
	g_LetterExToChar[AX_GAP] = u;
	g_AlignChar[u] = u;
	g_UnalignChar[u] = u;
	}

static void SetAlphaNucleo()
	{
	Res('A', NX_A)
	Res('C', NX_C)
	Res('G', NX_G)
	Res('U', NX_U)
	Res('T', NX_T)

	Wild('N', NX_N)
	Wild('R', NX_R)
	Wild('Y', NX_Y)
	}

static void SetAlphaAmino()
	{
	Res('A', AX_A)
	Res('C', AX_C)
	Res('D', AX_D)
	Res('E', AX_E)
	Res('F', AX_F)
	Res('G', AX_G)
	Res('H', AX_H)
	Res('I', AX_I)
	Res('K', AX_K)
	Res('L', AX_L)
	Res('M', AX_M)
	Res('N', AX_N)
	Res('P', AX_P)
	Res('Q', AX_Q)
	Res('R', AX_R)
	Res('S', AX_S)
	Res('T', AX_T)
	Res('V', AX_V)
	Res('W', AX_W)
	Res('Y', AX_Y)

	Wild('B', AX_B)
	Wild('X', AX_X)
	Wild('Z', AX_Z)
	}

void SetAlpha(ALPHA Alpha)
	{
	InitArrays();

	SetGapChar('.');
	SetGapChar('-');

	switch (Alpha)
		{
	case ALPHA_Amino:
		SetAlphaAmino();
		break;

	case ALPHA_Nucleo:
		SetAlphaNucleo();
		break;

	default:
		Quit("Invalid Alpha=%d", Alpha);
		}

	g_AlphaSize = GetAlphaSize(Alpha);
	g_Alpha = Alpha;
	}

char GetWildcardChar()
	{
	switch (g_Alpha)
		{
	case ALPHA_Amino:
		return 'X';

	case ALPHA_Nucleo:
		return 'N';

	default:
		Quit("Invalid Alpha=%d", g_Alpha);
		}
	return '?';
	}

bool IsNucleo(char c)
	{
	return strchr("ACGTURYNacgturyn", c) != 0;
	}
