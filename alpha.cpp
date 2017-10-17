#include "muscle.h"
#include <ctype.h>

// The 21st amino acid, Selenocysteine (Sel, U), is treated
// as an X because I have no substitution probabilities for it.

const char szAmino[]	= "ACDEFGHIKLMNPQRSTVWY";
const char szAminoEx[]	= "ACDEFGHIKLMNPQRSTVWYBZX";

static unsigned uReverseAmino[256];
unsigned CharToLetterAminoEx[256];

static void InitCharToLetterAminoEx()
	{
	memset(CharToLetterAminoEx, 0xff, sizeof(CharToLetterAminoEx));

	CharToLetterAminoEx['a'] = AX_A;
	CharToLetterAminoEx['c'] = AX_C;
	CharToLetterAminoEx['d'] = AX_D;
	CharToLetterAminoEx['e'] = AX_E;
	CharToLetterAminoEx['f'] = AX_F;
	CharToLetterAminoEx['g'] = AX_G;
	CharToLetterAminoEx['h'] = AX_H;
	CharToLetterAminoEx['i'] = AX_I;
	CharToLetterAminoEx['k'] = AX_K;
	CharToLetterAminoEx['l'] = AX_L;
	CharToLetterAminoEx['m'] = AX_M;
	CharToLetterAminoEx['n'] = AX_N;
	CharToLetterAminoEx['p'] = AX_P;
	CharToLetterAminoEx['q'] = AX_Q;
	CharToLetterAminoEx['r'] = AX_R;
	CharToLetterAminoEx['s'] = AX_S;
	CharToLetterAminoEx['t'] = AX_T;
	CharToLetterAminoEx['v'] = AX_V;
	CharToLetterAminoEx['w'] = AX_W;
	CharToLetterAminoEx['y'] = AX_Y;

	CharToLetterAminoEx['b'] = AX_B;
	CharToLetterAminoEx['z'] = AX_Z;
	CharToLetterAminoEx['x'] = AX_X;

	CharToLetterAminoEx['A'] = AX_A;
	CharToLetterAminoEx['C'] = AX_C;
	CharToLetterAminoEx['D'] = AX_D;
	CharToLetterAminoEx['E'] = AX_E;
	CharToLetterAminoEx['F'] = AX_F;
	CharToLetterAminoEx['G'] = AX_G;
	CharToLetterAminoEx['H'] = AX_H;
	CharToLetterAminoEx['I'] = AX_I;
	CharToLetterAminoEx['K'] = AX_K;
	CharToLetterAminoEx['L'] = AX_L;
	CharToLetterAminoEx['M'] = AX_M;
	CharToLetterAminoEx['N'] = AX_N;
	CharToLetterAminoEx['P'] = AX_P;
	CharToLetterAminoEx['Q'] = AX_Q;
	CharToLetterAminoEx['R'] = AX_R;
	CharToLetterAminoEx['S'] = AX_S;
	CharToLetterAminoEx['T'] = AX_T;
	CharToLetterAminoEx['V'] = AX_V;
	CharToLetterAminoEx['W'] = AX_W;
	CharToLetterAminoEx['Y'] = AX_Y;

	CharToLetterAminoEx['B'] = AX_B;
	CharToLetterAminoEx['Z'] = AX_Z;
	CharToLetterAminoEx['X'] = AX_X;

	CharToLetterAminoEx['.'] = AX_GAP;
	CharToLetterAminoEx['-'] = AX_GAP;
	CharToLetterAminoEx['~'] = AX_GAP;
	}

static bool InitReverseTables()
	{
	InitCharToLetterAminoEx();

	memset(uReverseAmino, 0xff, sizeof(uReverseAmino));

	for (unsigned n = 0; n < sizeof(szAmino); ++n)
		{
		unsigned char c = (unsigned char) szAmino[n];
		if (isalpha(c))
			{
			uReverseAmino[toupper(c)] = n;
			uReverseAmino[tolower(c)] = n;
			}
		else
			uReverseAmino[tolower(c)] = n;
		}

	return true;
	}
static bool bInitDone = InitReverseTables();

char LetterToCharAmino(unsigned uLetter)
	{
	assert(uLetter < 20);
	return szAmino[uLetter];
	}

char LetterToCharAminoEx(unsigned uLetter)
	{
 	assert(uLetter < 23);
	return szAminoEx[uLetter];
	}

unsigned CharToLetterAmino(char c)
	{
	unsigned uLetter = uReverseAmino[(unsigned char) c];
	assert(uLetter < 20);
	return uLetter;
	}

unsigned CharToLetterX(char c)
	{
	return CharToLetterAminoEx[c];
	}

bool IsValidAmino(char c)
	{
//	return 0 != strchr(szAmino, toupper(c));
	return 0xffffffff != uReverseAmino[(unsigned char) c];
	}

bool IsValidAminoEx(char c)
	{
//	return 0 != strchr(szAminoEx, toupper(c));
	return 0xffffffff != CharToLetterAminoEx[(unsigned char) c];
	}

bool IsGap(char c)
	{
	return '-' == c || '.' == c || c == '~';
	}

bool IsNonTerminalGap(char c)
	{
	return '-' == c || '.' == c;
	}

bool IsTerminalGap(char c)
	{
	return '~' == c;
	}

bool IsWildcard(char c)
	{
	c = toupper(c);
	return 'X' == c || 'B' == c || 'Z' == c || 'U' == c;
	}

bool StrHasAmino(const char *Str)
	{
	while (char c = *Str++)
		if (IsValidAmino(c))
			return true;
	return false;
	}

bool StrHasGap(const char *Str)
	{
	while (char c = *Str++)
		if (IsGap(c))
			return true;
	return false;
	}

char UnalignChar(char c)
	{
	if (isalpha(c))
		return tolower(c);
	else if ('-' == c || '.' == c || '~' == c)
		return '.';
	assert(false);
	return 0;
	}

bool IsAlignedChar(char c)
	{
	return (isalpha(c) && isupper(c)) || '-' == c || '~' == c;
	}

bool IsUnalignedChar(char c)
	{
	return !IsAlignedChar(c);
	}
