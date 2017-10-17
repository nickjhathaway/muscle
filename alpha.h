#ifndef	Alpha_h
#define	Alpha_h

unsigned CharToLetterX(char c);
unsigned CharToLetterAmino(char c);
unsigned CharToLetterAminoX(char c);
char LetterToCharAmino(unsigned uLetter);
char LetterToCharAminoEx(unsigned uLetter);
char LetterToCharNucleic(unsigned uLetter);
bool IsValidAmino(char c);
bool IsValidAminoEx(char c);
bool IsGap(char c);
bool IsTerminalGap(char c);
bool IsNonTerminalGap(char c);
bool IsWildcard(char c);
bool StrHasAmino(const char *Str);
bool StrHasGap(const char *Str);
char UnalignChar(char c);
bool IsAlignedChar(char c);
bool IsUnalignedChar(char c);

extern unsigned CharToLetterAminoEx[];

// AX=Amino alphabet with eXtensions (B, Z and X)
enum AX
	{
	AX_A,
	AX_C,
	AX_D,
	AX_E,
	AX_F,
	AX_G,
	AX_H,
	AX_I,
	AX_K,
	AX_L,
	AX_M,
	AX_N,
	AX_P,
	AX_Q,
	AX_R,
	AX_S,
	AX_T,
	AX_V,
	AX_W,
	AX_Y,

	AX_X,	// Unknown

	AX_B,	// D or N
	AX_Z,	// E or Q

	AX_GAP,
	};

const int AX_COUNT = 24;

#endif	// Alpha_h
