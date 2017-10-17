#include "muscle.h"
#include <stdio.h>
#include <ctype.h>
#include "msa.h"
#include "textfile.h"

static bool StripWeightFromName(char *strName, WEIGHT &w)
	{
	w = BTInsane;
	char *ptrLastBar = strrchr(strName, '|');
	if (0 == ptrLastBar)
		return false;
	if (0 != strncmp(ptrLastBar, "|weight=", 8))
		return false;
	w = (WEIGHT) atof(ptrLastBar+8);
	*ptrLastBar = 0;
	return true;
	}

void MSA::FromFASTAFile(TextFile &File, bool bFirstOnly)
	{
	Clear();

// Arbitrary choice for longest allowed input line
	const unsigned MAX_FASTA_LINE = 16000;
	char szLine[MAX_FASTA_LINE+1];

// Make two passes through the file, one to determine the
// number of sequences and the length of the sequence so
// that we can allocate exactly the right amount of
// memory, the second pass to store the data.
	TEXTFILEPOS Pos = File.GetPos();

// First pass.
	unsigned uLettersPerLine = 80;
	bool bLettersPerLineKnown = false;
	unsigned uSeqCount = 0;
	unsigned uColCount = 0;
	unsigned uCurrentSeqLength = 0;
	bool bAminoSeen = false;
	bool bGapSeen = false;
	for (;;)
		{
		bool bEOF = File.GetTrimLine(szLine, sizeof(szLine));
		if (bEOF || !strncmp(szLine, "//", 2))
			{
			if (uSeqCount > 1)
				{
				if (uColCount != uCurrentSeqLength)
					Quit("FASTA file %s, last sequence length=%d, others=%d",
					  File.GetFileName(), uCurrentSeqLength, uColCount);
				}
			else
				uColCount = uCurrentSeqLength;
			break;
			}
		if ('>' == szLine[0])
			{
			if (uSeqCount > 1)
				{
				if (uColCount != uCurrentSeqLength)
					Quit("FASTA file %s(%u), sequence %d is length=%d, others=%d",
					  File.GetFileName(), File.GetLineNr(), uSeqCount,
					  uCurrentSeqLength, uColCount);
				}
			else
				uColCount = uCurrentSeqLength;
			if (bFirstOnly && 1 == uSeqCount)
				goto SecondPass;
			++uSeqCount;
			uCurrentSeqLength = 0;
			}
		else
			{
			if (!bAminoSeen)
				bAminoSeen = StrHasAmino(szLine);
			if (!bGapSeen)
				bGapSeen = StrHasGap(szLine);
			unsigned uLetterCount = (unsigned) strlen(szLine);
			if (bLettersPerLineKnown)
				{
				if (uLetterCount > uLettersPerLine)
					Quit("FASTA file %s(%d), seq %d: Expecting <=%d letters, got %d",
					  File.GetFileName(), File.GetLineNr(), uSeqCount,
					  uLettersPerLine, uLetterCount);
				}
			else
				{
				uLettersPerLine = uLetterCount;
				bLettersPerLineKnown = true;
				}
			uCurrentSeqLength += uLetterCount;
			}
		}

SecondPass:
	if (bFirstOnly && uSeqCount > 1)
		Quit("Bug in MSA::FromFASTAFile (bFirstOnly)");

	SetSize(uSeqCount, uColCount);
	m_uColCount = uColCount;

	File.SetPos(Pos);
	unsigned uLinesPerSeq = (uColCount - 1)/uLettersPerLine + 1;
	unsigned uSeqIndex;
	for (uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		bool bEOF = File.GetTrimLine(szLine, sizeof(szLine));
		if (bEOF)
			Quit("Unexpected EOF in FASTA file %s (first line of seq %d)",
			  File.GetFileName(), uSeqIndex);
		if ('>' != szLine[0])
			Quit("Expected annotation line");

		m_szNames[uSeqIndex] = strdup(szLine+1);
		m_Weights[uSeqIndex] = BTInsane;
		if (0 == strlen(m_szNames[uSeqIndex]))
			sprintf(m_szNames[uSeqIndex], "[%d]", uSeqIndex);
		
		for (unsigned uLineIndex = 0; uLineIndex < uLinesPerSeq; ++uLineIndex)
			{
			char *pLineBase = m_szSeqs[uSeqIndex] + uLineIndex*uLettersPerLine;
		// Add one to expected number of chars to
		// allow for nul byte to terminate string
			unsigned uCharCount;
			if (uLinesPerSeq - 1 == uLineIndex)
				uCharCount = (uColCount - 1)%uLettersPerLine + 2;
			else
				uCharCount = uLettersPerLine + 1;
			bEOF = File.GetTrimLine(pLineBase, uCharCount);

			if (bEOF)
				Quit("Unexpected EOF in FASTA file %s, seq %d line %d",
				  File.GetFileName(), uSeqIndex+1, uLineIndex+1);

			if (uLinesPerSeq - 1 != uLineIndex &&
			  strlen(pLineBase) != uLettersPerLine)
				Quit("FASTA line should be %d characters file %s seq %d line %d",
				  uLettersPerLine, File.GetFileName(), uSeqIndex+1, uLineIndex+1);
			}
		}

	if (!bFirstOnly)
		{
	// Should have eaten all lines in the file...
		bool bEOF = File.GetTrimLine(szLine, sizeof(szLine));
		if (!bEOF && 0 != strcmp("//", szLine))
			Quit("Expected EOF in FASTA file %s", File.GetFileName());
		}

	for (uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		for (unsigned n = 0; n < GetColCount(); ++n)
			{
			char c = m_szSeqs[uSeqIndex][n];
			if (!::IsValidAminoEx(c) && !::IsGap(c))
				{
				Warning("Warning -- invalid character '%c' replaced by 'X'", c);
				c = 'X';
				}
			m_szSeqs[uSeqIndex][n] = c;
			}
	}

void MSA::ToFASTAFile(TextFile &File) const
	{
	const unsigned uColCount = GetColCount();
	assert(uColCount > 0);
	const unsigned uLettersPerLine = 69;
	const unsigned uLinesPerSeq = (GetColCount() - 1)/uLettersPerLine + 1;
	const unsigned uSeqCount = GetSeqCount();

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		File.PutString(">");
		File.PutString(GetSeqName(uSeqIndex));
		File.PutString("\n");

		unsigned n = 0;
		for (unsigned uLine = 0; uLine < uLinesPerSeq; ++uLine)
			{
			unsigned uLetters = uColCount - uLine*uLettersPerLine;
			if (uLetters > uLettersPerLine)
				uLetters = uLettersPerLine;
			for (unsigned i = 0; i < uLetters; ++i)
				{
				char c = GetChar(uSeqIndex, n);
				File.PutChar(c);
				++n;
				}
			File.PutChar('\n');
			}
		}
	}
