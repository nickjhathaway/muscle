#include "muscle.h"
#include <stdio.h>
#include <ctype.h>
#include "msa.h"
#include "textfile.h"

const unsigned uCharsPerLine = 50;
const unsigned uCharsPerBlock = 10;

void MSA::ToMSFFile(TextFile &File, const char *ptrComment) const
	{
	File.PutString("PileUp\n");

	if (0 != ptrComment)
		File.PutFormat("Comment: %s\n", ptrComment);
	else
		File.PutString("\n");

	File.PutFormat("MSF: %u\n\n", GetColCount());

	int iLongestNameLength = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		{
		const char *ptrName = GetSeqName(uSeqIndex);
		int i = (int) strlen(ptrName);
		if (i > iLongestNameLength)
			iLongestNameLength = i;
		}
	
	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		{
		const char *ptrName = GetSeqName(uSeqIndex);
		File.PutFormat(" Name: %s", ptrName);
		int i = iLongestNameLength - (int) strlen(ptrName);
		while (i > 0)
			{
			File.PutString(" ");
			--i;
			}
		File.PutFormat("  Len: %u\n", GetColCount());
		}
	File.PutString("\n//\n");
	if (0 == GetColCount())
		return;

	unsigned uLineCount = (GetColCount() - 1)/uCharsPerLine + 1;
	for (unsigned uLineIndex = 0; uLineIndex < uLineCount; ++uLineIndex)
		{
		File.PutString("\n");
		unsigned uStartColIndex = uLineIndex*uCharsPerLine;
		unsigned uEndColIndex = uStartColIndex + uCharsPerLine - 1;
		if (uEndColIndex >= GetColCount())
			uEndColIndex = GetColCount() - 1;
		for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
			{
			const char *ptrName = GetSeqName(uSeqIndex);
			File.PutFormat("%s   ", ptrName);
			int i = iLongestNameLength - (int) strlen(ptrName);
			while (i > 0)
				{
				File.PutString(" ");
				--i;
				}
			for (unsigned uColIndex = uStartColIndex; uColIndex <= uEndColIndex;
			  ++uColIndex)
				{
				if (0 == uColIndex%uCharsPerBlock)
					File.PutString(" ");
				File.PutFormat("%c", GetChar(uSeqIndex, uColIndex));
				}
			File.PutString("\n");
			}
		}
	}

void MSA::FromMSFFile(TextFile &File)
	{
	Clear();

// Arbitrary choice for longest allowed input line
	const unsigned MAX_MSF_LINE = 4096;
	char szLine[MAX_MSF_LINE+1];

// First pass: count sequences
	unsigned uSeqCount = 0;
	for (;;)
		{
		bool bEOF = File.GetTrimLine(szLine, sizeof(szLine));
		if (bEOF)
			Quit("Unexpected EOF in MSF file %s", File.GetFileName());
		if (strncmp(szLine, "Name:", 5) == 0)
			++uSeqCount;
		else if (strncmp(szLine, "//", 2) == 0)
			break;
		}
	File.Rewind();

// Second pass
	SetSize(uSeqCount, 1);
	unsigned *SeqLength = new unsigned[uSeqCount];
	memset(SeqLength, 0, uSeqCount*sizeof(unsigned));
	unsigned uSeqIndex = 0;
	for (;;)
		{
		bool bEOF = File.GetTrimLine(szLine, sizeof(szLine));
		if (bEOF)
			Quit("Unexpected EOF in MSF file %s", File.GetFileName());
		if (strncmp(szLine, "Name:", 5) == 0)
			{
			char *ptrName = szLine + 5;
			for (;;)
				{
				char c = *ptrName;
				if (0 == c)
					Quit("Unexpected EOL in MSF file %s\n", File.GetFileName());
				if (!isspace(c))
					break;
				++ptrName;
				}
			char *ptrAfterName = ptrName + 1;
			for (;;)
				{
				char c = *ptrAfterName;
				if (0 == c)
					break;
				if (isspace(c))
					{
					*ptrAfterName = 0;
					break;
					}
				++ptrAfterName;
				}
			SetSeqName(uSeqIndex, ptrName);
			++uSeqIndex;
			}
		else if (strncmp(szLine, "//", 2) == 0)
			break;
		}
	for (;;)
		{
		bool bEOF = File.GetTrimLine(szLine, sizeof(szLine));
		if (bEOF)
			break;
		if (0 == strlen(szLine))
			continue;
		char *ptrAfterName = szLine;
		for (;;)
			{
			char c = *ptrAfterName;
			if (0 == c)
				Quit("Unexpected EOL in MSF file %s\n", File.GetFileName());
			if (isspace(c))
				break;
			++ptrAfterName;
			}
		*ptrAfterName = 0;
		unsigned uSeqIndex;
		bool bFound = GetSeqIndex(szLine, &uSeqIndex);
		assert(bFound);
		char *ptrLetters = ptrAfterName + 1;
		unsigned uCol = SeqLength[uSeqIndex];
		while (char c = *ptrLetters++)
			{
			if (' ' != c)
				SetChar(uSeqIndex, uCol++, c);
			}
		SeqLength[uSeqIndex] = uCol;
		}
	delete[] SeqLength;
	}
