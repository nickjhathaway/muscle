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
