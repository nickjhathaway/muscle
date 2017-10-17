#include "muscle.h"
#include "seqvect.h"
#include "textfile.h"
#include "msa.h"

const size_t MAX_FASTA_LINE = 16000;

SeqVect::~SeqVect()
	{
	Clear();
	}

void SeqVect::Clear()
	{
	for (size_t n = 0; n < size(); ++n)
		delete (*this)[n];
	}

void SeqVect::ToFASTAFile(TextFile &File) const
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->ToFASTAFile(File);
		}
	}

void SeqVect::FromFASTAFile(TextFile &File)
	{
	Clear();

	char szLine[MAX_FASTA_LINE];
	Seq *ptrSeq = 0;
	for (;;)
		{
		bool bEof = File.GetLine(szLine, sizeof(szLine));
		if (bEof)
			return;
		if ('>' == szLine[0])
			{
			ptrSeq = new Seq;
			push_back(ptrSeq);
			size_t n = strlen(szLine);
			if (1 == n)
				Quit("Missing annotation following '>' in FASTA file %s line %u",
				File.GetFileName(), File.GetLineNr());

			ptrSeq->SetName(szLine + 1);
			}
		else
			{
			const char *ptrChar = szLine;
			while (char c = *ptrChar++)
				{
				if (isspace(c))
					continue;
				if (IsGap(c))
					continue;
				if (!IsValidAminoEx(c))
					{
					if (isprint(c))
						{
						Warning("Invalid amino acid '%c' in FASTA file %s line %d, replaced by 'X'",
						  c, File.GetFileName(), File.GetLineNr());
						c = 'X';
						}
					else
						Quit("Invalid byte hex %02x in FASTA file %s line %d",
						  (unsigned char) c, File.GetFileName(), File.GetLineNr());
					}
				c = toupper(c);
				if (0 == ptrSeq)
					Quit("Invalid FASTA file, missing '>'");
				ptrSeq->push_back(c);
				}
			}
		}
	}

void SeqVect::PadToMSA(MSA &msa)
	{
	unsigned uSeqCount = Length();
	if (0 == uSeqCount)
		{
		msa.Clear();
		return;
		}

	unsigned uLongestSeqLength = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		unsigned uColCount = ptrSeq->Length();
		if (uColCount > uLongestSeqLength)
			uLongestSeqLength = uColCount;
		}
	msa.SetSize(uSeqCount, uLongestSeqLength);
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		msa.SetSeqName(uSeqIndex, ptrSeq->GetName());
		unsigned uColCount = ptrSeq->Length();
		unsigned uColIndex;
		for (uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			char c = ptrSeq->at(uColIndex);
			msa.SetChar(uSeqIndex, uColIndex, c);
			}
		while (uColIndex < uLongestSeqLength)
			msa.SetChar(uSeqIndex, uColIndex++, '.');
		}
	}

void SeqVect::Copy(const SeqVect &rhs)
	{
	clear();
	unsigned uSeqCount = rhs.Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = rhs.at(uSeqIndex);
		Seq *ptrSeqCopy = new Seq;
		ptrSeqCopy->Copy(*ptrSeq);
		push_back(ptrSeqCopy);
		}
	}

void SeqVect::StripGaps()
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->StripGaps();
		}
	}

void SeqVect::StripGapsAndWhitespace()
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->StripGapsAndWhitespace();
		}
	}

void SeqVect::ToUpper()
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->ToUpper();
		}
	}

bool SeqVect::FindName(const char *ptrName, unsigned *ptruIndex) const
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const Seq *ptrSeq = at(uSeqIndex);
		if (0 == stricmp(ptrSeq->GetName(), ptrName))
			{
			*ptruIndex = uSeqIndex;
			return true;
			}
		}
	return false;
	}

void SeqVect::AppendSeq(const Seq &s)
	{
	Seq *ptrSeqCopy = new Seq;
	ptrSeqCopy->Copy(s);
	push_back(ptrSeqCopy);
	}

void SeqVect::LogMe() const
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->LogMe();
		}
	}

const char *SeqVect::GetSeqName(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < size());
	const Seq *ptrSeq = at(uSeqIndex);
	return ptrSeq->GetName();
	}

unsigned SeqVect::GetSeqId(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < size());
	const Seq *ptrSeq = at(uSeqIndex);
	return ptrSeq->GetId();
	}

Seq &SeqVect::GetSeq(unsigned uSeqIndex)
	{
	assert(uSeqIndex < size());
	return *at(uSeqIndex);
	}

const Seq &SeqVect::GetSeq(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < size());
	return *at(uSeqIndex);
	}

void SeqVect::SetSeqId(unsigned uSeqIndex, unsigned uId)
	{
	assert(uSeqIndex < size());
	Seq *ptrSeq = at(uSeqIndex);
	return ptrSeq->SetId(uId);
	}
