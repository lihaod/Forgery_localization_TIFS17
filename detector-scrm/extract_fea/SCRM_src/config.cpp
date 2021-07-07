#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>

#ifndef CONFIG_H_
#define CONFIG_H_

class Config
{
public:
	bool enable_s, enable_c;
	bool verbose;
	int T_s, T_c;
	int order_s, order_c;
	int *SPAMsymmCoord_s, *SPAMsymmCoord_c;
	int SPAMsymmDim_s, SPAMsymmDim_c;
	int *MINMAXsymmCoord_s, *MINMAXsymmCoord_c;
	int MINMAXsymmDim_s, MINMAXsymmDim_c;
	int numInMap;
	bool mergeSpams;
	bool eraseLSB;
	bool parity;
	bool roundup5;          // if roundup5==true, round(0.5)=1, else round(0.5)=0. In conventional SRM roundup5=false.
	bool mapMode;
	std::vector<float> qs_1;
	std::vector<float> qs_others;

	Config(bool verbose, bool enable_s, bool enable_c, int T_s, int T_c, int order_s, int order_c, std::vector<float> qs_1, std::vector<float> qs_others, bool symmSign, bool symmReverse, bool symmMinMax, bool mergeSpams, bool eraseLSB, bool parity, bool roundup5, bool mapMode = true, int numInMap = 1)
	{
		this->verbose = verbose;
		this->enable_s = enable_s;
		this->enable_c = enable_c;
		this->T_s = T_s;
		this->order_s = order_s;
		this->T_c = T_c;
		this->order_c = order_c;
		this->mergeSpams = mergeSpams;
		this->eraseLSB = eraseLSB;
		this->parity = parity;
		this->qs_1 = qs_1;
		this->qs_others = qs_others;
		this->roundup5 = roundup5;
		this->mapMode = mapMode;
		this->numInMap = numInMap;

		if (enable_s)
			GetSymmCoords_s(symmSign, symmReverse, symmMinMax);
		if (enable_c)
			GetSymmCoords_c(symmSign, symmReverse, symmMinMax);
	}

	~Config()
	{
		if (enable_s)
		{
			delete[] SPAMsymmCoord_s;
			delete[] MINMAXsymmCoord_s;
		}
		if (enable_c)
		{
			delete[] SPAMsymmCoord_c;
			delete[] MINMAXsymmCoord_c;
		}
	}

private:

	void GetSymmCoords_s(bool symmSign, bool symmReverse, bool symmMinMax)
	{
		// Preparation of inSymCoord matrix for co-occurrence and symmetrization
		int B = 2 * this->T_s + 1;
		int fullDim = (int)std::pow((float)B, this->order_s);

		int alreadyUsed;

		// MINMAX
		alreadyUsed = 0;
		MINMAXsymmCoord_s = new int[2 * fullDim]; // [0, fullDim-1] = min; [fullDim, 2*fullDim-1] = max
		for (int i = 0; i<2 * fullDim; i++) MINMAXsymmCoord_s[i] = -1;

		for (int numIter = 0; numIter < fullDim; numIter++)
		{
			if (MINMAXsymmCoord_s[numIter] == -1)
			{
				int coordReverse = 0;
				int num = numIter;
				for (int i = 0; i<this->order_s; i++)
				{
					coordReverse += (num % B) * ((int)std::pow((float)B, order_s - i - 1));
					num = num / B;
				}
				// To the same bin: min(X), max(-X), min(Xreverse), max(-Xreverse)
				if (MINMAXsymmCoord_s[numIter] == -1)
				{
					MINMAXsymmCoord_s[numIter] = alreadyUsed; // min(X)
					if (symmMinMax) MINMAXsymmCoord_s[2 * fullDim - numIter - 1] = alreadyUsed; // max(-X)
					if (symmReverse) MINMAXsymmCoord_s[coordReverse] = alreadyUsed; // min(Xreverse)
					if ((symmMinMax) && (symmReverse)) MINMAXsymmCoord_s[2 * fullDim - coordReverse - 1] = alreadyUsed; // max(-Xreverse)
					alreadyUsed++;
				}
			}
		}
		for (int numIter = 0; numIter < fullDim; numIter++)
		{
			if (MINMAXsymmCoord_s[fullDim + numIter] == -1)
			{
				int coordReverse = 0;
				int num = numIter;
				for (int i = 0; i<this->order_s; i++)
				{
					coordReverse += (num % B) * ((int)std::pow((float)B, order_s - i - 1));
					num = num / B;
				}
				// To the same bin: max(X), min(-X), max(Xreverse), min(-Xreverse)
				if (MINMAXsymmCoord_s[fullDim + numIter] == -1)
				{
					MINMAXsymmCoord_s[fullDim + numIter] = alreadyUsed; // max(X)
					if (symmMinMax) MINMAXsymmCoord_s[fullDim - numIter - 1] = alreadyUsed; // min(-X)
					if (symmReverse) MINMAXsymmCoord_s[fullDim + coordReverse] = alreadyUsed; // max(Xreverse)
					if ((symmMinMax) && (symmReverse)) MINMAXsymmCoord_s[fullDim - coordReverse - 1] = alreadyUsed; // min(-Xreverse)
					alreadyUsed++;
				}
			}
		}
		MINMAXsymmDim_s = alreadyUsed;

		// SPAM
		alreadyUsed = 0;
		SPAMsymmCoord_s = new int[fullDim];
		for (int i = 0; i<fullDim; i++) SPAMsymmCoord_s[i] = -1;
		for (int numIter = 0; numIter < fullDim; numIter++)
		{
			if (SPAMsymmCoord_s[numIter] == -1)
			{
				int coordReverse = 0;
				int num = numIter;
				for (int i = 0; i<this->order_s; i++)
				{
					coordReverse += (num % B) * ((int)std::pow((float)B, order_s - i - 1));
					num = num / B;
				}
				// To the same bin: X, -X, Xreverse, -Xreverse
				SPAMsymmCoord_s[numIter] = alreadyUsed; // X
				if (symmSign) SPAMsymmCoord_s[fullDim - numIter - 1] = alreadyUsed; // -X
				if (symmReverse) SPAMsymmCoord_s[coordReverse] = alreadyUsed; // Xreverse
				if ((symmSign) && (symmReverse)) SPAMsymmCoord_s[fullDim - coordReverse - 1] = alreadyUsed; // -Xreverse
				alreadyUsed++;
			}
		}
		SPAMsymmDim_s = alreadyUsed;
		// In order to have the same order of the features as the matlab SRM - shift +1
		for (int i = 0; i<fullDim; i++)
		{
			if (SPAMsymmCoord_s[i] == alreadyUsed - 1) SPAMsymmCoord_s[i] = 0;
			else SPAMsymmCoord_s[i]++;
		}
	}

	void GetSymmCoords_c(bool symmSign, bool symmReverse, bool symmMinMax)
	{
		// Preparation of inSymCoord matrix for co-occurrence and symmetrization
		int B = 2 * this->T_c + 1;
		int fullDim = (int)std::pow((float)B, this->order_c);

		int alreadyUsed;

		// MINMAX
		alreadyUsed = 0;
		MINMAXsymmCoord_c = new int[2 * fullDim]; // [0, fullDim-1] = min; [fullDim, 2*fullDim-1] = max
		for (int i = 0; i<2 * fullDim; i++) MINMAXsymmCoord_c[i] = -1;

		for (int numIter = 0; numIter < fullDim; numIter++)
		{
			if (MINMAXsymmCoord_c[numIter] == -1)
			{
				int coordReverse = 0;
				int num = numIter;
				for (int i = 0; i<this->order_c; i++)
				{
					coordReverse += (num % B) * ((int)std::pow((float)B, order_c - i - 1));
					num = num / B;
				}
				// To the same bin: min(X), max(-X), min(Xreverse), max(-Xreverse)
				if (MINMAXsymmCoord_c[numIter] == -1)
				{
					MINMAXsymmCoord_c[numIter] = alreadyUsed; // min(X)
					if (symmMinMax) MINMAXsymmCoord_c[2 * fullDim - numIter - 1] = alreadyUsed; // max(-X)
					if (symmReverse) MINMAXsymmCoord_c[coordReverse] = alreadyUsed; // min(Xreverse)
					if ((symmMinMax) && (symmReverse)) MINMAXsymmCoord_c[2 * fullDim - coordReverse - 1] = alreadyUsed; // max(-Xreverse)
					alreadyUsed++;
				}
			}
		}
		for (int numIter = 0; numIter < fullDim; numIter++)
		{
			if (MINMAXsymmCoord_c[fullDim + numIter] == -1)
			{
				int coordReverse = 0;
				int num = numIter;
				for (int i = 0; i<this->order_c; i++)
				{
					coordReverse += (num % B) * ((int)std::pow((float)B, order_c - i - 1));
					num = num / B;
				}
				// To the same bin: max(X), min(-X), max(Xreverse), min(-Xreverse)
				if (MINMAXsymmCoord_c[fullDim + numIter] == -1)
				{
					MINMAXsymmCoord_c[fullDim + numIter] = alreadyUsed; // max(X)
					if (symmMinMax) MINMAXsymmCoord_c[fullDim - numIter - 1] = alreadyUsed; // min(-X)
					if (symmReverse) MINMAXsymmCoord_c[fullDim + coordReverse] = alreadyUsed; // max(Xreverse)
					if ((symmMinMax) && (symmReverse)) MINMAXsymmCoord_c[fullDim - coordReverse - 1] = alreadyUsed; // min(-Xreverse)
					alreadyUsed++;
				}
			}
		}
		MINMAXsymmDim_c = alreadyUsed;

		// SPAM
		alreadyUsed = 0;
		SPAMsymmCoord_c = new int[fullDim];
		for (int i = 0; i<fullDim; i++) SPAMsymmCoord_c[i] = -1;
		for (int numIter = 0; numIter < fullDim; numIter++)
		{
			if (SPAMsymmCoord_c[numIter] == -1)
			{
				int coordReverse = 0;
				int num = numIter;
				for (int i = 0; i<this->order_c; i++)
				{
					coordReverse += (num % B) * ((int)std::pow((float)B, order_c - i - 1));
					num = num / B;
				}
				// To the same bin: X, -X, Xreverse, -Xreverse
				SPAMsymmCoord_c[numIter] = alreadyUsed; // X
				if (symmSign) SPAMsymmCoord_c[fullDim - numIter - 1] = alreadyUsed; // -X
				if (symmReverse) SPAMsymmCoord_c[coordReverse] = alreadyUsed; // Xreverse
				if ((symmSign) && (symmReverse)) SPAMsymmCoord_c[fullDim - coordReverse - 1] = alreadyUsed; // -Xreverse
				alreadyUsed++;
			}
		}
		SPAMsymmDim_c = alreadyUsed;
		// In order to have the same order of the features as the matlab SRM - shift +1
		for (int i = 0; i<fullDim; i++)
		{
			if (SPAMsymmCoord_c[i] == alreadyUsed - 1) SPAMsymmCoord_c[i] = 0;
			else SPAMsymmCoord_c[i]++;
		}
	}
};

#endif