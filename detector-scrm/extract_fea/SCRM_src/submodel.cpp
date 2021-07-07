#include "submodel.h"
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include "mat2D.h"
#include <math.h>
#include "config.cpp"
#include <limits>
#include <string.h>

Submodel::Submodel(float q)
{
	this->q = q;
	mergeInto = "";
}

Submodel::~Submodel()
{
	for (int i=0; i<(int)this->fea.size(); i++) delete [] this->fea[i];
	delete [] this->multi;
}

std::vector<float *> Submodel::ReturnFea()
{
	return this->fea;
}

void Submodel::AddFea(std::vector<std::vector<mat2D<int>** > > OpVect)
{	
	float *cooc = new float[this->symmDim*this->numInMap]; 
	for (int i = 0; i<this->symmDim*this->numInMap; i++) cooc[i] = 0;
	// if MINMAX-like submodel
	if (this->minmax)
	{
		for (int OpIndex=0; OpIndex < (int)OpVect.size(); OpIndex++)
		{
			mat2D<int> *min[3] = { new mat2D<int>(OpVect[OpIndex][0][0]->rows, OpVect[OpIndex][0][0]->cols),
				new mat2D<int>(OpVect[OpIndex][0][1]->rows, OpVect[OpIndex][0][1]->cols),
				new mat2D<int>(OpVect[OpIndex][0][2]->rows, OpVect[OpIndex][0][2]->cols) };
			mat2D<int> *max[3] = { new mat2D<int>(OpVect[OpIndex][0][0]->rows, OpVect[OpIndex][0][0]->cols),
				new mat2D<int>(OpVect[OpIndex][0][1]->rows, OpVect[OpIndex][0][1]->cols),
				new mat2D<int>(OpVect[OpIndex][0][2]->rows, OpVect[OpIndex][0][2]->cols) };

			// Merge all the residuals into MIN and MAX residual
			this->GetMinMax(OpVect[OpIndex], min, max);
			
			// Get the co-occurrence
			float *tempCooc = this->GetCooc(min, max, this->coocDirs[OpIndex]);
			delete min[0]; delete min[1]; delete min[2];
			delete max[0]; delete max[1]; delete max[2];

			// Add the co-occurrences to the previous ones
			for (int i = 0; i<this->symmDim*this->numInMap; i++) cooc[i] += tempCooc[i];

			delete [] tempCooc;
		}
	}
	// if SPAM-like submodel
	else
	{
		for (int OpIndex=0; OpIndex < (int)OpVect.size(); OpIndex++)
		{
			// Get the co-occurrence
			float *tempCooc = this->GetCooc(OpVect[OpIndex][0], NULL, this->coocDirs[OpIndex]);

			// Add the co-occurrences to the previous ones
			for (int i = 0; i<this->symmDim*this->numInMap; i++) cooc[i] += tempCooc[i];

			delete [] tempCooc;
		}
	}
	// Normalize coocs to sum(Cooc)==1
	float correctionFactor = 3;
	if (feaTypeIsColor) correctionFactor = 1;
	for (int j = 0; j < this->numInMap; j++)
	{
		int totalCoocs = 0;
		int offset = j*this->symmDim;
		for (int i=0; i<this->symmDim; i++)
			totalCoocs += (int)cooc[i + offset];
		if (this->minmax) for (int i = 0; i<this->symmDim; i++) cooc[i + offset] /= (totalCoocs / this->normalizedFactor / correctionFactor);
		else for (int i = 0; i<this->symmDim; i++) cooc[i + offset] /= (totalCoocs / this->normalizedFactor / correctionFactor);
	}
	// Add the image features to the FEA vector
	fea.push_back(cooc);
}

void Submodel::GetMinMax(std::vector<mat2D<int> **> residuals, mat2D<int> **outMin, mat2D<int> **outMax)
{
	for (int channel = 0; channel < 3; channel++)
	{
		for (int r = 0; r < outMin[channel]->rows; r++)
		{
			for (int c = 0; c < outMin[channel]->cols; c++)
			{
				int min = std::numeric_limits<int>::max();
				int max = std::numeric_limits<int>::min();
				for (int i = 0; i < (int)residuals.size(); i++)
				{
					int resVal = residuals[i][channel]->Read(r, c);
					if (resVal < min) min = resVal;
					if (resVal > max) max = resVal;
				}
				outMin[channel]->Write(r, c, min);
				outMax[channel]->Write(r, c, max);
			}
		}
	}
}

float* Submodel::GetCooc(mat2D<int>** SPAM_Min_Residual, mat2D<int>** Max_Residual, std::string coocDir)
{
	float *cooc = new float[this->symmDim*this->numInMap];
	for (int i = 0; i<this->symmDim*this->numInMap; i++) cooc[i] = 0;

	// Horizontal co-occurrence
	if ((coocDir.compare("hor") == 0) || (coocDir.compare("horver") == 0))
	{
		for (int channel = 0; channel < 3; channel++)
		{
			for (int r = 0; r < SPAM_Min_Residual[channel]->rows; r++)
			{
				for (int c = 0; c < SPAM_Min_Residual[channel]->cols - this->order + 1; c++)
				{
					// MINMAX features
					if (this->minmax)
					{
						// MIN
						int feaNumber = 0;
						for (int o = 0; o < order; o++)
							feaNumber += (SPAM_Min_Residual[channel]->Read(r, c + o) + this->T)*this->multi[o];

						int dimToIncrease = this->inSymmCoord[feaNumber];
						if (this->mapMode)
							cooc[dimToIncrease] += map->Read(r, c);
						else
							cooc[dimToIncrease + ((int)map->Read(r, c))*this->symmDim] ++;

						// MAX
						feaNumber = 0;
						for (int o = 0; o < order; o++)
							feaNumber += (Max_Residual[channel]->Read(r, c + o) + this->T)*this->multi[o];

						feaNumber += fullDim;
						dimToIncrease = this->inSymmCoord[feaNumber];
						if (this->mapMode)
							cooc[dimToIncrease] += map->Read(r, c);
						else
							cooc[dimToIncrease + ((int)map->Read(r, c))*this->symmDim] ++;
					}
					// SPAM-like features
					else
					{
						int feaNumber = 0;
						for (int o = 0; o < order; o++)
							feaNumber += (SPAM_Min_Residual[channel]->Read(r, c + o) + this->T)*this->multi[o];
						int dimToIncrease = this->inSymmCoord[feaNumber];
						if (this->mapMode)
							cooc[dimToIncrease] += map->Read(r, c);
						else
							cooc[dimToIncrease + ((int)map->Read(r, c))*this->symmDim] ++;
					}
				}
			}
		}
	}
	// Vertical co-occurrence
	if ((coocDir.compare("ver") == 0) || (coocDir.compare("horver") == 0))
	{
		for (int channel = 0; channel < 3; channel++)
		{
			for (int r = 0; r < SPAM_Min_Residual[channel]->rows - this->order + 1; r++)
			{
				for (int c = 0; c < SPAM_Min_Residual[channel]->cols; c++)
				{
					// MINMAX features
					if (this->minmax)
					{
						// MIN
						int feaNumber = 0;
						for (int o = 0; o < order; o++)
							feaNumber += (SPAM_Min_Residual[channel]->Read(r + o, c) + this->T)*this->multi[o];

						int dimToIncrease = this->inSymmCoord[feaNumber];
						if (this->mapMode)
							cooc[dimToIncrease] += map->Read(r, c);
						else
							cooc[dimToIncrease + ((int)map->Read(r, c))*this->symmDim] ++;

						// MAX
						feaNumber = 0;
						for (int o = 0; o < order; o++)
							feaNumber += (Max_Residual[channel]->Read(r + o, c) + this->T)*this->multi[o];

						feaNumber += fullDim;
						dimToIncrease = this->inSymmCoord[feaNumber];
						if (this->mapMode)
							cooc[dimToIncrease] += map->Read(r, c);
						else
							cooc[dimToIncrease + ((int)map->Read(r, c))*this->symmDim] ++;

					}
					// SPAM-like features
					else
					{
						int feaNumber = 0;
						for (int o = 0; o < order; o++)
							feaNumber += (SPAM_Min_Residual[channel]->Read(r + o, c) + this->T)*this->multi[o];
						int dimToIncrease = this->inSymmCoord[feaNumber];
						if (this->mapMode)
							cooc[dimToIncrease] += map->Read(r, c);
						else
							cooc[dimToIncrease + ((int)map->Read(r, c))*this->symmDim] ++;
					}
				}
			}
		}
	}
	// Inter channel co-occurrence
	if (coocDir.compare("col")==0)
	{
		for (int r=0; r<SPAM_Min_Residual[0]->rows; r++)
		{
			for (int c=0; c < SPAM_Min_Residual[0]->cols; c++)
			{
				// MINMAX features
				if (this->minmax)
				{
					// MIN
					int feaNumber = 0;
					for (int o=0; o<order;o++)
						feaNumber += (SPAM_Min_Residual[o]->Read(r, c) + this->T)*this->multi[o];
	
					int dimToIncrease = this->inSymmCoord[feaNumber];
    				if (this->mapMode)
						cooc[dimToIncrease] += map->Read(r, c);
					else
						cooc[dimToIncrease + ((int)map->Read(r, c))*this->symmDim] ++;

					// MAX
					feaNumber = 0;
					for (int o=0; o<order;o++)
						feaNumber += (Max_Residual[o]->Read(r, c) + this->T)*this->multi[o];
	
					feaNumber += fullDim;
					dimToIncrease = this->inSymmCoord[feaNumber];
    				if (this->mapMode)
						cooc[dimToIncrease] += map->Read(r, c);
					else
						cooc[dimToIncrease + ((int)map->Read(r, c))*this->symmDim] ++;
				}
				// SPAM-like features
				else
				{
					int feaNumber = 0;
					for (int o=0; o<order;o++)
						feaNumber += (SPAM_Min_Residual[o]->Read(r, c) + this->T)*this->multi[o];
					int dimToIncrease = this->inSymmCoord[feaNumber];
    				if (this->mapMode)
						cooc[dimToIncrease] += map->Read(r, c);
					else
						cooc[dimToIncrease + ((int)map->Read(r, c))*this->symmDim] ++;
				}
			}	
		}
	}

	return cooc;
}

std::string Submodel::GetName()
{
    std::stringstream qs (std::stringstream::in | std::stringstream::out);
    qs << q;
	std::string qstring = qs.str();
	int i;
	for (i=0; (i<(int)qstring.length()) && (qstring[i]!='.'); i++);
	if (i < (int)qstring.length()) qstring.erase(i, 1);
	return this->modelName + "_q" + qstring;
}

void Submodel::Initialize(Config *config, const char *feaType)
{
	if (strcmp("color", feaType)==0)
	{
		this->order = config->order_c;
		this->T = config->T_c;
		this->feaTypeIsColor = true;
	}
	else
	{
		this->order = config->order_s;
		this->T = config->T_s;
		this->feaTypeIsColor = false;
	}
	this->numInMap = config->numInMap;
	this->mapMode = config->mapMode;

	this->multi = new int[this->order];
	for (int multiIndex=0; multiIndex < this->order; multiIndex++)
		multi[multiIndex] = (int)std::pow((float)2*this->T+1, multiIndex);
	int B = 2*this->T+1;
	fullDim = (int)std::pow((float)B, this->order);
	if (this->feaTypeIsColor)
	{
		if (minmax)
		{
			this->inSymmCoord = config->MINMAXsymmCoord_c;
			this->symmDim = config->MINMAXsymmDim_c;
			this->normalizedFactor = 4;
		}
		else
		{
			this->inSymmCoord = config->SPAMsymmCoord_c;
			this->symmDim = config->SPAMsymmDim_c;
			this->normalizedFactor = 2;
		}
	}
	else
	{
		if (minmax)
		{
			this->inSymmCoord = config->MINMAXsymmCoord_s;
			this->symmDim = config->MINMAXsymmDim_s;
			this->normalizedFactor = 4;
		}
		else
		{
			this->inSymmCoord = config->SPAMsymmCoord_s;
			this->symmDim = config->SPAMsymmDim_s;
			this->normalizedFactor = 2;
		}
	}
}
