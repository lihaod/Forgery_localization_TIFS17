#include "mat2D.h"
#include "submodel.h"
#include "s.h"
#include "config.cpp"
#include <vector>
#include <math.h>
#include <cmath>

s::s(std::vector<float> qs, Config *config)
{
	this->qs = qs;
	this->config = config;
}

s::~s()
{
	for (int i=0; i<(int)this->submodels.size(); i++) 
		for (int j=0; j<(int)this->submodels[i].size(); j++)
			delete submodels[i][j];
	delete[] endOfSubmodel_s;
	delete[] endOfSubmodel_c;
}

void s::GetResidual(mat2D<int> **img, mat2D<int>* kernel, mat2D<int> **residual)
{
	for (int channel = 0; channel < 3; channel++)
	{
		residual[channel] = new mat2D<int>(img[channel]->rows - kernel->rows + 1, img[channel]->cols - kernel->cols + 1);
		for (int ir = 0; ir < (img[channel]->rows - kernel->rows + 1); ir++)
		{
			for (int ic = 0; ic < (img[channel]->cols - kernel->cols + 1); ic++)
			{
				int convVal = 0;
				for (int kr = 0; kr < kernel->rows; kr++)
				{
					for (int kc = 0; kc < kernel->cols; kc++)
					{
						convVal = convVal + img[channel]->Read(ir + kr, ic + kc) * kernel->Read(kr, kc);
					}
				}
				residual[channel]->Write(ir, ic, convVal);
			}
		}
	}
}

void s::Quantize(mat2D<int>** residual, float totalKernelQ, mat2D<int> **quantizedResidual)
{
	for (int channel = 0; channel < 3; channel++)
	{
		quantizedResidual[channel] = new mat2D<int>(residual[channel]->rows, residual[channel]->cols);
		for (int r = 0; r < residual[channel]->rows; r++)
		{
			for (int c = 0; c<residual[channel]->cols; c++)
			{
				float tempValF = (float)residual[channel]->Read(r, c);
				tempValF = tempValF / totalKernelQ;
				int tempValI = 0;
				// rounding
				if (config->roundup5)
					tempValI = (int)((tempValF>0.0) ? floor(tempValF + 0.5) : ceil(tempValF - 0.5));
				else
					tempValI = (int)((tempValF - floor(tempValF) > 0.5) ? ceil(tempValF) : floor(tempValF));
				quantizedResidual[channel]->Write(r, c, tempValI);
			}
		}
	}
}

void s::Truncate(mat2D<int>** quantizedResidual, int T, mat2D<int> **QTResidual)
{
	for (int channel = 0; channel < 3; channel++)
	{
		QTResidual[channel] = new mat2D<int>(quantizedResidual[channel]->rows, quantizedResidual[channel]->cols);
		for (int r = 0; r < quantizedResidual[channel]->rows; r++)
		{
			for (int c = 0; c < quantizedResidual[channel]->cols; c++)
			{
				int tempVal = quantizedResidual[channel]->Read(r, c);
				tempVal = tempVal > T ? T : tempVal;
				tempVal = tempVal < -T ? -T : tempVal;
				QTResidual[channel]->Write(r, c, tempVal);
			}
		}
	}
}

void s::MultiplyByParity(std::vector<mat2D<int> **> QResVect, mat2D<int> **parity)
{
	for (int resIndex=0; resIndex < (int)QResVect.size(); resIndex++)
	{
		for (int channel = 0; channel < 3; channel++)
		{
			mat2D<int> *currentQRes = QResVect[resIndex][channel];

			for (int r = 0; r < currentQRes->rows; r++)
				for (int c = 0; c < currentQRes->cols; c++)
				{
					int iValue = currentQRes->Read(r, c) * parity[channel]->Read(r + this->cutEdgesForParityBy, c + this->cutEdgesForParityBy);
					currentQRes->Write(r, c, iValue);
				}
		}
	}
}
