#ifndef S_H_
#define S_H_

#include "mat2D.h"
#include "submodel.h"
#include "config.cpp"
#include <vector>

class s
{
public:
	std::vector<std::vector<Submodel *> > submodels;

	s(std::vector<float> qs, Config *config);
	~s();

	virtual void ComputeImage(mat2D<int> **img, mat2D<double> *map, mat2D<int> **parity) = 0;

protected:
	Config *config;
	std::vector<float> qs;
	int quantMultiplier;
	int cutEdgesForParityBy;
	int *endOfSubmodel_s;
	int *endOfSubmodel_c;

	void GetResidual(mat2D<int> **img, mat2D<int>* kernel, mat2D<int> **residual);
	void Quantize(mat2D<int>** residual, float totalKernelQ, mat2D<int> **quantizedResidual);
	void Truncate(mat2D<int>** quantizedResidual, int T, mat2D<int> **QTResidual);
	void MultiplyByParity(std::vector<mat2D<int> **> QResVect, mat2D<int> **parity);

};

#endif
