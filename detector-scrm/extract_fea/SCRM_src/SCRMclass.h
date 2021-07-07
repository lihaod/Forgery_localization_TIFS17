#ifndef CRMCLASS_H_
#define CRMCLASS_H_


#include "submodel.h"
#include "mat2D.h"
#include <vector>
#include "s.h"
#include "config.cpp"

class SCRMclass
{
public:
	std::vector<s *> submodelClasses;

	SCRMclass(Config *config);
	~SCRMclass();
	void ComputeFeatures(void);

	void ComputeFeatures(mat2D<int> **I, mat2D<double> * map);
	std::vector<Submodel *> GetSubmodels();

private:
	// bool verbose;
	Config *config;
	std::vector<Submodel *> AddedMergesSpams;

	std::vector<Submodel *> PostProcessing(std::vector<Submodel *> submodels);
	mat2D<int> *GetParity(mat2D<int> *I);
	void clearAddedMergesSpams();
};

#endif