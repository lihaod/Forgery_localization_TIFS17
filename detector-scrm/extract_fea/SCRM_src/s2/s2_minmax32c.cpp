#include "../submodel.h"
#include "../mat2D.h"
#include "../config.cpp"

class s2_minmax32c: public Submodel
{
public:
	s2_minmax32c(float q, Config *config) : Submodel(q) 
	{
		this->modelName =  "s2_minmax32c";
		this->minmax = true;

		this->coocDirs.push_back("col");
		this->coocDirs.push_back("col");

		Initialize(config, "color");
		this->normalizedFactor = 2;
	}

	~s2_minmax32c()
	{
	}

	virtual void ComputeFea(std::vector<mat2D<int> **> QResVect)
	{
		std::vector<std::vector<mat2D<int> **> > OpVect;

		// [0] - Horizontal, [1] - Vertical, [2] - Diagonal, [3] - Minor Diagonal

		// Horizontal + Vertical + Diagonal
		std::vector<mat2D<int> **> HVD = std::vector<mat2D<int> **>();
		HVD.push_back(QResVect[0]);HVD.push_back(QResVect[1]);HVD.push_back(QResVect[2]);
		OpVect.push_back(HVD);

		// Horizontal + Vertical + Minor Diagonal
		std::vector<mat2D<int> **> HVM = std::vector<mat2D<int> **>();
		HVM.push_back(QResVect[0]);HVM.push_back(QResVect[1]);HVM.push_back(QResVect[3]);
		OpVect.push_back(HVM);

		this->AddFea(OpVect);
	}
};
