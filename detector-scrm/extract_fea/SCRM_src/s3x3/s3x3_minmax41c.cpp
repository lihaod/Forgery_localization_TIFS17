#include "../submodel.h"
#include "../mat2D.h"
#include "../config.cpp"

class s3x3_minmax41c: public Submodel
{
public:
	s3x3_minmax41c(float q, Config *config) : Submodel(q) 
	{
		this->modelName =  "s3x3_minmax41c";
		this->minmax = true;

		this->coocDirs.push_back("col");

		Initialize(config, "color");
		this->normalizedFactor = 2;
	}

	~s3x3_minmax41c()
	{
	}

	virtual void ComputeFea(std::vector<mat2D<int> **> QResVect)
	{
		std::vector<std::vector<mat2D<int> **> > OpVect;

		// [0] - Right, [1] - Left, [2] - Up, [3] - Down
		// [4] - All

		// Twice the same, vertical and horizontal co-occurrence
		// Right Up Left Down
		std::vector<mat2D<int> **> RLUD = std::vector<mat2D<int> **>();
		RLUD.push_back(QResVect[0]);RLUD.push_back(QResVect[1]);RLUD.push_back(QResVect[2]);RLUD.push_back(QResVect[3]);
		OpVect.push_back(RLUD);

		this->AddFea(OpVect);
	}
};
