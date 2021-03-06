#include "../submodel.h"
#include "../mat2D.h"
#include "../config.cpp"

class s1_minmax41c: public Submodel
{
public:
	s1_minmax41c(float q, Config *config) : Submodel(q) 
	{
		this->modelName =  "s1_minmax41c";
		this->minmax = true;

		this->coocDirs.push_back("col");

		Initialize(config, "color");
		this->normalizedFactor = 2;
	}

	~s1_minmax41c()
	{
	}

	virtual void ComputeFea(std::vector<mat2D<int> **> QResVect)
	{
		std::vector<std::vector<mat2D<int> **> > OpVect;

		// [0] - Right, [1] - Left, [2] - Up, [3] - Down
		// [4] - Right_Up, [5] - Right_Down, [6] - Left_Up, [7] - Left_Down

		// Twice the same, vertical and horizontal co-occurrence
		// Right Up Left Down
		std::vector<mat2D<int> **> RLUD = std::vector<mat2D<int> **>();
		RLUD.push_back(QResVect[0]);RLUD.push_back(QResVect[1]);RLUD.push_back(QResVect[2]);RLUD.push_back(QResVect[3]);
		OpVect.push_back(RLUD);

		this->AddFea(OpVect);
	}
};
