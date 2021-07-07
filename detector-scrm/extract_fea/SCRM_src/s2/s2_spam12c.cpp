#include "../submodel.h"
#include "../mat2D.h"
#include "../config.cpp"

class s2_spam12c: public Submodel
{
public:
	s2_spam12c(float q, Config *config) : Submodel(q) 
	{
		this->modelName =  "s2_spam12c";
		this->minmax = false;
		this->coocDirs.push_back("col");
		this->coocDirs.push_back("col");

		Initialize(config, "color");
	}

	~s2_spam12c()
	{
	}

	virtual void ComputeFea(std::vector<mat2D<int> **> QResVect)
	{
		std::vector<std::vector<mat2D<int> **> > OpVect;

		// [0] - Horizontal, [1] - Vertical, [2] - Diagonal, [3] - Minor Diagonal

		// Horizontal
		std::vector<mat2D<int> **> H = std::vector<mat2D<int> **>();
		H.push_back(QResVect[0]);
		OpVect.push_back(H);

		// Vertical
		std::vector<mat2D<int> **> V = std::vector<mat2D<int> **>();
		V.push_back(QResVect[1]);
		OpVect.push_back(V);

		this->AddFea(OpVect);
	}
};
