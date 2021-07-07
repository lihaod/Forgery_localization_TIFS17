#include "../submodel.h"
#include "../mat2D.h"
#include "../config.cpp"

class s3x3_minmax22v: public Submodel
{
public:
	s3x3_minmax22v(float q, Config *config) : Submodel(q) 
	{
		this->modelName =  "s3x3_minmax22v";
		this->minmax = true;

		this->coocDirs.push_back("ver");
		this->coocDirs.push_back("hor");

		Initialize(config);
	}

	~s3x3_minmax22v()
	{
	}

	virtual void ComputeFea(std::vector<mat2D<int> **> QResVect)
	{
		std::vector<std::vector<mat2D<int> **> > OpVect;

		// [0] - Right, [1] - Left, [2] - Up, [3] - Down
		// [4] - All

		// Up + Down
		std::vector<mat2D<int> **> UD = std::vector<mat2D<int> **>();
		UD.push_back(QResVect[2]);UD.push_back(QResVect[3]);
		OpVect.push_back(UD);

		// Right + Left
		std::vector<mat2D<int> **> RL = std::vector<mat2D<int> **>();
		RL.push_back(QResVect[0]);RL.push_back(QResVect[1]);
		OpVect.push_back(RL);

		this->AddFea(OpVect);
	}
};
