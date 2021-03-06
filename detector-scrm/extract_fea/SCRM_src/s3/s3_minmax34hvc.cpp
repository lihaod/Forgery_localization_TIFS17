#include "../submodel.h"
#include "../mat2D.h"
#include "../config.cpp"

class s3_minmax34hvc: public Submodel
{
public:
	s3_minmax34hvc(float q, Config *config) : Submodel(q) 
	{
		this->modelName =  "s3_minmax34hvc";
		this->minmax = true;

		this->coocDirs.push_back("col");
		this->coocDirs.push_back("col");
		this->coocDirs.push_back("col");
		this->coocDirs.push_back("col");

		Initialize(config, "color");
	}

	~s3_minmax34hvc()
	{
	}

	virtual void ComputeFea(std::vector<mat2D<int> **> QResVect)
	{
		std::vector<std::vector<mat2D<int> **> > OpVect;

		// [0] - Right, [1] - Left, [2] - Up, [3] - Down
		// [4] - Right_Up, [5] - Right_Down, [6] - Left_Up, [7] - Left_Down

		// Right Up Left
		std::vector<mat2D<int> **> RUL = std::vector<mat2D<int> **>();
		RUL.push_back(QResVect[0]);RUL.push_back(QResVect[2]);RUL.push_back(QResVect[1]);
		OpVect.push_back(RUL);

		// Right Down Left
		std::vector<mat2D<int> **> RDL = std::vector<mat2D<int> **>();
		RDL.push_back(QResVect[0]);RDL.push_back(QResVect[3]);RDL.push_back(QResVect[1]);
		OpVect.push_back(RDL);

		// Up Right Down
		std::vector<mat2D<int> **> URD = std::vector<mat2D<int> **>();
		URD.push_back(QResVect[2]);URD.push_back(QResVect[0]);URD.push_back(QResVect[3]);
		OpVect.push_back(URD);

		// Up Left Down
		std::vector<mat2D<int> **> ULD = std::vector<mat2D<int> **>();
		ULD.push_back(QResVect[2]);ULD.push_back(QResVect[1]);ULD.push_back(QResVect[3]);
		OpVect.push_back(ULD);

		this->AddFea(OpVect);
	}
};
