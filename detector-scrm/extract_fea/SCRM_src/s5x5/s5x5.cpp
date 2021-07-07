/*
 This function outputs co-occurrences of ALL 3rd-order residuals
 listed in Figure 1 in our journal HUGO paper (version from June 14), 
 including the naming convention.
 */

#include "../mat2D.h"
#include "../submodel.h"
#include "../config.cpp"
#include "../s.h"

#include "s5x5_spam11.cpp"
#include "s5x5_spam14h.cpp"
#include "s5x5_spam14v.cpp"
#include "s5x5_minmax22h.cpp"
#include "s5x5_minmax22v.cpp"
#include "s5x5_minmax24.cpp"
#include "s5x5_minmax41.cpp"
#include "s5x5_spam11c.cpp"
#include "s5x5_spam14c.cpp"
#include "s5x5_minmax22c.cpp"
#include "s5x5_minmax24c.cpp"
#include "s5x5_minmax41c.cpp"

class s5x5 : s
{
public:
	void CreateKernels()
	{
		mat2D<int> *temp;
		cutEdgesForParityBy = 2;

		// Right Kernel
		temp = new mat2D<int>(5, 5);
		temp->Write(0, 0, 0); temp->Write(0, 1, 0); temp->Write(0, 2, -2); temp->Write(0, 3, 2); temp->Write(0, 4,-1);
		temp->Write(1, 0, 0); temp->Write(1, 1, 0); temp->Write(1, 2,  8); temp->Write(1, 3,-6); temp->Write(1, 4, 2);
		temp->Write(2, 0, 0); temp->Write(2, 1, 0); temp->Write(2, 2,-12); temp->Write(2, 3, 8); temp->Write(2, 4,-2);
		temp->Write(3, 0, 0); temp->Write(3, 1, 0); temp->Write(3, 2,  8); temp->Write(3, 3,-6); temp->Write(3, 4, 2);
		temp->Write(4, 0, 0); temp->Write(4, 1, 0); temp->Write(4, 2, -2); temp->Write(4, 3, 2); temp->Write(4, 4,-1);
		kerR = temp;

		// Left Kernel
		temp = new mat2D<int>(5, 5);
		temp->Write(0, 0,-1); temp->Write(0, 1, 2); temp->Write(0, 2, -2); temp->Write(0, 3, 0); temp->Write(0, 4, 0);
		temp->Write(1, 0, 2); temp->Write(1, 1,-6); temp->Write(1, 2,  8); temp->Write(1, 3, 0); temp->Write(1, 4, 0);
		temp->Write(2, 0,-2); temp->Write(2, 1, 8); temp->Write(2, 2,-12); temp->Write(2, 3, 0); temp->Write(2, 4, 0);
		temp->Write(3, 0, 2); temp->Write(3, 1,-6); temp->Write(3, 2,  8); temp->Write(3, 3, 0); temp->Write(3, 4, 0);
		temp->Write(4, 0,-1); temp->Write(4, 1, 2); temp->Write(4, 2, -2); temp->Write(4, 3, 0); temp->Write(4, 4, 0);
		kerL = temp;

		// Up Kernel
		temp = new mat2D<int>(5, 5);
		temp->Write(0, 0,-1); temp->Write(0, 1, 2); temp->Write(0, 2, -2); temp->Write(0, 3, 2); temp->Write(0, 4,-1);
		temp->Write(1, 0, 2); temp->Write(1, 1,-6); temp->Write(1, 2,  8); temp->Write(1, 3,-6); temp->Write(1, 4, 2);
		temp->Write(2, 0,-2); temp->Write(2, 1, 8); temp->Write(2, 2,-12); temp->Write(2, 3, 8); temp->Write(2, 4,-2);
		temp->Write(3, 0, 0); temp->Write(3, 1, 0); temp->Write(3, 2,  0); temp->Write(3, 3, 0); temp->Write(3, 4, 0);
		temp->Write(4, 0, 0); temp->Write(4, 1, 0); temp->Write(4, 2,  0); temp->Write(4, 3, 0); temp->Write(4, 4, 0);
		kerU = temp;

		// Down Kernel
		temp = new mat2D<int>(5, 5);
		temp->Write(0, 0, 0); temp->Write(0, 1, 0); temp->Write(0, 2,  0); temp->Write(0, 3, 0); temp->Write(0, 4, 0);
		temp->Write(1, 0, 0); temp->Write(1, 1, 0); temp->Write(1, 2,  0); temp->Write(1, 3, 0); temp->Write(1, 4, 0);
		temp->Write(2, 0,-2); temp->Write(2, 1, 8); temp->Write(2, 2,-12); temp->Write(2, 3, 8); temp->Write(2, 4,-2);
		temp->Write(3, 0, 2); temp->Write(3, 1,-6); temp->Write(3, 2,  8); temp->Write(3, 3,-6); temp->Write(3, 4, 2);
		temp->Write(4, 0,-1); temp->Write(4, 1, 2); temp->Write(4, 2, -2); temp->Write(4, 3, 2); temp->Write(4, 4,-1);
		kerD = temp;

		// All Kernel
		temp = new mat2D<int>(5, 5);
		temp->Write(0, 0,-1); temp->Write(0, 1, 2); temp->Write(0, 2, -2); temp->Write(0, 3, 2); temp->Write(0, 4,-1);
		temp->Write(1, 0, 2); temp->Write(1, 1,-6); temp->Write(1, 2,  8); temp->Write(1, 3,-6); temp->Write(1, 4, 2);
		temp->Write(2, 0,-2); temp->Write(2, 1, 8); temp->Write(2, 2,-12); temp->Write(2, 3, 8); temp->Write(2, 4,-2);
		temp->Write(3, 0, 2); temp->Write(3, 1,-6); temp->Write(3, 2,  8); temp->Write(3, 3,-6); temp->Write(3, 4, 2);
		temp->Write(4, 0,-1); temp->Write(4, 1, 2); temp->Write(4, 2, -2); temp->Write(4, 3, 2); temp->Write(4, 4,-1);
		kerAll = temp;
	}

	s5x5(std::vector<float> qs, Config *config) : s(qs, config)
	{
		this->CreateKernels();
		quantMultiplier = 12;
		endOfSubmodel_s = new int[(int)qs.size()];
		endOfSubmodel_c = new int[(int)qs.size()];

		for (int qIndex=0; qIndex < (int)qs.size(); qIndex++)
		{
			float q = qs[qIndex];
			std::vector<Submodel *> submodelsForQ;

			endOfSubmodel_c[qIndex] = 0;
			if (config->enable_c)
			{
				submodelsForQ.push_back(new s5x5_spam11c(q, config));
				submodelsForQ.push_back(new s5x5_spam14c(q, config));
				submodelsForQ.push_back(new s5x5_minmax24c(q, config));
				submodelsForQ.push_back(new s5x5_minmax22c(q, config));
				submodelsForQ.push_back(new s5x5_minmax41c(q, config));
				endOfSubmodel_c[qIndex] = (int)submodelsForQ.size();
			}

			endOfSubmodel_s[qIndex] = endOfSubmodel_c[qIndex];
			if (config->enable_s)
			{
				submodelsForQ.push_back(new s5x5_minmax24(q, config));
				submodelsForQ.push_back(new s5x5_minmax22h(q, config));
				submodelsForQ.push_back(new s5x5_minmax22v(q, config));
				submodelsForQ.push_back(new s5x5_minmax41(q, config));
				submodelsForQ.push_back(new s5x5_spam14h(q, config));
				submodelsForQ.push_back(new s5x5_spam14v(q, config));
				submodelsForQ.push_back(new s5x5_spam11(q, config));
				endOfSubmodel_s[qIndex] = (int)submodelsForQ.size();
			}

			this->submodels.push_back(submodelsForQ);
		}
	}

	~s5x5()
	{
		delete kerR; delete kerL; delete kerU; delete kerD;
		delete kerAll;
	}

	void ComputeImage(mat2D<int> **img, mat2D<double> * map, mat2D<int> **parity)
	{
		mat2D<int> *R[3] = { NULL, NULL, NULL }; GetResidual(img, kerR, R);
		mat2D<int> *L[3] = { NULL, NULL, NULL }; GetResidual(img, kerL, L);
		mat2D<int> *U[3] = { NULL, NULL, NULL }; GetResidual(img, kerU, U);
		mat2D<int> *D[3] = { NULL, NULL, NULL }; GetResidual(img, kerD, D);
		mat2D<int> *All[3] = { NULL, NULL, NULL }; GetResidual(img, kerAll, All);

		mat2D<double> * pMap = new mat2D<double>(img[0]->rows - 4, img[0]->cols - 4);
		for (int i = 0; i<img[0]->rows - 4; i++)
			for (int j = 0; j<img[0]->cols - 4; j++)
				pMap->Write(i, j, map->Read(i + 2, j + 2));
		for (int i = 0; i<this->submodels.size(); i++)
			for (int j = 0; j<this->submodels[i].size(); j++)
				this->submodels[i][j]->map = pMap;

		for (int qIndex = 0; qIndex < (int)submodels.size(); qIndex++)
		{
			float q = qs[qIndex] * quantMultiplier;

			mat2D<int> *QRes_R[3] = { NULL, NULL, NULL }; Quantize(R, q, QRes_R);
			mat2D<int> *QRes_L[3] = { NULL, NULL, NULL }; Quantize(L, q, QRes_L);
			mat2D<int> *QRes_U[3] = { NULL, NULL, NULL }; Quantize(U, q, QRes_U);
			mat2D<int> *QRes_D[3] = { NULL, NULL, NULL }; Quantize(D, q, QRes_D);
			mat2D<int> *QRes_All[3] = { NULL, NULL, NULL }; Quantize(All, q, QRes_All);

			int submodelIndex = 0;
			std::vector<mat2D<int> **> QTResVect;
			mat2D<int> *QTRes_R[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_R);
			mat2D<int> *QTRes_L[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_L);
			mat2D<int> *QTRes_U[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_U);
			mat2D<int> *QTRes_D[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_D);
			mat2D<int> *QTRes_All[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_All);

			if (config->enable_c)
			{
				Truncate(QRes_R, config->T_c, QTRes_R);
				Truncate(QRes_L, config->T_c, QTRes_L);
				Truncate(QRes_U, config->T_c, QTRes_U);
				Truncate(QRes_D, config->T_c, QTRes_D);
				Truncate(QRes_All, config->T_c, QTRes_All);
				if (config->parity) MultiplyByParity(QTResVect, parity);
				while (submodelIndex < endOfSubmodel_c[qIndex])
				{
					submodels[qIndex][submodelIndex]->ComputeFea(QTResVect);
					submodelIndex++;
				}

				for (int i = 0; i < (int)QTResVect.size(); i++)
					for (int channel = 0; channel < 3; channel++)
						delete QTResVect[i][channel];
			}

			if (config->enable_s)
			{
				Truncate(QRes_R, config->T_s, QTRes_R);
				Truncate(QRes_L, config->T_s, QTRes_L);
				Truncate(QRes_U, config->T_s, QTRes_U);
				Truncate(QRes_D, config->T_s, QTRes_D);
				Truncate(QRes_All, config->T_s, QTRes_All);
				if (config->parity) MultiplyByParity(QTResVect, parity);
				while (submodelIndex < endOfSubmodel_s[qIndex])
				{
					submodels[qIndex][submodelIndex]->ComputeFea(QTResVect);
					submodelIndex++;
				}

				for (int i = 0; i < (int)QTResVect.size(); i++)
					for (int channel = 0; channel < 3; channel++)
						delete QTResVect[i][channel];
			}

			for (int channel = 0; channel < 3; channel++)
			{
				delete QRes_R[channel]; delete QRes_L[channel]; delete QRes_U[channel]; delete QRes_D[channel];
				delete QRes_All[channel];
			}
		}

		for (int channel = 0; channel < 3; channel++)
		{
			delete R[channel]; delete L[channel]; delete U[channel]; delete D[channel];
			delete All[channel];
		}
		delete pMap;
	}

private:
	mat2D<int> *kerR;
	mat2D<int> *kerL;
	mat2D<int> *kerU;
	mat2D<int> *kerD;
	mat2D<int> *kerAll;
};