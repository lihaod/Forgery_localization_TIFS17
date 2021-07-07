/*
 This function outputs co-occurrences of ALL 3rd-order residuals
 listed in Figure 1 in our journal HUGO paper (version from June 14), 
 including the naming convention.
 List of outputted features:

 1a) spam14h
 1b) spam14v (orthogonal-spam)
 1c) minmax22v
 1d) minmax24
 1e) minmax34v
 1f) minmax41
 1g) minmax34
 1h) minmax48h
 1i) minmax54

 Naming convention:

 name = {type}{f}{sigma}{scan}
 type \in {spam, minmax}
 f \in {1,2,3,4,5} number of filters that are "minmaxed"
 sigma \in {1,2,3,4,8} symmetry index
 scan \in {h,v,\emptyset} scan of the cooc matrix (empty = sum of both 
 h and v scans).
 */

#include "../mat2D.h"
#include "../submodel.h"
#include "../config.cpp"
#include "../s.h"

#include "s3_spam14h.cpp"
#include "s3_spam14v.cpp"
#include "s3_minmax22h.cpp"
#include "s3_minmax22v.cpp"
#include "s3_minmax24.cpp"
#include "s3_minmax34.cpp"
#include "s3_minmax34h.cpp"
#include "s3_minmax34v.cpp"
#include "s3_minmax41.cpp"
#include "s3_minmax48h.cpp"
#include "s3_minmax48v.cpp"
#include "s3_minmax54.cpp"
#include "s3_spam14c.cpp"
#include "s3_minmax22c.cpp"
#include "s3_minmax24c.cpp"
#include "s3_minmax34c.cpp"
#include "s3_minmax34hvc.cpp"
#include "s3_minmax41c.cpp"
#include "s3_minmax48c.cpp"
#include "s3_minmax54c.cpp"

class s3 : s
{
public:
	void CreateKernels()
	{
		mat2D<int> *temp;
		cutEdgesForParityBy = 2;
		
		// Right Kernel
		temp = new mat2D<int>(5, 5);
		temp->Write(0, 0, 0); temp->Write(0, 1, 0); temp->Write(0, 2, 0); temp->Write(0, 3, 0); temp->Write(0, 4, 0);
		temp->Write(1, 0, 0); temp->Write(1, 1, 0); temp->Write(1, 2, 0); temp->Write(1, 3, 0); temp->Write(1, 4, 0);
		temp->Write(2, 0, 0); temp->Write(2, 1, 1); temp->Write(2, 2,-3); temp->Write(2, 3, 3); temp->Write(2, 4,-1);
		temp->Write(3, 0, 0); temp->Write(3, 1, 0); temp->Write(3, 2, 0); temp->Write(3, 3, 0); temp->Write(3, 4, 0);
		temp->Write(4, 0, 0); temp->Write(4, 1, 0); temp->Write(4, 2, 0); temp->Write(4, 3, 0); temp->Write(4, 4, 0);
		kerR = temp;

		// Left Kernel
		temp = new mat2D<int>(5, 5);
		temp->Write(0, 0, 0); temp->Write(0, 1, 0); temp->Write(0, 2, 0); temp->Write(0, 3, 0); temp->Write(0, 4, 0);
		temp->Write(1, 0, 0); temp->Write(1, 1, 0); temp->Write(1, 2, 0); temp->Write(1, 3, 0); temp->Write(1, 4, 0);
		temp->Write(2, 0,-1); temp->Write(2, 1, 3); temp->Write(2, 2,-3); temp->Write(2, 3, 1); temp->Write(2, 4, 0);
		temp->Write(3, 0, 0); temp->Write(3, 1, 0); temp->Write(3, 2, 0); temp->Write(3, 3, 0); temp->Write(3, 4, 0);
		temp->Write(4, 0, 0); temp->Write(4, 1, 0); temp->Write(4, 2, 0); temp->Write(4, 3, 0); temp->Write(4, 4, 0);
		kerL = temp;

		// Up Kernel
		temp = new mat2D<int>(5, 5);
		temp->Write(0, 0, 0); temp->Write(0, 1, 0); temp->Write(0, 2,-1); temp->Write(0, 3, 0); temp->Write(0, 4, 0);
		temp->Write(1, 0, 0); temp->Write(1, 1, 0); temp->Write(1, 2, 3); temp->Write(1, 3, 0); temp->Write(1, 4, 0);
		temp->Write(2, 0, 0); temp->Write(2, 1, 0); temp->Write(2, 2,-3); temp->Write(2, 3, 0); temp->Write(2, 4, 0);
		temp->Write(3, 0, 0); temp->Write(3, 1, 0); temp->Write(3, 2, 1); temp->Write(3, 3, 0); temp->Write(3, 4, 0);
		temp->Write(4, 0, 0); temp->Write(4, 1, 0); temp->Write(4, 2, 0); temp->Write(4, 3, 0); temp->Write(4, 4, 0);
		kerU = temp;

		// Down Kernel
		temp = new mat2D<int>(5, 5);
		temp->Write(0, 0, 0); temp->Write(0, 1, 0); temp->Write(0, 2, 0); temp->Write(0, 3, 0); temp->Write(0, 4, 0);
		temp->Write(1, 0, 0); temp->Write(1, 1, 0); temp->Write(1, 2, 1); temp->Write(1, 3, 0); temp->Write(1, 4, 0);
		temp->Write(2, 0, 0); temp->Write(2, 1, 0); temp->Write(2, 2,-3); temp->Write(2, 3, 0); temp->Write(2, 4, 0);
		temp->Write(3, 0, 0); temp->Write(3, 1, 0); temp->Write(3, 2, 3); temp->Write(3, 3, 0); temp->Write(3, 4, 0);
		temp->Write(4, 0, 0); temp->Write(4, 1, 0); temp->Write(4, 2,-1); temp->Write(4, 3, 0); temp->Write(4, 4, 0);
		kerD = temp;

		// Right Up Kernel
		temp = new mat2D<int>(5, 5);
		temp->Write(0, 0, 0); temp->Write(0, 1, 0); temp->Write(0, 2, 0); temp->Write(0, 3, 0); temp->Write(0, 4,-1);
		temp->Write(1, 0, 0); temp->Write(1, 1, 0); temp->Write(1, 2, 0); temp->Write(1, 3, 3); temp->Write(1, 4, 0);
		temp->Write(2, 0, 0); temp->Write(2, 1, 0); temp->Write(2, 2,-3); temp->Write(2, 3, 0); temp->Write(2, 4, 0);
		temp->Write(3, 0, 0); temp->Write(3, 1, 1); temp->Write(3, 2, 0); temp->Write(3, 3, 0); temp->Write(3, 4, 0);
		temp->Write(4, 0, 0); temp->Write(4, 1, 0); temp->Write(4, 2, 0); temp->Write(4, 3, 0); temp->Write(4, 4, 0);
		kerRU = temp;

		// Right Down Kernel
		temp = new mat2D<int>(5, 5);
		temp->Write(0, 0, 0); temp->Write(0, 1, 0); temp->Write(0, 2, 0); temp->Write(0, 3, 0); temp->Write(0, 4, 0);
		temp->Write(1, 0, 0); temp->Write(1, 1, 1); temp->Write(1, 2, 0); temp->Write(1, 3, 0); temp->Write(1, 4, 0);
		temp->Write(2, 0, 0); temp->Write(2, 1, 0); temp->Write(2, 2,-3); temp->Write(2, 3, 0); temp->Write(2, 4, 0);
		temp->Write(3, 0, 0); temp->Write(3, 1, 0); temp->Write(3, 2, 0); temp->Write(3, 3, 3); temp->Write(3, 4, 0);
		temp->Write(4, 0, 0); temp->Write(4, 1, 0); temp->Write(4, 2, 0); temp->Write(4, 3, 0); temp->Write(4, 4,-1);
		kerRD = temp;

		// Left Up Kernel
		temp = new mat2D<int>(5, 5);
		temp->Write(0, 0,-1); temp->Write(0, 1, 0); temp->Write(0, 2, 0); temp->Write(0, 3, 0); temp->Write(0, 4, 0);
		temp->Write(1, 0, 0); temp->Write(1, 1, 3); temp->Write(1, 2, 0); temp->Write(1, 3, 0); temp->Write(1, 4, 0);
		temp->Write(2, 0, 0); temp->Write(2, 1, 0); temp->Write(2, 2,-3); temp->Write(2, 3, 0); temp->Write(2, 4, 0);
		temp->Write(3, 0, 0); temp->Write(3, 1, 0); temp->Write(3, 2, 0); temp->Write(3, 3, 1); temp->Write(3, 4, 0);
		temp->Write(4, 0, 0); temp->Write(4, 1, 0); temp->Write(4, 2, 0); temp->Write(4, 3, 0); temp->Write(4, 4, 0);
		kerLU = temp;

		// Left Down Kernel
		temp = new mat2D<int>(5, 5);
		temp->Write(0, 0, 0); temp->Write(0, 1, 0); temp->Write(0, 2, 0); temp->Write(0, 3, 0); temp->Write(0, 4, 0);
		temp->Write(1, 0, 0); temp->Write(1, 1, 0); temp->Write(1, 2, 0); temp->Write(1, 3, 1); temp->Write(1, 4, 0);
		temp->Write(2, 0, 0); temp->Write(2, 1, 0); temp->Write(2, 2,-3); temp->Write(2, 3, 0); temp->Write(2, 4, 0);
		temp->Write(3, 0, 0); temp->Write(3, 1, 3); temp->Write(3, 2, 0); temp->Write(3, 3, 0); temp->Write(3, 4, 0);
		temp->Write(4, 0,-1); temp->Write(4, 1, 0); temp->Write(4, 2, 0); temp->Write(4, 3, 0); temp->Write(4, 4, 0);
		kerLD = temp;
	}

	s3(std::vector<float> qs, Config *config) : s(qs, config)
	{
		this->CreateKernels();
		quantMultiplier = 3;
		endOfSubmodel_s = new int[(int)qs.size()];
		endOfSubmodel_c = new int[(int)qs.size()];

		for (int qIndex = 0; qIndex < (int)qs.size(); qIndex++)
		{
			float q = qs[qIndex];
			std::vector<Submodel *> submodelsForQ;

			endOfSubmodel_c[qIndex] = 0;
			if (config->enable_c)
			{
				submodelsForQ.push_back(new s3_minmax22c(q, config));
				submodelsForQ.push_back(new s3_spam14c(q, config));
				submodelsForQ.push_back(new s3_minmax24c(q, config));
				submodelsForQ.push_back(new s3_minmax34hvc(q, config));
				submodelsForQ.push_back(new s3_minmax41c(q, config));
				submodelsForQ.push_back(new s3_minmax34c(q, config));
				submodelsForQ.push_back(new s3_minmax48c(q, config));
				submodelsForQ.push_back(new s3_minmax54c(q, config));
				endOfSubmodel_c[qIndex] = (int)submodelsForQ.size();
			}

			endOfSubmodel_s[qIndex] = endOfSubmodel_c[qIndex];
			if (config->enable_s)
			{
				submodelsForQ.push_back(new s3_minmax22h(q, config));
				submodelsForQ.push_back(new s3_minmax22v(q, config));
				submodelsForQ.push_back(new s3_minmax24(q, config));
				submodelsForQ.push_back(new s3_minmax34h(q, config));
				submodelsForQ.push_back(new s3_minmax34v(q, config));
				submodelsForQ.push_back(new s3_minmax41(q, config));
				submodelsForQ.push_back(new s3_minmax34(q, config));
				submodelsForQ.push_back(new s3_minmax48h(q, config));
				submodelsForQ.push_back(new s3_minmax48v(q, config));
				submodelsForQ.push_back(new s3_minmax54(q, config));
				submodelsForQ.push_back(new s3_spam14h(q, config));
				submodelsForQ.push_back(new s3_spam14v(q, config));
				endOfSubmodel_s[qIndex] = (int)submodelsForQ.size();
			}

			this->submodels.push_back(submodelsForQ);
		}
	}

	~s3()
	{
		delete kerR; delete kerL; delete kerU; delete kerD;
		delete kerRU; delete kerRD; delete kerLU; delete kerLD;
	}

	void ComputeImage(mat2D<int> **img, mat2D<double> * map, mat2D<int> **parity)
	{
		mat2D<int> *R[3] = { NULL, NULL, NULL }; GetResidual(img, kerR, R);
		mat2D<int> *L[3] = { NULL, NULL, NULL }; GetResidual(img, kerL, L);
		mat2D<int> *U[3] = { NULL, NULL, NULL }; GetResidual(img, kerU, U);
		mat2D<int> *D[3] = { NULL, NULL, NULL }; GetResidual(img, kerD, D);
		mat2D<int> *RU[3] = { NULL, NULL, NULL }; GetResidual(img, kerRU, RU);
		mat2D<int> *RD[3] = { NULL, NULL, NULL }; GetResidual(img, kerRD, RD);
		mat2D<int> *LU[3] = { NULL, NULL, NULL }; GetResidual(img, kerLU, LU);
		mat2D<int> *LD[3] = { NULL, NULL, NULL }; GetResidual(img, kerLD, LD);

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
			mat2D<int> *QRes_RU[3] = { NULL, NULL, NULL }; Quantize(RU, q, QRes_RU);
			mat2D<int> *QRes_RD[3] = { NULL, NULL, NULL }; Quantize(RD, q, QRes_RD);
			mat2D<int> *QRes_LU[3] = { NULL, NULL, NULL }; Quantize(LU, q, QRes_LU);
			mat2D<int> *QRes_LD[3] = { NULL, NULL, NULL }; Quantize(LD, q, QRes_LD);

			int submodelIndex = 0;
			std::vector<mat2D<int> **> QTResVect;
			mat2D<int> *QTRes_R[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_R);
			mat2D<int> *QTRes_L[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_L);
			mat2D<int> *QTRes_U[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_U);
			mat2D<int> *QTRes_D[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_D);
			mat2D<int> *QTRes_RU[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_RU);
			mat2D<int> *QTRes_RD[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_RD);
			mat2D<int> *QTRes_LU[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_LU);
			mat2D<int> *QTRes_LD[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_LD);

			if (config->enable_c)
			{
				Truncate(QRes_R, config->T_c, QTRes_R);
				Truncate(QRes_L, config->T_c, QTRes_L);
				Truncate(QRes_U, config->T_c, QTRes_U);
				Truncate(QRes_D, config->T_c, QTRes_D);
				Truncate(QRes_RU, config->T_c, QTRes_RU);
				Truncate(QRes_RD, config->T_c, QTRes_RD);
				Truncate(QRes_LU, config->T_c, QTRes_LU);
				Truncate(QRes_LD, config->T_c, QTRes_LD);
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
				Truncate(QRes_RU, config->T_s, QTRes_RU);
				Truncate(QRes_RD, config->T_s, QTRes_RD);
				Truncate(QRes_LU, config->T_s, QTRes_LU);
				Truncate(QRes_LD, config->T_s, QTRes_LD);
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
				delete QRes_RU[channel]; delete QRes_RD[channel]; delete QRes_LU[channel]; delete QRes_LD[channel];
			}
		}

		for (int channel = 0; channel < 3; channel++)
		{
			delete R[channel]; delete L[channel]; delete U[channel]; delete D[channel];
			delete RU[channel]; delete RD[channel]; delete LU[channel]; delete LD[channel];
		}
		delete pMap;
	}

private:
	mat2D<int> *kerR;
	mat2D<int> *kerL;
	mat2D<int> *kerU;
	mat2D<int> *kerD;
	mat2D<int> *kerRU;
	mat2D<int> *kerRD;
	mat2D<int> *kerLU;
	mat2D<int> *kerLD;
};