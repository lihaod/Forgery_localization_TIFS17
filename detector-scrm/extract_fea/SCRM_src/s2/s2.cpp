/* 
 This class computes and contains co-occurrences of ALL 2nd-order residuals
 listed in Figure 1 in our journal HUGO paper (version from June 14), 
 including the naming convention.

 List of outputted features:

 1a) spam12h
 1b) spam12v (orthogonal-spam)
 1c) minmax21
 1d) minmax41
 1e) minmax24h (24v is also outputted but not listed in Figure 1)
 1f) minmax32

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

#include "s2_spam12h.cpp"
#include "s2_spam12v.cpp"
#include "s2_minmax21.cpp"
#include "s2_minmax24h.cpp"
#include "s2_minmax24v.cpp"
#include "s2_minmax32.cpp"
#include "s2_minmax41.cpp"
#include "s2_spam12c.cpp"
#include "s2_minmax21c.cpp"
#include "s2_minmax24c.cpp"
#include "s2_minmax32c.cpp"
#include "s2_minmax41c.cpp"

class s2 : s
{
public:
	void CreateKernels()
	{
		mat2D<int> *temp;
		cutEdgesForParityBy = 1;
		
		// Horizontal Kernel
		temp = new mat2D<int>(3, 3);
		temp->Write(0, 0, 0); temp->Write(0, 1, 0); temp->Write(0, 2, 0);
		temp->Write(1, 0, 1); temp->Write(1, 1,-2); temp->Write(1, 2, 1);
		temp->Write(2, 0, 0); temp->Write(2, 1, 0); temp->Write(2, 2, 0);
		kerH = temp;

		// Vertical Kernel
		temp = new mat2D<int>(3, 3);
		temp->Write(0, 0, 0); temp->Write(0, 1, 1); temp->Write(0, 2, 0);
		temp->Write(1, 0, 0); temp->Write(1, 1,-2); temp->Write(1, 2, 0);
		temp->Write(2, 0, 0); temp->Write(2, 1, 1); temp->Write(2, 2, 0);
		kerV = temp;

		// Diagonal Kernel
		temp = new mat2D<int>(3, 3);
		temp->Write(0, 0, 1); temp->Write(0, 1, 0); temp->Write(0, 2, 0);
		temp->Write(1, 0, 0); temp->Write(1, 1,-2); temp->Write(1, 2, 0);
		temp->Write(2, 0, 0); temp->Write(2, 1, 0); temp->Write(2, 2, 1);
		kerD = temp;

		// Minor diagonal Kernel
		temp = new mat2D<int>(3, 3);
		temp->Write(0, 0, 0); temp->Write(0, 1, 0); temp->Write(0, 2, 1);
		temp->Write(1, 0, 0); temp->Write(1, 1,-2); temp->Write(1, 2, 0);
		temp->Write(2, 0, 1); temp->Write(2, 1, 0); temp->Write(2, 2, 0);
		kerM = temp;
	}

	s2(std::vector<float> qs, Config *config) : s(qs, config)
	{
		this->CreateKernels();
		quantMultiplier = 2;
		endOfSubmodel_s = new int[(int)qs.size()];
		endOfSubmodel_c = new int[(int)qs.size()];

		for (int qIndex=0; qIndex < (int)qs.size(); qIndex++)
		{
			float q = qs[qIndex];
			std::vector<Submodel *> submodelsForQ;

			endOfSubmodel_c[qIndex] = 0;
			if (config->enable_c)
			{
				submodelsForQ.push_back(new s2_spam12c(q, config));
				submodelsForQ.push_back(new s2_minmax21c(q, config));
				submodelsForQ.push_back(new s2_minmax41c(q, config));
				submodelsForQ.push_back(new s2_minmax32c(q, config));
				submodelsForQ.push_back(new s2_minmax24c(q, config));
				endOfSubmodel_c[qIndex] = (int)submodelsForQ.size();
			}

			endOfSubmodel_s[qIndex] = endOfSubmodel_c[qIndex];
			if (config->enable_s)
			{
				
				submodelsForQ.push_back(new s2_minmax21(q, config));
				submodelsForQ.push_back(new s2_minmax41(q, config));
				submodelsForQ.push_back(new s2_minmax32(q, config));
				submodelsForQ.push_back(new s2_minmax24h(q, config));
				submodelsForQ.push_back(new s2_minmax24v(q, config));
				submodelsForQ.push_back(new s2_spam12h(q, config));
				submodelsForQ.push_back(new s2_spam12v(q, config));
				endOfSubmodel_s[qIndex] = (int)submodelsForQ.size();
			}

			this->submodels.push_back(submodelsForQ);
		}
	}

	~s2()
	{
		delete kerH; delete kerV; delete kerD; delete kerM;
	}

	void ComputeImage(mat2D<int> **img, mat2D<double> * map, mat2D<int> **parity)
	{
		mat2D<int> *H[3] = { NULL, NULL, NULL }; GetResidual(img, kerH, H);
		mat2D<int> *V[3] = { NULL, NULL, NULL }; GetResidual(img, kerV, V);
		mat2D<int> *D[3] = { NULL, NULL, NULL }; GetResidual(img, kerD, D);
		mat2D<int> *M[3] = { NULL, NULL, NULL }; GetResidual(img, kerM, M);

		mat2D<double> * pMap = new mat2D<double>(img[0]->rows - 2, img[0]->cols - 2);
		for (int i = 0; i<img[0]->rows - 2; i++)
			for (int j = 0; j<img[0]->cols - 2; j++)
				pMap->Write(i, j, map->Read(i + 1, j + 1));
		for (int i = 0; i<this->submodels.size(); i++)
			for (int j = 0; j<this->submodels[i].size(); j++)
				this->submodels[i][j]->map = pMap;

		for (int qIndex=0; qIndex < (int)submodels.size(); qIndex++)
		{
			float q = qs[qIndex] * quantMultiplier;

			mat2D<int> *QRes_H[3] = { NULL, NULL, NULL }; Quantize(H, q, QRes_H);
			mat2D<int> *QRes_V[3] = { NULL, NULL, NULL }; Quantize(V, q, QRes_V);
			mat2D<int> *QRes_D[3] = { NULL, NULL, NULL }; Quantize(D, q, QRes_D);
			mat2D<int> *QRes_M[3] = { NULL, NULL, NULL }; Quantize(M, q, QRes_M);


			int submodelIndex = 0;
			std::vector<mat2D<int> **> QTResVect;
			mat2D<int> *QTRes_H[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_H);
			mat2D<int> *QTRes_V[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_V);
			mat2D<int> *QTRes_D[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_D);
			mat2D<int> *QTRes_M[3] = { NULL, NULL, NULL }; QTResVect.push_back(QTRes_M);

			if (config->enable_c)
			{
				Truncate(QRes_H, config->T_c, QTRes_H);
				Truncate(QRes_V, config->T_c, QTRes_V);
				Truncate(QRes_D, config->T_c, QTRes_D);
				Truncate(QRes_M, config->T_c, QTRes_M);
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
				Truncate(QRes_H, config->T_s, QTRes_H);
				Truncate(QRes_V, config->T_s, QTRes_V);
				Truncate(QRes_D, config->T_s, QTRes_D);
				Truncate(QRes_M, config->T_s, QTRes_M);
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
				delete QRes_H[channel]; delete QRes_V[channel]; delete QRes_D[channel]; delete QRes_M[channel];
			}
		}

		for (int channel = 0; channel < 3; channel++)
		{
			delete H[channel]; delete V[channel]; delete D[channel]; delete M[channel];
		}
		delete pMap;
	}

private:
	mat2D<int> *kerH;
	mat2D<int> *kerV;
	mat2D<int> *kerD;
	mat2D<int> *kerM;
};