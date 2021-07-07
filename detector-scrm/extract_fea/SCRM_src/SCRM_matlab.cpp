//#define CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>

#ifdef GLNXA64
#include <stdio.h>
#include <cstring>
#endif

#include <vector>
#include "submodel.h"
#include "SCRMclass.h"

#include <mex.h>


/*
prhs[0] - cell array of image paths
prhs[1] - struct config
config.T			- int32		- default 2				- residual threshold
config.order		- int32		- default 4				- co-occurrence order
config.q1			- vector	- default [1, 2]		- quantization steps for 1st order
config.qothers		- vector	- default [1, 1.5, 2]	- quantization steps for other orders
config.merge_spams	- logical	- default true			- if true then spam features are merged
config.symm_sign	- logical	- default true			- if true then spam symmetry is used
config.symm_reverse	- logical	- default true			- if true then reverse symmetry is used
config.symm_minmax	- logical	- default true			- if true then minmax symmetry is used
config.eraseLSB		- logical	- default false			- if true then all LSB are erased from the image
config.parity		- logical	- default false			- if true then parity residual is applied
config.roundup5		- logical	- default false			- if true then 0.5 will be rounded to 1.0
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

	if ((nrhs != 1) && (nrhs != 2) && (nrhs != 3))
		mexErrMsgTxt("One or two inputs are required.\ninput #1 - [uint8_3-D_array Image]\ninput #2 - [struct config]");

	const mxArray *imageSet = prhs[0];
	const char *classofimageSet = mxGetClassName(imageSet);
	const mwSize imgDimNum = mxGetNumberOfDimensions(imageSet);
	const mwSize *imgDim = mxGetDimensions(imageSet);

	if (strcmp(classofimageSet, "uint8") != 0 || imgDimNum != 3 || imgDim[2] != 3)
		mexErrMsgTxt("The first input must be a uint8 COLOR image matrix.");

	// Default config
	bool enable_c = true;
	bool enable_s = true;
	int T_s = 2;
	int order_s = 4;
	int T_c = 3;
	int order_c = 3;
	bool mergeSpams = true;
	bool ss = true, sr = true, sm = true;
	bool eraseLSB = false, parity = false;
	bool roundup5 = true;
	bool mapMode = true;
	int numInMap = 1;
	mat2D<double> *c_map = NULL;
	std::vector<float> q_1 = std::vector<float>();
	std::vector<float> q_others = std::vector<float>();

	int nfields = mxGetNumberOfFields(prhs[1]);
	if (nfields == 0) mexErrMsgTxt("The config structure is empty.");
	for (int fieldIndex = 0; fieldIndex<nfields; fieldIndex++)
	{
		const char *fieldName = mxGetFieldNameByNumber(prhs[1], fieldIndex);
		const mxArray *fieldContent = mxGetFieldByNumber(prhs[1], 0, fieldIndex);

		if (strcmp(fieldName, "q1") == 0)
			if (mxIsDouble(fieldContent))
				for (int i = 0; i < mxGetNumberOfElements(fieldContent); i++)
					q_1.push_back((float)mxGetPr(fieldContent)[i]);
			else mexErrMsgTxt("'config.q1' must be a double scalar or vector");
		else if (strcmp(fieldName, "qothers") == 0)
			if (mxIsDouble(fieldContent))
				for (int i = 0; i < mxGetNumberOfElements(fieldContent); i++)
					q_others.push_back((float)mxGetPr(fieldContent)[i]);
			else mexErrMsgTxt("'config.qothers' must be a double scalar or vector");
		else if (strcmp(fieldName, "map") == 0)
			if (mxIsDouble(fieldContent))
			{
				int rows = (int)mxGetM(fieldContent);
				int cols = (int)mxGetN(fieldContent);
				c_map = new mat2D<double>(rows, cols);
				double *map_array = (double *)mxGetData(fieldContent);
				for (int c = 0; c<cols; c++)
					for (int r = 0; r<rows; r++)
						c_map->Write(r, c, map_array[r + c*rows]);
			}
			else mexErrMsgTxt("'config.map' must be a double matrix");
		// if a field is not scalar
		else if ((mxGetM(fieldContent) != 1) || (mxGetN(fieldContent) != 1))
			mexErrMsgTxt("All config fields must be scalars, except for 'q1', 'qothers' and 'map'.");
		// if every field is scalar
		else if (strcmp(fieldName, "enable_s") == 0)
			if (mxIsLogical(fieldContent)) enable_s = mxIsLogicalScalarTrue(fieldContent);
			else mexErrMsgTxt("'config.enable_s' must be of type 'logical'");
		else if (strcmp(fieldName, "enable_c") == 0)
			if (mxIsLogical(fieldContent)) enable_c = mxIsLogicalScalarTrue(fieldContent);
			else mexErrMsgTxt("'config.enable_c' must be of type 'logical'");
		else if (strcmp(fieldName, "T_s") == 0)
			if (mxIsClass(fieldContent, "int32")) T_s = (int)mxGetScalar(fieldContent);
			else mexErrMsgTxt("'config.T_s' must be of type 'int32'");
		else if (strcmp(fieldName, "order_s") == 0)
			if (mxIsClass(fieldContent, "int32")) order_s = (int)mxGetScalar(fieldContent);
			else mexErrMsgTxt("'config.order_s' must be of type 'int32'");
		else if (strcmp(fieldName, "T_c") == 0)
			if (mxIsClass(fieldContent, "int32")) T_c = (int)mxGetScalar(fieldContent);
			else mexErrMsgTxt("'config.T_c' must be of type 'int32'");
		else if (strcmp(fieldName, "order_c") == 0)
			if (mxIsClass(fieldContent, "int32")) order_c = (int)mxGetScalar(fieldContent);
			else mexErrMsgTxt("'config.order_c' must be of type 'int32'");
		else if (strcmp(fieldName, "merge_spams") == 0)
			if (mxIsLogical(fieldContent)) mergeSpams = mxIsLogicalScalarTrue(fieldContent);
			else mexErrMsgTxt("'config.mergeSpams' must be of type 'logical'");
		else if (strcmp(fieldName, "symm_sign") == 0)
			if (mxIsLogical(fieldContent)) ss = mxIsLogicalScalarTrue(fieldContent);
			else mexErrMsgTxt("'config.symm_sign' must be of type 'logical'");
		else if (strcmp(fieldName, "symm_reverse") == 0)
			if (mxIsLogical(fieldContent)) sr = mxIsLogicalScalarTrue(fieldContent);
			else mexErrMsgTxt("'config.symm_reverse' must be of type 'logical'");
		else if (strcmp(fieldName, "symm_minmax") == 0)
			if (mxIsLogical(fieldContent)) sm = mxIsLogicalScalarTrue(fieldContent);
			else mexErrMsgTxt("'config.symm_minmax' must be of type 'logical'");
		else if (strcmp(fieldName, "eraseLSB") == 0)
			if (mxIsLogical(fieldContent)) eraseLSB = mxIsLogicalScalarTrue(fieldContent);
			else mexErrMsgTxt("'config.eraseLSB' must be of type 'logical'");
		else if (strcmp(fieldName, "parity") == 0)
			if (mxIsLogical(fieldContent)) parity = mxIsLogicalScalarTrue(fieldContent);
			else mexErrMsgTxt("'config.parity' must be of type 'logical'");
		else if (strcmp(fieldName, "roundup5") == 0)
			if (mxIsLogical(fieldContent)) roundup5 = mxIsLogicalScalarTrue(fieldContent);
			else mexErrMsgTxt("'config.roundup5' must be of type 'logical'");
		else if (strcmp(fieldName, "mapIsWeights") == 0)
			if (mxIsLogical(fieldContent)) mapMode = mxIsLogicalScalarTrue(fieldContent);
			else mexErrMsgTxt("'config.mapModeIsWeights' must be of type 'logical', 1 for 'weights', 0 for 'segments'");
		else if (strcmp(fieldName, "numInMap") == 0)
			if (mxIsClass(fieldContent, "int32")) numInMap = (int)mxGetScalar(fieldContent);
			else mexErrMsgTxt("'config.numInMap' must be of type 'int32'");
	}

	// if no map is given, set mapMode = true
	if (NULL == c_map) mapMode = true;
	// if map is weights, set numInMap as 1
	if (mapMode) numInMap = 1;

	// set default quantization steps if they are not specified
	if (q_1.size() == 0){ q_1.push_back(1); q_1.push_back(2); }
	if (q_others.size() == 0){ q_others.push_back(1); q_others.push_back(1.5); q_others.push_back(2); }

	// create config object
	Config *config = new Config(false, enable_s, enable_c, T_s, T_c, order_s, order_c, q_1, q_others, ss, sr, sm, mergeSpams, eraseLSB, parity, roundup5, mapMode, numInMap);

	// create object with all the submodels and compute the features
	SCRMclass *SCRMobj = new SCRMclass(config);

	
	// convert mxArray *imageSet to mat2D
	mat2D<int> *c_image[3];
	const unsigned char *image_array = (unsigned char *)mxGetData(imageSet);
	for (int channel = 0; channel < 3; channel++)
	{
		c_image[channel] = new mat2D<int>((int)imgDim[0], (int)imgDim[1]);
		for (int c = 0; c < imgDim[1]; c++)
			for (int r = 0; r < imgDim[0]; r++)
				c_image[channel]->Write(r, c, image_array[channel*imgDim[0] * imgDim[1] + c*imgDim[0] + r]);
	}

	// remove LSBs from the image if in config
	if (eraseLSB)
	for (int channel = 0; channel < 3; channel++)
		for (int r = 0; r<imgDim[0]; r++)
			for (int c = 0; c<imgDim[1]; c++)
			{
				int iValue = c_image[channel]->Read(r, c);
				iValue = iValue - (iValue % 2);
				c_image[channel]->Write(r, c, iValue);
			}

	if (NULL == c_map)
	{
		c_map = new mat2D<double>((int)imgDim[0], (int)imgDim[1]);
		for (int c = 0; c<imgDim[1]; c++)
			for (int r = 0; r<imgDim[0]; r++)
				c_map->Write(r, c, 1);
	}

	// Run the feature computation
	SCRMobj->ComputeFeatures(c_image, c_map);

	for (int channel = 0; channel < 3; channel++) delete c_image[channel];
	delete c_map;

	std::vector<Submodel *> submodels = SCRMobj->GetSubmodels();
	const char **submodelNames = new const char*[submodels.size()];
	std::string *submodelNamesAsString = new std::string[submodels.size()];
	for (int submodelIndex = 0; submodelIndex < (int)submodels.size(); submodelIndex++)
	{
		submodelNamesAsString[submodelIndex] = submodels[submodelIndex]->GetName();
		submodelNames[submodelIndex] = submodelNamesAsString[submodelIndex].c_str();
	}
	mwSize structSize[2];
	structSize[0] = 1;
	structSize[1] = 1;
	plhs[0] = mxCreateStructArray(1, structSize, (int)submodels.size(), submodelNames);
	for (int submodelIndex = 0; submodelIndex < submodels.size(); submodelIndex++)
	{
		Submodel *currentSubmodel = submodels[submodelIndex];
		mwSize feaSize[3];
		feaSize[0] = (int)currentSubmodel->ReturnFea().size();
		feaSize[1] = currentSubmodel->symmDim;
		feaSize[2] = currentSubmodel->numInMap;
		mxArray *fea = mxCreateNumericArray(3, feaSize, mxSINGLE_CLASS, mxREAL);
		for (int r = 0; r < (int)feaSize[0]; r++)
			for (int c = 0; c < currentSubmodel->symmDim; c++)
				for (int h = 0; h < currentSubmodel->numInMap; h++)
					((float*)mxGetPr(fea))[(h*(int)feaSize[0] * (int)feaSize[1]) + (c*(int)feaSize[0]) + r] = (currentSubmodel->ReturnFea())[r][c + h*(int)feaSize[1]];
		mxSetFieldByNumber(plhs[0], 0, submodelIndex, fea);
	}

	delete[] submodelNamesAsString;
	delete[] submodelNames;

	delete config;
	delete SCRMobj;
}
