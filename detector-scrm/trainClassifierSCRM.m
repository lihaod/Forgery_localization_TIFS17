function trainClassifierSCRM(img_list,fea_path,trained_model_path)
% TRAINCLASSIFIERSCRM(trainImg_path,fea_path,trained_model_path) extracts
% SCRM features from training images and trains a detector.
%   fea_path - the directory for saving extracted features.
%   trained_model_path - the directory for saving trained model.

if ~exist(fea_path,'dir')
    mkdir(fea_path);
end

pool_size = 80; 
blk_sz = 32; 
stp_sz = 16; 
extFeaFromBlock(img_list,fea_path,blk_sz,stp_sz,pool_size,true);

minMdfRate = 0.1;
maxMdfRate = 0.9;
mkTrainingSet(fea_path,minMdfRate,maxMdfRate)

addpath(fea_path);

settings.mode = 1;
settings.trnRate = 1;
settings.isShuffle = false;
settings.verbose = 1;
settings.saveModel = true;
settings.saveModelPath = trained_model_path;
settings.saveModelName = 'scrm';
settings.saveResult = false;

ensembleTrnTst({'au.mat'},{'sp.mat'},settings);

rmpath(fea_path);

