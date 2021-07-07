%%
% compile scripts will run 'clear all'!!!

cd SCRM_src
compile
copyfile(['SCRM.' mexext],'..');
cd('..')

% cd SRM_src
% compile
% copyfile(['SRM.' mexext],'..');
% cd('..')
% 
% cd CRM_src
% compile
% copyfile(['CRM.' mexext],'..');
% cd('..')
