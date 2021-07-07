%%
clear
addpath('detector-scrm\','detector-scrm\ensemble\','detector-scrm\extract_fea\',...
    'detector-pm\','detector-pm\patchmatch-2.1\');
imgName = 'aa61a96b0a18b8dbc65fd20af3644958.png';
I = imread(['img\' imgName]);
GT = imread(['img\mask\' imgName]);
I_GT = uint8(cat(3,double(I(:,:,1))+128*(GT==0),double(I(:,:,2))-128*(GT==0),double(I(:,:,3))-128*(GT==0)));

M_fea = detectForgerySCRM(I, 'trained_ensemble_scrm.mat');
M_pm  = detectForgeryPM(I);
M_fus = fusionMap(M_fea.map,M_pm.map);
if ~isempty(GT)
    F1 = getMeasure(M_fus.map,GT);
end

figure
subplot(2,2,1),imshow(I); title('Input image');
subplot(2,2,2),imshow(I_GT); title('Ground truth');
subplot(2,3,4),imshow(M_fea.map); title('M^{Fea}');
subplot(2,3,5),imshow(M_pm.map); title('M^{PM}');
subplot(2,3,6),imshow(M_fus.map); title('Localization map');