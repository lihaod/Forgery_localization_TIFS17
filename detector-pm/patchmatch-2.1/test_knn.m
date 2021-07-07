clear
close all
warning('off','images:initSize:adjustingMag');

cores = 4;

b=imread('test.png');
% c=imread('c.png');

% % kNN without enrichment
% tic;
% cnn=nnmex(b, b, 'cputiled', 7, 16, [], [], [], [], cores, [], [], [], [], [], 4);
% toc
% imshow(cnn(:,:,1,1), []); figure
% imshow(cnn(:,:,1,2), []); figure
% imshow(cnn(:,:,1,3), []); figure
% imshow(cnn(:,:,1,4), []); figure
% D = sqrt(double(cnn(:,:,3,:)));
% format long;
% disp(['Average dist (no enrichment):', num2str(mean(D(:)))]);

% kNN with enrichment -- both images must be the same. Enrichment requires the number of NN iterations to be even -- if not it will round down to the next even number.
patch_w = 7;
tic;
cnn=nnmex(b, b, 'enrich', patch_w, 6, [], [], [], [], cores, [], [], [], [], [], 4);
toc
%%
coord_x = repmat(1:size(cnn,2)-patch_w+1,size(cnn,1)-patch_w+1,1);
coord_y = repmat((1:size(cnn,1)-patch_w+1).',1,size(cnn,2)-patch_w+1);

cnn_map = zeros([size(coord_x) size(cnn,4)]);
for k = 1:size(cnn,4)
    cnn_map(:,:,k) = sqrt((double(cnn(1:end-patch_w+1,1:end-patch_w+1,1,k))-coord_x).^2+(double(cnn(1:end-patch_w+1,1:end-patch_w+1,2,k))-coord_y).^2);
end
cnn_map = sum(cnn_map,3);
imshow(cnn_map,[]);
%%
D = sqrt(double(cnn(:,:,3,:)));
disp(['Average dist (enrichment):', num2str(mean(D(:)))]);
warning('on','images:initSize:adjustingMag');