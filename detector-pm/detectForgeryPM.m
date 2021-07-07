function response_map = detectForgeryPM(image,iter,cores,patch_w,knn,shrink_sz,wsize,thr_relOffset,thr_absOffset)

narginchk(1,9);

warning('off','images:initSize:adjustingMag');

if ischar(image)
    image = imread(image);
end

if nargin < 9; thr_absOffset = 8; end
if nargin < 8; thr_relOffset = 2; end
if nargin < 7; wsize = 25; end
if nargin < 6; shrink_sz = 1024; end
if nargin < 5; knn = 4; end
if nargin < 4; patch_w = 7; end
if nargin < 3; cores = 10; end
if nargin < 2; iter = 6; end
    
if size(image,3)==1; image = repmat(image,1,1,3); end
I_h = size(image,1); I_w = size(image,2);
I_maxsz = max(I_h,I_w);
if I_maxsz > shrink_sz
    image = imresize(image,shrink_sz/I_maxsz,'bicubic'); % shrink the image to prevent the PatchMatch algorithm crashing
end
height = size(image,1)-patch_w+1;
width = size(image,2)-patch_w+1;
cnn=nnmex(image, image, 'enrich', patch_w, iter, [], [], [], [], cores, [], [], [], [], [], knn);
[coord_x,coord_y] = meshgrid(0:width-1,0:height-1);
cnn_map = zeros([size(coord_x),2,size(cnn,4)]);
r = floor(wsize/2);

D = zeros(size(cnn_map,1),size(cnn_map,2));
for nn = 1:knn
    cnn_map(:,:,1,nn) = double(cnn(1:height,1:width,1,nn))-coord_x;
    cnn_map(:,:,2,nn) = double(cnn(1:height,1:width,2,nn))-coord_y;
    dx = cnn_map(:,:,1,nn);
    dy = cnn_map(:,:,2,nn);
    dx = padarray(dx,[r r]);
    dy = padarray(dy,[r r]);
    shield = (abs(dx)+abs(dy))<thr_absOffset;
    
    D = D + procCNNMap(dx,dy,shield,r,thr_relOffset,thr_absOffset);
end

D = D/nn;
map = padarray(D,[patch_w-1 patch_w-1],0,'post');
if I_maxsz > shrink_sz
    map = imresize(map,[I_h I_w],'bilinear');
end

response_map = struct();
response_map.map = 1-map;

