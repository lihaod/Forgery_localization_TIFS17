function response_map = fusionMap(map_fea,map_pm)

lambda1 = 0.39;
lambda2 = 4.26;
thr = 0.48;
fusionFun = @(a,b)(imfilter(mat2gray(a).^lambda1,ones(64)/64^2,'symmetric').*imfilter((b).^lambda2,ones(64)/64^2,'symmetric'));

response_map = struct();
response_map.map = fusionFun(map_fea,map_pm)>=thr;