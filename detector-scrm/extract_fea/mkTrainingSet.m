function mkTrainingSet(feaPath,minMdfRate,maxMdfRate)
% find out the patches with minMdfRate<block_modify_rate<maxMdfRate in an image, and ramdomly slect
% the same number of non-modified patches.

narginchk(1,3);
if nargin < 2; minMdfRate = 0.1; end
if nargin < 3; maxMdfRate = 0.9; end

fileList = dir(fullfile(feaPath,'*.mat'));


F_sp = [];
F_au = [];
shortfallNum = zeros(length(fileList),1);
surplusNum = zeros(length(fileList),1);
surplus = cell(length(fileList),1);
for k = 1:length(fileList)
    load(fullfile(feaPath,fileList(k).name),'block_fea','block_modify_rate');
    [instNum,feaDim] = size(block_fea);
    zeroMdfRate = find(block_modify_rate == 0);
    fea_zero = block_fea(zeroMdfRate,:);
    properInx_bin = block_modify_rate>=minMdfRate & block_modify_rate<=maxMdfRate;
    properInx = find(properInx_bin);
    F_sp = [F_sp;block_fea(properInx,:)];
    randInx = randperm(length(zeroMdfRate),min(length(zeroMdfRate),length(properInx)));
    shortfallNum(k) = length(properInx)-length(randInx);
    F_au = [F_au; fea_zero(randInx,:); zeros(shortfallNum(k),feaDim)];
    surplusNum(k) = length(zeroMdfRate)-length(randInx);
    surplus{k} = setdiff(zeroMdfRate,randInx);
    fprintf('#%4d: %s %6d %6d %6d %6d \n', k, fileList(k).name, instNum, length(properInx),shortfallNum(k), surplusNum(k));
end
clear zeroMdfRate fea_zero properInx randInx

fprintf('Fill the shortfall features...\n');
shortfallFeaInx = find(sum(F_au,2)==0);
if length(shortfallFeaInx)~= sum(shortfallNum)
    error('shortfall number inconsistent.');
end
selectInx = randperm(sum(surplusNum),sum(shortfallNum));
F_sf = zeros(sum(shortfallNum),feaDim,'single');
cumSurplusNum = [0;cumsum(surplusNum)];
for k = 1:length(fileList)
    if surplusNum(k)==0
        continue;
    end
    thisSelectInx = find(selectInx>cumSurplusNum(k)&selectInx<=cumSurplusNum(k+1));
    if ~isempty(thisSelectInx)
        load(fullfile(feaPath,fileList(k).name),'block_fea','block_modify_rate');
        for j = 1:length(thisSelectInx)
            F_sf(thisSelectInx(j),:) = block_fea(surplus{k}(selectInx(thisSelectInx(j))-cumSurplusNum(k)),:);
        end
        fprintf('#%4d: %s %6d\n', k, fileList(k).name, length(thisSelectInx));
    end
end
F_au(shortfallFeaInx,:) = F_sf;

F = F_au;
save(fullfile(feaPath,'au.mat'),'F','-v7.3');
F = F_sp;
save(fullfile(feaPath,'sp.mat'),'F','-v7.3');
fprintf('Finished\n');