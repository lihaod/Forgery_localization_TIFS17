function extFeaFromBlock(img_list,fea_path,blk_sz,stp_sz,pool_size,require_gt)

global nimble_17_data

narginchk(2,6);

if nargin < 6; require_gt = false; end
if nargin < 5; pool_size = 10; end
if nargin < 4; stp_sz = 16; end
if nargin < 3; blk_sz = 64; end

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj) || poolobj.NumWorkers~=pool_size
    delete(gcp('nocreate'));
    parpool(pool_size);
end

if ~exist(fea_path,'dir'); mkdir(fea_path); end

if ~exist('img_list', 'var') || isempty(img_list)
    img_list = listTestImages();
end
for index = 1:length(img_list)
    if exist(fullfile(fea_path, [img_list{index} '.mat']),'file')
        continue;
    end
    m = metaImage(img_list{index});
    m = m.loadSideImages();
    if require_gt && isempty(m.groundTruth)
        continue;
    end
    img = m.pixmap;
    [a,b,~] = size(img);
    
    h_num = floor((a-blk_sz)/stp_sz)+1;
    w_num = floor((b-blk_sz)/stp_sz)+1;
    
    block_fea = zeros(h_num*w_num, 18157);
    img_block = zeros(blk_sz, blk_sz, 3);
    if require_gt
        img_mask = m.groundTruthBin;
        block_modify_rate = zeros(h_num*w_num, 1);
    end
    
    parfor index_fea = 1:h_num*w_num
        % for index_fea = 1:h_num*w_num
        i = floor((index_fea-1)/w_num) + 1;
        j = mod((index_fea-1),w_num) + 1;
        img_block = img((i-1)*stp_sz+(1:blk_sz),(j-1)*stp_sz+(1:blk_sz), :);
        block_fea(index_fea,:) = scrmq18157(img_block);
        if require_gt
            img_mask_block = img_mask((i-1)*stp_sz+(1:blk_sz),(j-1)*stp_sz+(1:blk_sz));
            block_modify_rate(index_fea, 1) = sum(sum(img_mask_block==1))/numel(img_mask_block);
        end
    end
    
    if require_gt
        save(fullfile(fea_path, [img_list{index} '.mat']),'block_modify_rate', 'block_fea','-v7.3'); 
    else
        save(fullfile(fea_path, [img_list{index} '.mat']), 'block_fea','-v7.3');
    end
    
end

delete(gcp('nocreate'))