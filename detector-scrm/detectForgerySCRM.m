function response_map = detectForgerySCRM(image, trained_model, blk_sz, stp_sz)
% Output the response map using sliding window with blk_sz*blk_sz and stp_sz,
% trained_model  - the path of trained classifier.
% blk_sz[=64] - the size of sliding window.
% stp_sz[=16] - the size of step.

narginchk(2,4);

if ischar(image)
    image = imread(image);
end
if nargin < 4; stp_sz = 16; end
if nargin < 3; blk_sz = 64; end

trained_ensemble = load(trained_model,'trained_ensemble');
trained_ensemble = trained_ensemble.trained_ensemble;
maxVote = length(trained_ensemble);

a = size(image,1);
b = size(image,2);
vote_map = zeros(a, b);
correct_array = zeros(a,b);

h_num = floor((a-blk_sz)/stp_sz)+1;
w_num = floor((b-blk_sz)/stp_sz)+1;
    
block_fea = zeros(h_num*w_num, 18157);
    
for index_fea = 1:h_num*w_num
    i = floor((index_fea-1)/w_num) + 1;
    j = mod((index_fea-1),w_num) + 1;
    img_block = image((i-1)*stp_sz+(1:blk_sz),(j-1)*stp_sz+(1:blk_sz), :);
    block_fea(index_fea,:) = scrmq18157(img_block);
end

test_result = ensemble_testing(block_fea,trained_ensemble);
votes = test_result.votes;

blk_index = 0;
for i = 1:stp_sz:a-blk_sz+1
    for j = 1:stp_sz:b-blk_sz+1
        blk_index = blk_index+1;
        vote_map(i:i+blk_sz-1,j:j+blk_sz-1) = vote_map(i:i+blk_sz-1,j:j+blk_sz-1) - votes(blk_index);
        correct_array(i:i+blk_sz-1,j:j+blk_sz-1) = correct_array(i:i+blk_sz-1,j:j+blk_sz-1) +1;
    end
end

vote_map = vote_map./correct_array;
vote_map = (vote_map+maxVote)/maxVote/2;
vote_map(isnan(vote_map))=1;

response_map = struct();
response_map.map = vote_map;
end