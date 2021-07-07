function Fea = scrmq18157(IMAGE)

if ischar(IMAGE)
    I = imread(IMAGE);
else
    I = IMAGE;
end

config.enable_s = true;
config.T_s = int32(2);
config.order_s = int32(4);
config.enable_c = true;
config.T_c = int32(3);
config.order_c = int32(3);
config.q1 = 1;
config.qothers = 1;
config.roundup5 = true;

fea_struct = SCRM(I,config);
fea_cell = struct2cell(fea_struct);
fea_cell = fea_cell([1:8 19:23 29:36 47:50 55:58 66 9:18 63 24:28 64 37:46 65 51:54 67 59:62 68 69]);

Fea = cell2mat(transpose(fea_cell));
