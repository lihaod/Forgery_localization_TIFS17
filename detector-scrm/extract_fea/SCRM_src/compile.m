% compile mex
clear all
outputName = 'SCRM';
debugFlag = 0;

try
    delete([outputName '.' mexext]);
    delete([outputName '.' mexext '.pdb']);
catch exception
    rethrow(exception)
end

if debugFlag
    mex('-O', '-g', '-largeArrayDims', '-output', outputName, ...
        'SCRM_matlab.cpp', 'SCRMclass.cpp', 'submodel.cpp', 's.cpp', 'config.cpp', ...
        '-I../include', ['-D' computer]);
else
    mex('-O', '-largeArrayDims', '-output', outputName, ...
        'SCRM_matlab.cpp', 'SCRMclass.cpp', 'submodel.cpp', 's.cpp', 'config.cpp', ...
        '-I../include', ['-D' computer]);
end

