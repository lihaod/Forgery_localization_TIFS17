%%
mex procCNNMap.cpp
cd ./patchmatch-2.1
if strcmp('PCWIN',computer)||strcmp('PCWIN64',computer)
    mex OPTIMFLAGS="/Ox /Oi /Oy /DNDEBUG /fp:fast /arch:SSE2 /DMEX_MODE /openmp" mexutil.cpp nn.cpp nnmex.cpp patch.cpp vecnn.cpp simnn.cpp allegro_emu.cpp knn.cpp -output nnmex_mex
    mex OPTIMFLAGS="/Ox /Oi /Oy /DNDEBUG /fp:fast /arch:SSE2 /DMEX_MODE /openmp" mexutil.cpp nn.cpp votemex.cpp patch.cpp vecnn.cpp simnn.cpp allegro_emu.cpp knn.cpp -output votemex
elseif strcmp('GLNXA64',computer)
    mex CXXOPTIMFLAGS='-O6 -w -s -ffast-math -fomit-frame-pointer -fstrength-reduce -fopenmp -msse2 -funroll-loops -fPIC' CXXFLAGS='-DNDEBUG -DUNIX_MODE -DMEXMODE -fopenmp' CXXLIBS='${CXXLIBS} -Wl,--export-dynamic -Wl,-e,mexFunction -shared -lgomp' knn.cpp mexutil.cpp nn.cpp nnmex.cpp patch.cpp vecnn.cpp simnn.cpp allegro_emu.cpp -output nnmex -output nnmex
    mex CXXOPTIMFLAGS='-O6 -w -s -ffast-math -fomit-frame-pointer -fstrength-reduce -fopenmp -msse2 -funroll-loops -fPIC' CXXFLAGS='-DNDEBUG -DUNIX_MODE -DMEXMODE -fopenmp' CXXLIBS='${CXXLIBS} -Wl,--export-dynamic -Wl,-e,mexFunction -shared -lgomp' knn.cpp mexutil.cpp nn.cpp votemex.cpp patch.cpp vecnn.cpp simnn.cpp allegro_emu.cpp -output nnmex -output votemex
elseif strcmp('MACI64',computer)
    mex CXXOPTIMFLAGS='-O6 -w -s -ffast-math -fomit-frame-pointer -fstrength-reduce -fopenmp -msse2 -funroll-loops -fPIC' CXXFLAGS='-DNDEBUG -DUNIX_MODE -DMEXMODE -fopenmp' CXXLIBS='${CXXLIBS} -Wl,--export-dynamic -Wl,-e,mexFunction -shared -lgomp' knn.cpp mexutil.cpp nn.cpp nnmex.cpp patch.cpp vecnn.cpp simnn.cpp allegro_emu.cpp -output nnmex -output nnmex
mex CXXOPTIMFLAGS='-O6 -w -s -ffast-math -fomit-frame-pointer -fstrength-reduce -fopenmp -msse2 -funroll-loops -fPIC' CXXFLAGS='-DNDEBUG -DUNIX_MODE -DMEXMODE -fopenmp' CXXLIBS='${CXXLIBS} -Wl,--export-dynamic -Wl,-e,mexFunction -shared -lgomp' knn.cpp mexutil.cpp nn.cpp votemex.cpp patch.cpp vecnn.cpp simnn.cpp allegro_emu.cpp -output nnmex -output votemex
end
cd ..