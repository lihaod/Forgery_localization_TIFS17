%%
mex OPTIMFLAGS="/Ox /Oi /Oy /DNDEBUG /fp:fast /arch:SSE2 /DMEX_MODE /openmp" mexutil.cpp nn.cpp nnmex.cpp patch.cpp vecnn.cpp simnn.cpp allegro_emu.cpp knn.cpp -output nnmex_mex
mex OPTIMFLAGS="/Ox /Oi /Oy /DNDEBUG /fp:fast /arch:SSE2 /DMEX_MODE /openmp" mexutil.cpp nn.cpp votemex.cpp patch.cpp vecnn.cpp simnn.cpp allegro_emu.cpp knn.cpp -output votemex
