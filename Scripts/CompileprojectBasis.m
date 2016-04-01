switch computer
    case 'PCWIN64'
        mex -v projectBasis.cpp -largeArrayDims COMPFLAGS="\$COMPFLAGS /openmp /DNDEBUG" OPTIMFLAGS="\$OPTIMFLAGS /openmp /DNDEBUG"
    case {'MACI64','GLNXA64'}
        mex -v -largeArrayDims projectBasis.cpp CXXFLAGS="\$CXXFLAGS -fopenmp -DNDEBUG" LDFLAGS="\$LDFLAGS -fopenmp -DNDEBUG"
end