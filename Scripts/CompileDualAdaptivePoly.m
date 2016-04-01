function CompileDualAdaptivePoly()

switch computer
    case 'PCWIN64'
        blaslib = fullfile(matlabroot, ...
            'extern', 'lib', 'win64', 'microsoft', 'libmwblas.lib');
        lapacklib = fullfile(matlabroot, ...
            'extern', 'lib', 'win64', 'microsoft', 'libmwlapack.lib');
        % mex('-v','-g', '-largeArrayDims', 'dualAdaptivePoly_mex.c', lapacklib,blaslib) %Debugger
        mex('-v','-O','-largeArrayDims', 'dualAdaptivePoly_mex.c','-Ifdlibm','COMPFLAGS ="$COMPFLAGS -D_IEEE_LIBM -D__GW32__ -D_LARGEFILE_SOURCE=1 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -Wall -O3 -fms-extensions -mms-bitfields  -fno-exceptions -fomit-frame-pointer -march=i386 -mcpu=i686"','-D__LITTLE_ENDIAN=1',lapacklib,blaslib) %Optimized
        % ['mex CXX="gcc" LD="gcc" CXXOPTIMFLAGS="$CXXOPTIMFLAGS -fopenmp -O2 -DNDEBUG" LDCXXOPTIMFLAGS="$LDCXXOPTIMFLAGS -fopenmp -O2" -lgomp ' name ' -v -lstdc++'];
        
    case {'MACI64','GLNXA64'}
        mex('-v','-O','-largeArrayDims','dualAdaptivePoly_mex.c','-Ifdlibm','COMPFLAGS ="$COMPFLAGS -D_IEEE_LIBM -D__GW32__ -D_LARGEFILE_SOURCE=1 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64"','-D__LITTLE_ENDIAN=1', '-lmwlapack','-lmwblas') %Optimized
        
        
        
end
end