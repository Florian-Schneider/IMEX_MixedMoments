%% Paths

p=genpath(cd);      %List of all subfolders
addpath(p,1);       %add them to the path


%% Compile mex files
if isempty(which('dualAdaptivePoly_mex'))
    CompileDualAdaptivePoly;
end
if isempty(which('projectBasis'))
    CompileProjectBasis;
end



clear p;   



clc
