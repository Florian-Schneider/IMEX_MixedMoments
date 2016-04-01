function [ R ] = loadRKDG( path )
%LOADRKDG Summary of this function goes here
%   Detailed explanation goes here
%s = load(path);
%R = RKDG.loadobj(s.RKDGData);

if ~(path(end)=='/') && ~(path(end)=='\')
    path = [path,'/'];
end
R = RKDG('load_flag',true,'Path',path,'Filename','RKDGData.mat');

end

