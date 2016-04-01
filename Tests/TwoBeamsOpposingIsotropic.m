%% Initialization
startup
nx = 20; ny = 20;
%% Model Generator
Model{1} = Mn('MomentOrder',1,'QuadratureOrder',113);
Model{2} = MK1();

MinModConstant = 0; %Obsolete for first order
%% TestCase and Grid
SpatialOrder = 1;
T = TwoBeamsOpposing('t_final',1.2);
G = T.generateGrid(nx,ny,SpatialOrder,Model,MinModConstant);

%% filename for saving and run

for i=2%1:length(Model)
    foldername = ['Results/TwoBeams/TwoBeams - isotropic - ' num2str(nx) ' - ' char(datetime('now','TimeZone','local','Format','y-dd-MM_HH-mm-ss')) '/'];    
    R{i} = RKDG('save_flag',true,'Path',foldername,'Filename','RKDGData.mat','MomentModel',Model{i},'Grid',G,'Case',T,'SpatialOrder',G.order,'CollisionOperator','isotropic'); %isotropic and Laplace-Beltrami are the same here since sigma_s = 0.
    
    fprintf('Running TwoBeamsOpposing test with model class %s of order %d\n',class(Model{i}),Model{i}.MomentOrder)
    R{i}.run;
    
    %To load the solution again use e.g.
    %R{i} = loadRKDG(foldername);
end

%% plotting
i = 2;
R{i}.plot('Scaling','lin','Time',[0,1.2]) %Plot the i-th solution over the whole time period with a fixed scale of z and colour-values
