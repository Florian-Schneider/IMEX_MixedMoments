%% Initialization
startup
nx = 100; ny = 100;
%% Model Generator
Model{1} = Mn('MomentOrder',1,'QuadratureOrder',113);
Model{2} = Mn('MomentOrder',2,'QuadratureOrder',113);
Model{3} = Mn('MomentOrder',3,'QuadratureOrder',113);
Model{4} = MMn('MomentOrder',1,'QuadratureOrder',30);
Model{5} = MMn('MomentOrder',2,'QuadratureOrder',30);
Model{6} = MMn('MomentOrder',3,'QuadratureOrder',30);
Model{7} = MK1();
Model{8} = MK1Table();

MinModConstant = 0; %Obsolete for first order
%% TestCase and Grid
SpatialOrder = 1;
T = Linesource('t_final',0.45);
G = T.generateGrid(nx,ny,SpatialOrder,Model,MinModConstant);

%% filename for saving and run

for i=1:length(Model)
    foldername = ['Results/Linesource/Linesource - LaplaceBeltrami - ' num2str(nx) ' - ' char(datetime('now','TimeZone','local','Format','y-dd-MM_HH-mm-ss')) '/'];    
    R{i} = RKDG('save_flag',true,'Path',foldername,'Filename','RKDGData.mat','MomentModel',Model{i},'Grid',G,'Case',T,'SpatialOrder',G.order,'CollisionOperator','isotropic');
    
    fprintf('Running Linesource test with model class %s of order %d\n',class(Model{i}),Model{i}.MomentOrder)
    R{i}.run;
    
    %To load the solution again use e.g.
    %R{i} = loadRKDG(foldername);
end

%% plotting
i = 7;
R{i}.plot('Scaling','nlin','Time',[0,0.45]) %Plot the i-th solution over the whole time period with a rescaling of z and colour-values in every frame
