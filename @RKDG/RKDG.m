classdef RKDG < handle
    %RKDG: Master class for the IMEX scheme
    % For examples see also the Tests folder    
    % 
    % (c) Florian Schneider (schneider@mathematik.uni-kl.de) 2016
    % See licence.txt for more information about usage and redistribution
    
    properties
        Case;
        Grid;
        MomentModel;
        SpatialOrder;
        
        
        dt;
        CFL;
        Frames;
        Frame_cnt;
        psi;
        t;
        t_cnt;
        t_Frame;
        t_Frame_ind;
        u;
        u_Frame;
        u_0_rec;
        u0_Frame;
        runtime;
        
        CollisionOperator;
        
        save_flag;
        filename;
        restartCode;
        destroytime;
        
        NumThreads; %For parallel processing in OMP
        
        Statistics;
        
        Path
        Filename
    end
    
    methods
        
        function obj = RKDG(varargin)
            
            p = inputParser;
            p.addParameter('SpatialOrder',1);
            p.addParameter('MomentModel',[]);
            p.addParameter('Grid',[]);
            p.addParameter('Case',[]);
            p.addParameter('Frames',20);            
            p.addParameter('restartCode','');
            p.addParameter('destroytime',inf);
            p.addParameter('CollisionOperator','LaplaceBeltrami');
            p.KeepUnmatched = false;
            
            p.addParameter('Path',[]);
            p.addParameter('Filename',[]);
            p.addParameter('RKDGData',[]);
            p.addParameter('save_flag',false)
            p.addParameter('load_flag',false)
            
            p.parse(varargin{:});
            
            obj.Path = p.Results.Path;
            obj.Filename = p.Results.Filename;
            
            if ~p.Results.load_flag || (isempty(p.Results.Filename) && isempty(p.Results.RKDGData))
                
                obj.SpatialOrder = p.Results.SpatialOrder;
                obj.CollisionOperator = p.Results.CollisionOperator;
                obj.restartCode = p.Results.restartCode;
                obj.destroytime = p.Results.destroytime;
                obj.MomentModel = p.Results.MomentModel;
                
                
                if ~isempty(p.Results.Filename) && p.Results.save_flag
                     obj.save_flag = true;
                     if not(exist(p.Results.Path,'dir')==7)
                         mkdir(p.Results.Path)
                     end
                     if (exist([p.Results.Path p.Results.Filename],'file')==2)
                         warning('File already exists, will be overwritten');
                     end
                else
                    obj.save_flag = false;
                end
                
                
                obj.Case = p.Results.Case;
                obj.Frames = p.Results.Frames;
                obj.Grid= p.Results.Grid;

                if ~isempty(p.Results.Case)
                    obj.CalculateCFL();
                    obj.Frames = min(obj.Frames,length(obj.t));
                    fprintf('Calculate Initial conditions... ');
                    if isIsotropic(obj.Case.InitialCondition)
                        obj.psi = repmat(obj.Grid.projectToBasis(obj.Case.InitialCondition,obj.Grid.InteriorCells),1,1,obj.MomentModel.nq);
                    else
                        obj.psi = zeros(length(obj.Grid.InteriorCells),size(obj.Grid.Basis,3),obj.MomentModel.nq);
                        for i=1:size(obj.MomentModel.QuadraturePoints,2)
                            obj.psi(:,:,i) = obj.Grid.projectToBasis(@(x,y) obj.Case.InitialCondition(obj.MomentModel.QuadraturePoints(1,i),obj.MomentModel.QuadraturePoints(2,i),x,y),obj.Grid.InteriorCells); %Very slow...
                        end
                    end
                    obj.u = obj.MomentModel.ProjectToBasis(obj.psi);
                    fprintf('Done!\n');
                    
                    
                    
                    obj.Frame_cnt = 1;
                    obj.t_Frame_ind = unique([1:ceil(length(obj.t)/(obj.Frames-1)):length(obj.t),length(obj.t)]);
                    obj.t_Frame = obj.t(obj.t_Frame_ind);
                    obj.Frames = length(obj.t_Frame_ind);
                    obj.u_Frame = repmat(obj.u,[1,1,1,obj.Frames]);
                    obj.u_Frame(:,:,:,2:end) = NaN;
                    
                    obj.Frame_cnt = 1;
                    
                    obj.runtime = 0;
                    obj.NumThreads = 3; %Needs to be changed for every machine
                end
                
            else
                if isempty(p.Results.RKDGData)
                    obj.load(p.Results.Path,p.Results.Filename);
                else
                    obj.load(p.Results.RKDGData);
                end
            end
        end
        
        function CalculateCFL(obj)
            dx = obj.Grid.dx;
            dy = obj.Grid.dy;
            obj.CFL = 0.99;
            dteq = dx*dy/(dx+dy);
            obj.dt = obj.CFL*min(dteq(:));
            nsteps = ceil(obj.Case.t_final/obj.dt);
            obj.dt = obj.Case.t_final/nsteps;
            obj.t = 0:obj.dt:obj.Case.t_final;
            obj.t_cnt = 0;
        end
        
        function [U,J] = Collision(obj,psi,U)
            switch obj.CollisionOperator
                case 'LaplaceBeltrami'
                    [U,J] = obj.MomentModel.LaplaceBeltrami(psi,U);
                case 'isotropic'
                    [U,J] = obj.MomentModel.NeutronKernel(psi,U);
                otherwise
                    error('Operator not implemented')
            end
            
        end
        
        function RKDGData = save(obj,varargin)
            warning('off','MATLAB:structOnObject');
            RKDGData = struct(obj);
            RKDGData.Name = class(obj);
            warning('on','MATLAB:structOnObject');
            if nargin==3
                path = varargin{1}; 
                filename = varargin{2};  %#ok<*PROP>
                if strcmp(filename(end-3:end),'.mat')
                    filename = filename(1:end-4);
                end
                RKDGData.Case = obj.Case.save(path,[filename '_Case.mat']);
                RKDGData.MomentModel = obj.MomentModel.save(path,[filename '_MomentModel.mat']);
                RKDGData.Grid = obj.Grid.save(path,[filename '_Grid.mat']);
                save([path filename '.mat'],'RKDGData')
            else
                RKDGData.Case = obj.Case.save();
                RKDGData.MomentModel = obj.MomentModel.save();
                RKDGData.Grid = obj.Grid.save();            
            end
        end
        
        function load(obj,varargin)
            if nargin==2
                RKDGData = varargin{1}; 
            else
                path = varargin{1}; 
                filename = varargin{2}; 
                load([path filename]); %gives TestCaseData
            end
            if strcmp(class(obj),RKDGData.Name)
                RKDGData = rmfield(RKDGData,'Name');
                F = fieldnames(RKDGData);
                %try
                for i=1:length(F)
                    obj.(F{i}) = RKDGData.(F{i});
                end
                
                eval(sprintf('obj.Case = %s(''TestCaseData'',RKDGData.Case);',RKDGData.Case.Name));
                eval(sprintf('obj.MomentModel = %s(''ModelData'',RKDGData.MomentModel);',RKDGData.MomentModel.Name));
                eval(sprintf('obj.Grid = %s(''GridData'',RKDGData.Grid);',RKDGData.Grid.Name));
            else
                error(sprintf('It appears that you wanted to load the data for class "%s" into class "%s"\n',RKDGData.Name,class(obj))); %#ok<SPERR>
            end
            
        end
        
        function S = saveobj(obj)
            S = obj.save();
        end
        
        
        %% functions in own .m files
        du = DifferentialOperator(obj,U,t);
        U = DifferentialOperatorRHS(obj,U,t);
        plot(obj,varargin);
        run(obj);
        v = evaluateDensity(obj,x,y,t);
    end
    
    methods (Static)
        function obj = loadobj(S)
            if isstruct(S)
                obj = RKDG('load_flag',true,'RKDGData',S);
            end
        end
    end
    
end

