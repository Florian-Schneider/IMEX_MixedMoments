classdef MK1Table < MK1
    %MK1Table: Mixed-moment Kershaw closure of first order with tabulated
    %QM1 Laplace-Beltrami operator
    
    properties
        Table_pp;
        Table_mp;
        Table_mm;
        Table_pm;
    end
    
    methods
        
        function obj = MK1Table(varargin)
            obj = obj@MK1(varargin{:});            
            obj.generateQuarterSpaceTables();
        end
        
        
        function [LB,JLB] = LaplaceBeltrami(obj,~,U)
            tol = 1e-10;
            JLB = [];
            
            E = U(1,:,:,:);
            Fxp = U(2,:,:,:);
            Fxm = U(3,:,:,:);
            Fyp = U(4,:,:,:);
            Fym = U(5,:,:,:);
            
            fxp = Fxp./E;
            fxm = Fxm./E;
            fyp = Fyp./E;
            fym = Fym./E;
            
            lambda1 = fyp./(fyp-fym);
            lambda1(fyp-fym<tol) = 0;
            lambda1 = max(min(lambda1,1),0);
            lambda2 = fxp./(fxp-fxm);
            lambda2(fxp-fxm<tol) = 0;
            lambda2 = max(min(lambda2,1),0);
            
            fppx = (fxp-fxm);
            fppy = (fyp-fym);
            fpmx = (fxp-fxm);
            fpmy = -(fyp-fym);
            fmpx = -(fxp-fxm);
            fmpy = (fyp-fym);
            fmmx = -(fxp-fxm);
            fmmy = -(fyp-fym);
            
            int1pp = obj.Table_pp.Interp_intlow(fppx,fppy);
            int1pm = obj.Table_pm.Interp_intup(fpmx,fpmy);
            
            int2pp = obj.Table_pp.Interp_intup(fppx,fppy);
            int2mp = obj.Table_mp.Interp_intlow(fmpx,fmpy);
            
            int3mm = obj.Table_mm.Interp_intlow(fmmx,fmmy);
            int3mp = obj.Table_mp.Interp_intup(fmpx,fmpy);
            
            int4mm = obj.Table_mm.Interp_intup(fmmx,fmmy);
            int4pm = obj.Table_pm.Interp_intlow(fpmx,fpmy);
            
            int1 = E/2.*(lambda1.*lambda2.*(int1pp)+(1-lambda1).*lambda2.*int1pm); %E/2 since tables are normalized with respect to Sij instead of halfspaces
            int3 = E/2.*(lambda1.*(1-lambda2).*int3mp + (1-lambda1).*(1-lambda2).*int3mm);
            int2 = E/2.*(lambda1.*lambda2.*int2pp + lambda1.*(1-lambda2).*int2mp);
            int4 = E/2.*((1-lambda1).*lambda2.*int4pm + (1-lambda1).*(1-lambda2).*int4mm);
            
            LB = 0*U; %initializes u0
            lmax = obj.MomentOrder;
            l=1;
            LB(l+1,:,:,:) = -l*(l+1)*U(l+1,:,:,:) + int2+int4; %Sxp
            LB(lmax+l+1,:,:,:) = -l*(l+1)*U(lmax+l+1,:,:,:)-int2-int4; %Sxm
            LB(2*lmax+l+1,:,:,:) = -l*(l+1)*U(2*lmax+l+1,:,:,:)+int1+int3; %Syp
            LB(3*lmax+l+1,:,:,:) = -l*(l+1)*U(3*lmax+l+1,:,:,:)-int1-int3; %Sym
            
            LB = LB/2;
        end
        
        
        
        function generateQuarterSpaceTables(obj)
            
            Tables = load('MQK1Tables');
            obj.Table_pp = Tables.Tablepp;
            obj.Table_mp = Tables.Tablemp;
            obj.Table_mm = Tables.Tablemm;
            obj.Table_pm = Tables.Tablepm;
            
            obj.Table_pp.Interp_u2xx = scatteredInterpolant(obj.Table_pp.u1x(:),obj.Table_pp.u1y(:),obj.Table_pp.u2xx(:),'linear','nearest');
            obj.Table_pp.Interp_u2xy = scatteredInterpolant(obj.Table_pp.u1x(:),obj.Table_pp.u1y(:),obj.Table_pp.u2xy(:),'linear','nearest');
            obj.Table_pp.Interp_u2yy = scatteredInterpolant(obj.Table_pp.u1x(:),obj.Table_pp.u1y(:),obj.Table_pp.u2yy(:),'linear','nearest');
            obj.Table_pp.Interp_intlow = scatteredInterpolant(obj.Table_pp.u1x(:),obj.Table_pp.u1y(:),obj.Table_pp.int0(:),'linear','nearest');
            obj.Table_pp.Interp_intup = scatteredInterpolant(obj.Table_pp.u1x(:),obj.Table_pp.u1y(:),obj.Table_pp.int1(:),'linear','nearest');
            
            obj.Table_mp.Interp_u2xx = scatteredInterpolant(obj.Table_mp.u1x(:),obj.Table_mp.u1y(:),obj.Table_mp.u2xx(:),'linear','nearest');
            obj.Table_mp.Interp_u2xy = scatteredInterpolant(obj.Table_mp.u1x(:),obj.Table_mp.u1y(:),obj.Table_mp.u2xy(:),'linear','nearest');
            obj.Table_mp.Interp_u2yy = scatteredInterpolant(obj.Table_mp.u1x(:),obj.Table_mp.u1y(:),obj.Table_mp.u2yy(:),'linear','nearest');
            obj.Table_mp.Interp_intlow = scatteredInterpolant(obj.Table_mp.u1x(:),obj.Table_mp.u1y(:),obj.Table_mp.int0(:),'linear','nearest');
            obj.Table_mp.Interp_intup = scatteredInterpolant(obj.Table_mp.u1x(:),obj.Table_mp.u1y(:),obj.Table_mp.int1(:),'linear','nearest');
            
            obj.Table_mm.Interp_u2xx = scatteredInterpolant(obj.Table_mm.u1x(:),obj.Table_mm.u1y(:),obj.Table_mm.u2xx(:),'linear','nearest');
            obj.Table_mm.Interp_u2xy = scatteredInterpolant(obj.Table_mm.u1x(:),obj.Table_mm.u1y(:),obj.Table_mm.u2xy(:),'linear','nearest');
            obj.Table_mm.Interp_u2yy = scatteredInterpolant(obj.Table_mm.u1x(:),obj.Table_mm.u1y(:),obj.Table_mm.u2yy(:),'linear','nearest');
            obj.Table_mm.Interp_intlow = scatteredInterpolant(obj.Table_mm.u1x(:),obj.Table_mm.u1y(:),obj.Table_mm.int0(:),'linear','nearest');
            obj.Table_mm.Interp_intup = scatteredInterpolant(obj.Table_mm.u1x(:),obj.Table_mm.u1y(:),obj.Table_mm.int1(:),'linear','nearest');
            
            obj.Table_pm.Interp_u2xx = scatteredInterpolant(obj.Table_pm.u1x(:),obj.Table_pm.u1y(:),obj.Table_pm.u2xx(:),'linear','nearest');
            obj.Table_pm.Interp_u2xy = scatteredInterpolant(obj.Table_pm.u1x(:),obj.Table_pm.u1y(:),obj.Table_pm.u2xy(:),'linear','nearest');
            obj.Table_pm.Interp_u2yy = scatteredInterpolant(obj.Table_pm.u1x(:),obj.Table_pm.u1y(:),obj.Table_pm.u2yy(:),'linear','nearest');
            obj.Table_pm.Interp_intlow = scatteredInterpolant(obj.Table_pm.u1x(:),obj.Table_pm.u1y(:),obj.Table_pm.int0(:),'linear','nearest');
            obj.Table_pm.Interp_intup = scatteredInterpolant(obj.Table_pm.u1x(:),obj.Table_pm.u1y(:),obj.Table_pm.int1(:),'linear','nearest');
            
            
        end
        
    end
    
    
end

