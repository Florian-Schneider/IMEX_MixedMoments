classdef QK1 < QMn
    
    
    properties
    end
    
    methods
        
        function obj = QK1(varargin)
            %Inputparser
            p = inputParser;
            p.KeepUnmatched=true;
            p.parse(varargin{:});
            Fields = p.Unmatched;
            Fields.MomentOrder = 1; %Force moment order.
            obj = obj@QMn(Fields);
            obj.positivityFlag = 1;
            NN = obj.NumberMoments/4;
            sizecap = @(s) s(2:end);
            obj.Density = @(U) reshape(U(1,:)+U(NN+1,:)+U(2*NN+1,:)+U(3*NN+1,:),[1,sizecap(size(U))]);
            [ obj.BasisFunctions,obj.BasisFunctionsAtQuadrature,obj.BasisFunctionsAtQuadratureSqueezed] = obj.generateBasisFunctions;
        end
        
        function u0 = ProjectToDensity(obj,U,c)
            %Moments are stored in the c-th component of U. Function
            %projects to u0
            Up = permute(U,[c,1:c-1,c+1:length(size(U))]);
            u0 = (Up(1,:)+Up(4,:)+Up(7,:)+Up(10,:))'/obj.BasisFunctionsAtQuadratureSqueezed(1);
            s = size(Up);
            u0 = reshape(u0,s(2:end));
            
        end
        
        function [psi] = ReconstructDistribution(obj,u)
            error('Kershaw closure...')
        end
        
               
        function  [FXP,FXM,FYP,FYM,Uc,Jc] = Closure(obj,U,R,scatteringtrue)
            tol = 1e-10;
            
            xi = -(2*pi-3*pi*sqrt((2))+4)/(3*pi*(sqrt((2))-1));
            a = 1;
            b = (-(2)^(1/2)*(a - xi));

            Epp = U(1,:,:,:);
            Fxpp = U(2,:,:,:);
            Fypp = U(3,:,:,:);
            Emp = U(4,:,:,:);
            Fxmp = U(5,:,:,:);
            Fymp = U(6,:,:,:);
            Emm = U(7,:,:,:);
            Fxmm = U(8,:,:,:);
            Fymm = U(9,:,:,:);
            Epm = U(10,:,:,:);
            Fxpm = U(11,:,:,:);
            Fypm = U(12,:,:,:);
            
            fppx = Fxpp./Epp;
            fppy = Fypp./Epp;
            fmpx = Fxmp./Emp;
            fmpy = Fymp./Emp;
            fmmx = Fxmm./Emm;
            fmmy = Fymm./Emm;
            fpmx = Fxpm./Epm;
            fpmy = Fypm./Epm;
            
            fppx(fppx<0 | Epp<=0) = 0;
            fppy(fppy<0 | Epp<=0) = 0;
            n = sqrt(fppx.^2+fppy.^2);
            fppx(n>=1) = fppx(n>=1)./n(n>=1);
            fppy(n>=1) = fppy(n>=1)./n(n>=1);
            
            fmpx(fmpx>0 | Emp<=0) = 0;
            fmpy(fmpy<0 | Emp<=0) = 0;
            n = sqrt(fmpx.^2+fmpy.^2);
            fmpx(n>=1) = fmpx(n>=1)./n(n>=1);
            fmpy(n>=1) = fmpy(n>=1)./n(n>=1);
            
            fmmx(fmmx>0 | Emm<=0) = 0;
            fmmy(fmmy>0 | Emm<=0) = 0;
            n = sqrt(fmmx.^2+fmmy.^2);
            fmmx(n>=1) = fmmx(n>=1)./n(n>=1);
            fmmy(n>=1) = fmmy(n>=1)./n(n>=1);
            
            fpmx(fpmx<0 | Epm<=0) = 0;
            fpmy(fpmy>0 | Epm<=0) = 0;
            n = sqrt(fpmx.^2+fpmy.^2);
            fpmx(n>=1) = fpmx(n>=1)./n(n>=1);
            fpmy(n>=1) = fpmy(n>=1)./n(n>=1);
            
            Fppx = Fxpp;
            Fppy = Fypp;
            Fmpx = Fxmp;
            Fmpy = Fymp;
            Fmmx = Fxmm;
            Fmmy = Fymm;
            Fpmx = Fxpm;
            Fpmy = Fypm;
            
            a1pp = (fppx.*fppy)*(16/(3*pi) - 8/3).*(fppx.^2+fppy.^2 - 1);
            n = sqrt(fppx.^2+fppy.^2);
            zeta = a+b*n;
            a2pp = -a1pp + zeta.*n.^2+(1-zeta).*n;
            nxx = fppx.*fppx./(fppx.^2+fppy.^2);
            nxy = fppx.*fppy./(fppx.^2+fppy.^2);
            nyy = fppy.*fppy./(fppx.^2+fppy.^2);
            nxx((fppx.^2+fppy.^2).^2<=tol) = 0;
            nxy((fppx.^2+fppy.^2).^2<=tol) = 0;
            nyy((fppx.^2+fppy.^2).^2<=tol) = 0;
            Dppxx = a1pp+a2pp.*nxx;
            Dppxy = a2pp.*nxy;
            Dppyy = a1pp+a2pp.*nyy;
            
            a1mp = -(fmpx.*fmpy)*(16/(3*pi) - 8/3).*(fmpx.^2+fmpy.^2 - 1);
            n = sqrt(fmpx.^2+fmpy.^2);
            zeta = a+b*n;
            a2mp = -a1mp + zeta.*n.^2+(1-zeta).*n;

            nxx = fmpx.*fmpx./(fmpx.^2+fmpy.^2);
            nxy = fmpx.*fmpy./(fmpx.^2+fmpy.^2);
            nyy = fmpy.*fmpy./(fmpx.^2+fmpy.^2);
            nxx((fmpx.^2+fmpy.^2).^2<=tol) = 0;
            nxy((fmpx.^2+fmpy.^2).^2<=tol) = 0;
            nyy((fmpx.^2+fmpy.^2).^2<=tol) = 0;
            Dmpxx = a1mp+a2mp.*nxx;
            Dmpxy = a2mp.*nxy;
            Dmpyy = a1mp+a2mp.*nyy;
            
            a1mm = (fmmx.*fmmy)*(16/(3*pi) - 8/3).*(fmmx.^2+fmmy.^2 - 1);
            n = sqrt(fmmx.^2+fmmy.^2);
            zeta = a+b*n;
            a2mm = -a1mm + zeta.*n.^2+(1-zeta).*n;

            nxx = fmmx.*fmmx./(fmmx.^2+fmmy.^2);
            nxy = fmmx.*fmmy./(fmmx.^2+fmmy.^2);
            nyy = fmmy.*fmmy./(fmmx.^2+fmmy.^2);
            nxx((fmmx.^2+fmmy.^2).^2<=tol) = 0;
            nxy((fmmx.^2+fmmy.^2).^2<=tol) = 0;
            nyy((fmmx.^2+fmmy.^2).^2<=tol) = 0;
            Dmmxx = a1mm+a2mm.*nxx;
            Dmmxy = a2mm.*nxy;
            Dmmyy = a1mm+a2mm.*nyy;
            
            a1pm = -(fpmx.*fpmy)*(16/(3*pi) - 8/3).*(fpmx.^2+fpmy.^2 - 1);
            n = sqrt(fpmx.^2+fpmy.^2);
            zeta = a+b*n;
            a2pm = -a1pm + zeta.*n.^2+(1-zeta).*n;

            nxx = fpmx.*fpmx./(fpmx.^2+fpmy.^2);
            nxy = fpmx.*fpmy./(fpmx.^2+fpmy.^2);
            nyy = fpmy.*fpmy./(fpmx.^2+fpmy.^2);
            nxx((fpmx.^2+fpmy.^2).^2<=tol) = 0;
            nxy((fpmx.^2+fpmy.^2).^2<=tol) = 0;
            nyy((fpmx.^2+fpmy.^2).^2<=tol) = 0;
            Dpmxx = a1pm+a2pm.*nxx;
            Dpmxy = a2pm.*nxy;
            Dpmyy = a1pm+a2pm.*nyy;
            
            Pppxx = Epp.*Dppxx;
            Pppxy = Epp.*Dppxy;
            Pppyy = Epp.*Dppyy;
            Pmpxx = Emp.*Dmpxx;
            Pmpxy = Emp.*Dmpxy;
            Pmpyy = Emp.*Dmpyy;
            Pmmxx = Emm.*Dmmxx;
            Pmmxy = Emm.*Dmmxy;
            Pmmyy = Emm.*Dmmyy;
            Ppmxx = Epm.*Dpmxx;
            Ppmxy = Epm.*Dpmxy;
            Ppmyy = Epm.*Dpmyy;
            
            FXP = cat(1,Fppx,Pppxx,Pppxy,0*Fmpx,0*Pmpxx,0*Pmpxy,0*Fmmx,0*Pmmxx,0*Pmmxy,Fpmx,Ppmxx,Ppmxy);
            FYP = cat(1,Fppy,Pppxy,Pppyy,Fmpy,Pmpxy,Pmpyy,0*Fmmy,0*Pmmxy,0*Pmmyy,0*Fpmy,0*Ppmxy,0*Ppmyy);
            FXM = cat(1,0*Fppx,0*Pppxx,0*Pppxy,Fmpx,Pmpxx,Pmpxy,Fmmx,Pmmxx,Pmmxy,0*Fpmx,0*Ppmxx,0*Ppmxy);
            FYM = cat(1,0*Fppy,0*Pppxy,0*Pppyy,0*Fmpy,0*Pmpxy,0*Pmpyy,Fmmy,Pmmxy,Pmmyy,Fpmy,Ppmxy,Ppmyy);
            
            if scatteringtrue
                [Uc,Jc] = R.Collision([],U);
                Uc = reshape(Uc,size(U));
            else
                Uc = 0*U;
                Jc = zeros(obj.NumberMoments,obj.NumberMoments,size(U,2));
            end
        end
        
    end
    
    
end

