classdef Mn < MomentModel
    %Mn: Full-moment minimum-entropy models of arbitrary order in spherical
    %harmonics basis
    
    properties
        alpha;
    end
    
    methods
        
        function obj = Mn(varargin)
            
            %Inputparser
            p = inputParser;
            p.KeepUnmatched=true;
            p.parse(varargin{:});
            obj = obj@MomentModel(p.Unmatched);
            obj.positivityFlag = 1;
        end
        
        function [ BasisFunctions,BasisFunctionsAtQuadrature,BasisFunctionsAtQuadratureSqueezed] = generateBasisFunctions(obj)
            
            BasisFunctions = [];
            
            mu = obj.QuadraturePoints(1,:);
            phi = obj.QuadraturePoints(2,:);
            lmax = obj.MomentOrder;
            
            P = cell(lmax+1);
            
            for l=0:lmax
                L = legendre(l,mu);
                for m=1:l+1
                    P{m,l+1} = L(m,:);
                end
            end
            C = cell(lmax,1);
            S = C;
            [C{:}] = deal(cos(phi));
            [S{:}] = deal(sin(phi));
            
            for m=2:lmax
                C{m} = C{1}.*C{m-1}-S{1}.*S{m-1};
                S{m} = S{1}.*C{m-1}+C{1}.*S{m-1};
            end
            
            NumberMoments = 1/2*(lmax^2+3*lmax)+1;
            
            BasisFunctionsAtQuadrature = zeros(NumberMoments,1,size(mu,2));
            [l,m] = meshgrid(0:lmax);
            N = sqrt((2*l+1)/4/pi.*factorial(max(l-m,0))./factorial(l+m));
            N(m>l) = 0;
            cnt = 1;
            for l=0:lmax
                for m=-l:l
                    if mod(l+m,2)==0
                        if m<0
                            BasisFunctionsAtQuadrature(cnt,1,:) = sqrt(2)*S{-m}.*N(-m+1,l+1).*P{-m+1,l+1};
                        end
                        if m==0
                            BasisFunctionsAtQuadrature(cnt,1,:) = N(m+1,l+1).*P{m+1,l+1};
                        end
                        if m>0
                            BasisFunctionsAtQuadrature(cnt,1,:) = sqrt(2)*C{m}.*N(m+1,l+1).*P{m+1,l+1};
                        end
                        cnt = cnt+1;
                    end
                end
            end
            
            BasisFunctionsAtQuadratureSqueezed = squeeze(BasisFunctionsAtQuadrature);
        end
        function [psi,out] = ReconstructDistribution(obj,u)
            persistent beta BM
            
            beta_size = size(u);
            if isempty(beta) || ~isequal(size(beta),beta_size)
                beta = repmat([log(1/4/pi);zeros(size(u,1)-1,1)],[1,beta_size(2)]);
                BM = repmat(eye(size(u,1)),[1 1 beta_size(2)]);
            end
            
            
            OptimPara = obj.OptimPara;
            psi = zeros([size(obj.BasisFunctionsAtQuadrature,3),size(u,2)]);
            
            for i=1:size(u,2);
                U = u(:,i);
                [alphaF,beta(:,i),BM(:,:,i)] = dualAdaptivePoly(U, beta(:,i),BM(:,:,i),OptimPara);
                psi(:,i) = U(1)*exp(alphaF'*obj.BasisFunctionsAtQuadratureSqueezed);
            end
            
            out = [];
        end
        
        function [UC,J] = NeutronKernel(obj,~,U)
            UC = -U;
            UC(1,:,:,:) = 0;
            J = repmat(diag(-ones(obj.NumberMoments,1)),[1 1 size(UC,2)]);
            J(1,1,:) = 0;
        end
        
        function [LB,J] = LaplaceBeltrami(obj,~,U)
            LB = U;
            lmax = obj.MomentOrder;
            cnt = 1;
            J = zeros(lmax,lmax,size(U,2));
            for l=0:lmax
                for m=-l:l
                    if mod(l+m,2)==0
                        LB(cnt,:) = -l*(l+1)*LB(cnt,:);
                        J(cnt,cnt,:) = -l*(l+1);
                        cnt = cnt+1;
                    end
                    
                end
            end
        end
        
        
    end
    
end

