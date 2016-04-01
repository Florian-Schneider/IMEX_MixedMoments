classdef Pn < MomentModel
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function obj = Pn(varargin)
            
            %Inputparser
            p = inputParser;
            p.KeepUnmatched=true;
            p.addParameter('QuadraturePoints',[]);
            p.addParameter('QuadratureWeights',[]);
            p.parse(varargin{:});
            
            Q = LebedevQuadrature(2*p.Unmatched.MomentOrder+1);
            quads.QuadraturePoints = Q(1:2,:);
            quads.QuadratureWeights = Q(3,:);
            fieldNames = fieldnames(p.Unmatched);
            
            for i = 1:size(fieldNames,1)
                quads.(fieldNames{i}) = p.Unmatched.(fieldNames{i});
            end
            obj = obj@MomentModel(quads);                        
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
            psi = obj.BasisFunctionsAtQuadratureSqueezed'*u; %Basis is orthonormal -> Fourier series
            out = [];
        end
        
        function UC = NeutronKernel(~,~,U)
            UC = -U;
            UC(1,:,:,:) = 0;
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

