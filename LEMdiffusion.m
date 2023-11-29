function [DiffusionMatrix] = LEMdiffusion(eN,mesh,GQ,order)
%Calculates the Local m-by-n Element Matrix for the linear diffusion operator
% This is done for any element eN in the finite element mesh.
%
% Input:
%  eN : Element number
%  mesh : Finite element mesh
%  GQ : Gaussian Quadrature parameters
%  order : weather the basis functions is linear or quadratic
% Return:
%  DiffusionMatrix : Local diffusion matrix
%
%Francesco Berteau (fb552) - November 2023
    
    %Jacobian for a given element
    J = mesh.elem(eN).J;
    %initialise matrix with zeros
    DiffusionMatrix = zeros(order+1);
    %initialise values of matrix with zeros
    Int00 = zeros(order +1);
    %diffusion coefficient values
    D = mesh.elem(eN).D;
     
    %if manual integration
    if GQ.switch == '0'
        %elements of local element matrix 2J = (x1-x0)
        Int00 = D/(2*J);
        %local 2x2 Element Matrix for diffusion
        DiffusionMatrix = [Int00 -Int00 ; -Int00 Int00];
    
    %if Gaussian Quadrature
    elseif GQ.switch == '1'
        %order of solver
        N = GQ.N; 
        %get Gaussian points and weights
        [GQ] = GQscheme(GQ);

        for i = 1:N
            w = GQ.gw(i);          %Gauss weights
            xipts = GQ.xipts(i);   %Gauss points
            %basis function at given point
            [~,dpsidxigrad] = EvalBasis(order,xipts);
            %local element diffusion matrix
            for m = 1:order+1
                for n = 1:order+1
                    Int00(m,n) = D*dpsidxigrad(m)*dpsidxigrad(n)/J;
                    DiffusionMatrix(m,n) = DiffusionMatrix(m,n) + w*(Int00(m,n));
                end
            end
        end
    end
end