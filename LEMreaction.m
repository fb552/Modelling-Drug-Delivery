function [ReactionMatrix] = LEMreaction(eN,mesh,GQ,order)
%Calculates the Local m-by-n Element Matrix for the linear reaction operator
% This is done for any element eN in the finite element mesh.
%
% Input:
%  eN : Element number
%  mesh : Finite element mesh
%  GQ : Gaussian Quadrature parameters
%  order : weather the basis functions is linear or quadratic
% Return:
%  ReactionMatrix : Local reaction matrix
%
%Francesco Berteau (fb552) - November 2023
    
    %Jacobian of given element
    J = mesh.elem(eN).J;
    %initialise matrix with zeros
    ReactionMatrix = zeros(order+1);
    %initialise values of matrix with zeros
    Int00 = zeros(order +1);
    %reaction coefficient values
    lambda = mesh.elem(eN).lambda;
    
    %if manual integration
    if GQ.switch == '0'
        %elements of local element matrix
        Int00 = 2*(lambda*J)/3;
        %local 2x2 Element Matrix for reaction
        ReactionMatrix = [Int00 Int00/2 ; Int00/2 Int00];

    %if gaussian quadrature
    elseif GQ.switch == '1'
        %number of gauss points
        N = GQ.npts; 
        %get Gaussian points and weights
        [GQ] = GQscheme(GQ);

        for i = 1:N
            w = GQ.gw(i);          %Gauss weights
            xipts = GQ.xipts(i);   %Gauss points
            %basis function value at given point
            [psi,~] = EvalBasis(order,xipts);
            %local element reaction matrix
            for m = 1:order+1
                for n = 1:order+1
                    Int00(m,n) = lambda*psi(m)*psi(n)*J;
                    ReactionMatrix(m,n) = ReactionMatrix(m,n) + w*(Int00(m,n));
                end
            end
        end
    end
end