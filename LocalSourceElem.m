function [LVsource] = LocalSourceElem(eN,mesh,GQ,order)
%Calculates the Local m-by-1 Element vector for the source term
% This is done for any element eN in the finite element mesh.
%
% Input:
%  eN : Element number
%  mesh : Finite element mesh
%  GQ : Gaussian Quadrature parameters
%  order : weather the basis functions is linear or quadratic
% Return:
%  LVsource : Local element source vector
%
%Francesco Berteau (fb552) - November 2023
    
    %Jacobian for a given element
    J = mesh.elem(eN).J;
    %initialise source vector with zeros
    LVsource = zeros(order+1,1);
    %initialise local element vector with zeros
    Int00 = zeros(order+1);
    %source term coefficient values
    f = mesh.elem(eN).f;

    %if manual integration
    if GQ.switch == '0'
        %value of elements in matrix 
        Int00 = f*J;
        %local 2-by-1 vector for corresponding eN
        LVsource = [Int00;Int00];

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
            %local element source vector
            for m = 1:order+1
                Int00(m) = f*psi(m)*J;
                LVsource(m) = LVsource(m) + w*(Int00(m));
            end
        end
    end
end