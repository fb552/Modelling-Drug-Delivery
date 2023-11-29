function [LVsource] = LocalSourceElem(eN,mesh,GQ,order)
%This function calculates the local 2-by-1 or 3-by-1 element vector for the
%source term depending if implementing manual integration or gaussian quadrature
%rule, for any element, eN, in the finite element mesh. This function also
%allows the user to select between linear and quadratic basis functions.
%
%Input arguments:
%eN - Element ID
%mesh - Finite element mesh
%GQ - Switch between manual integration and Gaussian quadrature rule
%order - Switch between linear and quadratic basis functions
%
%Return arguments:
%LVsource - Local 2-by-1 or 3-by-1 element source vector
    
    %Jacobian for a given element
    J = mesh.elem(eN).J;
    %initialise source vector with zeros
    LVsource = zeros(order+1,1);
    %initialise local element vector with zeros
    Int00 = zeros(order+1);
    %source term coefficient values
    f = mesh.elem(eN).f;

    %if manual integration
    if GQ.Switch == '0'
        %value of elements in matrix 
        Int00 = f*J;
        %local 2-by-1 vector for corresponding eN
        LVsource = [Int00;Int00];

    %if gaussian quadrature
    elseif GQ.Switch == '1'
        %order of solver
        N = GQ.N; 
        %get Gaussian points and weights
        [GQ] = GaussianQuadrature(GQ);

        for i = 1:N
            w = GQ.gw(i);          %Gauss weights
            xipts = GQ.xipts(i);   %Gauss points
            %basis function value at given point
            [psi,~] = EvalBasis(order,xipts);
            %local element source vector
            for m = 1:order+1
                Int00(m) = f*psi(m)*J;
                LVsource(m) = LVsource(m) + weight*(Int00(m));
            end
        end
    end
end