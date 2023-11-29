function [psi, dpsidxigrad] = EvalBasis(order,xipts)
%This function evaluates the basis functions implemented in the model.
%Basis functions are evaluated by calculating their values and gradients
%for a given local node (lnid) at a xipts.
%
%Input arguments:
%order - Switch between linear and quadratic basis functions
%xipts - Coordinate x between [-1,1]
%
%Return arguments:
%psi - Value of the basis function
%nvecnext - Gradient of the basis function
    %initialise psi as zeros array, with size order + 1
    psi = zeros(order+1,1);
    dpsidxigrad = zeros(order+1,1);

    %linear or quadratic basis functions
    switch order
        %linear
        case 1 
            %function values
            for lnid = 0:1
                sign = (-1)^(lnid+1);
                psi(lnid+1) = 0.5*(1+((sign*xipts)));
            end
            %function gradients
            for lnid = 1:order+1
                sign1 = (-1)^(lnid);
                dpsidxigrad(lnid) = 0.5 * sign1;
            end
        %quadratic
        case 2 
            %function values
            psi(1) = xipts*(xipts - 1)*0.5;
            psi(2) = (1 - xipts^2);
            psi(3) = xipts*(xipts + 1)*0.5;
            %function gradients
            dpsidxigrad(order-1) = xipts-0.5;
            dpsidxigrad(order) = -2*xipts;
            dpsidxigrad(order+1) = xipts+0.5;
    end
end