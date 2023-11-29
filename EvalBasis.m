function [psi, dpsidxigrad] = EvalBasis(order,xipts)
%Evaluates the basis functions used in the model
% The basis functions are evaluated by calculating their values and gradients
% for a given local node (lnid) at a xipts.
%
% Input:
%  order : weather the basis functions is linear or quadratic
%  xipts : x points between [-1,1]
% Return:
%  psi : the value of the basis function
%  nvecnext : the gradient of the basis function
%
%Francesco Berteau (fb552) - November 2023
    
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