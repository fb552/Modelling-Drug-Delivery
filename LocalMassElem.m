function [LMmass] = LocalMassElem(eN,mesh,GQ,order)
% This function computes the local mesh.ngn-by-mesh.ngn element mass matrix 
% for any element (eN) in the finite element mesh.
% In addition choose between manual integration (GQ=0) or Gaussian
% Quadrature (GQ=1), and between linear or quadratic basis function
% (orderbae)
%
% Francesco Berteau (fb552) - November 2023

    %Jacobian for a given element
    J = mesh.elem(eN).J;
    %initialise mass matrix with zeros
    LMmass = zeros(order+1);  
    %initialise local element matrix with zeros
    Int00 = zeros(order+1);
    
    %if manual integration
    if GQ.Switch == '0'
        %value of matrix at m = 0 & n = 0
        Int00 = 2*J/3;
        
        %2-by-2 element mass matrix
        LMmass = [Int00 Int00/2; Int00/2 Int00];
    
    %if Gaussian Quadrature
    elseif GQ.Switch == '1'
        N = GQ.N; %Access order of solver value
        [GQ] = GaussianQuadrature(GQ); %Generates Gaussian quadrature
        %Loops over number of Gauss points (xpoint)
        for i = 1:N
            weight = GQ.gaussweights(i); %Access Gauss weights
            xpoint = GQ.gaussxpoints(i); %Access Gauss points
            %Returns the corresponding value of the basis function at the given
            %Gauss point
            [psi,~] = BasisFuncEval(order,xpoint);
            %Creates the local element mass matrix multiplying together by
            %matching Gauss weights and adds to integral value
            for m = 1:order+1
                for n = 1:order+1
                    Int00(m,n) = psi(m)*psi(n)*J;
                    LMmass(m,n) = LMmass(m,n) + weight*(Int00(m,n));
                end
            end
        end
    end
end