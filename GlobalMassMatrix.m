function [GMmass] = GlobalMassMatrix(Ne,mesh,GQ,order)
%This function uses a FOR loop over all elements in the finite element mesh
%to calculate the local element mass matrices. Lastly, it places the local
%element matrices in the correct location in the global mass matrix.
%
%Input arguments:
%Ne - Number of elements
%mesh – Finite element mesh
%GQ - Switch between manual integration and Gaussian quadrature rule
%order - Switch between linear and quadratic basis functions
%
%Return arguments:
%GlobalMassMatrix - Global Mass matrix

    %initialize global mass matrix with zeros
    GMmass = zeros(mesh.ngn);
    
    %calculate and assembly the global mass matrix
    for eN = 1:Ne
    
        %local matrix for mass element
        [LMmass] = LocalMassElem(eN,mesh,GQ,order);
    
        %location in GMatrix
        i = order*eN - (order-1);
    
        GMmass(i:i+order,i:i+order) = GMmass(i:i+order,i:i+order) + LMmass;
    end
end
