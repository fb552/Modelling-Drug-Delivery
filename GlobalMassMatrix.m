function [GMmass] = GlobalMassMatrix(Ne,mesh,GQ,order)
%Assembles the single local mass elements into a global matrix
% The local mass elements are computed at each element in the finite
% element mesh. Based on their location they are then inserted in the 
% global mass matrix of size mesh.ngn-by-mesh.ngn. 
%
% Input:
%  Ne : Number of elements
%  mesh : Finite element mesh
%  GQ : Gaussian Quadrature parameters
%  order : weather the basis functions is linear or quadratic
% Return:
%  GMmass : Global Mass matrix
%
%Francesco Berteau (fb552) - November 2023

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