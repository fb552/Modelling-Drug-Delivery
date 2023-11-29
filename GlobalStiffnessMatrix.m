function [GMstiffness] = GlobalStiffnessMatrix(Ne,mesh,GQ,order)
%Assembles the single local diffusion and reaction elements into a global matrix
% The local diffusion and reaction elements are computed at each element 
% in the finite element mesh. Based on their location they are then 
% inserted in the global stiffness matrix of size mesh.ngn-by-mesh.ngn. 
%
% Input:
%  Ne : Number of elements
%  mesh : Finite element mesh
%  GQ : Gaussian Quadrature parameters
%  order : weather the basis functions is linear or quadratic
% Return:
%  GMstiffness : Global Stiffness matrix
%
%Francesco Berteau (fb552) - November 2023

    %initialize matrix with zeros
    GMstiffness = zeros(mesh.ngn);
    
    %calculate and assembly the global stiffness matrix
    for eN = 1:Ne
        %local Laplace Element Matrix for diffusion
        [LMdiffusion] = LEMdiffusion(eN, mesh,GQ,order);
        %local Laplace Element Matrix for reaction
        [LMreaction] = LEMreaction(eN, mesh,GQ,order);
    
        %difference between diffusion and reaction
        diff = LMdiffusion - LMreaction; 
        %location in GMatrix
        i = order*eN - (order-1);
        %add local matrices into global matrix
        GMstiffness(i:i+order,i:i+order) = GMstiffness(i:i+order,i:i+order) + diff;
    end
end