function [GMstiffness] = GlobalStiffnessMatrix(Ne,mesh,GQ,order)
%This function uses a FOR loop over all elements in the finite element mesh
%to calculate the local element diffusion and reaction matrices. Lastly, it
%places the local element matrices in the correct location in the global
%stiffness matrix.
%
%Input arguments:
%Ne - Number of elements
%mesh - Finite element mesh
%GQ - Switch between manual integration and Gaussian quadrature rule
%order - Switch between linear and quadratic basis functions
%
%Return arguments:
%GMstiffness - Global Stiffness matrix

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