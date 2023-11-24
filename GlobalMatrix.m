function [Gmatrix] = GlobalMatrix(D, lambda, mesh)
% This function computes the global mesh.ngn-by-mesh.ngn element matrix for the
% Global Matrix for any element (eN) in the finite element mesh. Combining
% the local matrices for diffusion (D) and reaction (lambda) terms.

    %initialize matrix with zeros
    Gmatrix = zeros(mesh.ngn);
    for eN = 1:mesh.ne
        
        %local Laplace Element Matrix for diffusion
        LEMdiff = LEMdiffusion(D,eN,mesh); 

        %local Laplace Element Matrix for reaction
        LEMreac = LEMreaction(lambda,eN,mesh); 
        
        %difference between diffusion and reaction
        diff = LEMdiff - LEMreac; 
        
        %add local matricies into Global
        Gmatrix(eN:eN+1,eN:eN+1) =  Gmatrix(eN:eN+1,eN:eN+1) + diff;
    end 
end