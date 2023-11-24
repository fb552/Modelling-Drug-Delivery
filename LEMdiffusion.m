function DiffusionMatrix = LEMdiffusion(D,eN,mesh)
% This function calculates the Local 2-by-2 Element Matrix for the linear
% diffusion operator, for any element, eN, in the finite element mesh.
% The diffusion coefficient D and the mesh.
    
    %Jacobian of given element
    J = mesh.elem(eN).J;

    %elements of local element matrix 2J = (x1-x0)
    SqMat = D/(2*J);

    %local 2x2 Element Matrix for diffusion
    DiffusionMatrix = [SqMat, -SqMat; -SqMat, SqMat];
end