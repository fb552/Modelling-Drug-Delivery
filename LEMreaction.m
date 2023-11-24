function ReactionMatrix = LEMreaction(lambda,eN,mesh)
% This function calculates the Local 2-by-2 Element Matrix for the linear
% reaction operator, for any element, eN, in the finite element mesh.
% The reaction coefficient lambda and the mesh.
    
    %Jacobian of given element
    J = mesh.elem(eN).J;

    %elements of local element matrix
    SqMat = (lambda*J)/3;

    %local 2x2 Element Matrix for reaction
    ReactionMatrix = [2*SqMat, SqMat; SqMat, 2*SqMat];
end