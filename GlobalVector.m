function [sourceTerm] = GlobalVector(f, mesh)
% Calculate the global mesh.ngn-by-1 element matrix for the source vector (f), 
% for any element (eN) in the finite element mesh.

    %initialize matrix with zeros
    sourceTerm = zeros(mesh.ngn, 1);
    for eN = 1:mesh.ne

        %Jacobian for a given element
        J = mesh.elem(eN).J;

        %value of elements in matrix 
        SqMat = f*J;

        %local 2-by-1 vector for corresponding eN
        localTerm = [SqMat; SqMat];

        %add local matricies into global vector
        sourceTerm(eN:eN+1) = sourceTerm(eN:eN+1) + localTerm;
    end
end