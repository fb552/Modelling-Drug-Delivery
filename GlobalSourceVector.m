function [GVsource] = GlobalSourceVector(Ne,mesh,GQ,order)
% Calculate the global mesh.ngn-by-1 element matrix for the source vector (f), 
% for any element (eN) in the finite element mesh.

    %initialize matrix with zeros
    GVsource = zeros(mesh.ngn, 1);    
    for eN = 1:Ne
        %local 2-by-1 vector for corresponding eN
        [LVsource] = LocalSourceElem(eN,mesh,GQ,order);

        %location in GVector
        i = orderbase*eN - (order-1);
        %add local vectors into global vector
        GVsurce(i:i+order,1) = GVsurce(i:i+order,1) + LVsource;
    end

end
