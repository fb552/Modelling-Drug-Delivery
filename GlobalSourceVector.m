function [GVsource] = GlobalSourceVector(Ne,mesh,GQ,order)
%Assembles the single local source elements into a global vector
% The local source elements are computed at each element in the finite
% element mesh. Based on their location they are then inserted in the 
% global souce vector of size mesh.ngn-by-1. 
%
% Input:
%  Ne : Number of elements
%  mesh : Finite element mesh
%  GQ : Gaussian Quadrature parameters
%  order : weather the basis functions is linear or quadratic
% Return:
%  GVsource : Global Source vector
%
%Francesco Berteau (fb552) - November 2023

    %initialize matrix with zeros
    GVsource = zeros(mesh.ngn, 1);    
    for eN = 1:Ne
        %local 2-by-1 vector for corresponding eN
        [LVsource] = LocalSourceElem(eN,mesh,GQ,order);

        %location in GVector
        i = order*eN - (order-1);
        %add local vectors into global vector
        GVsource(i:i+order,1) = GVsource(i:i+order,1) + LVsource;
    end

end