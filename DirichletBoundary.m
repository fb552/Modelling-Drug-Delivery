function [RHS,GlobalMatrix] = DirichletBoundary(boundary,RHS,GlobalMatrix)
%Applies Dirichlet Boundary Conditions to transient diffusion equation
% Sets the global matrix rows corresponding to boundary nodes to zero and
% the RHS vector to the boundaries values.
%
% Input:
% Boundary : Data structure containing all the boundary conditions
% RHS - RHS vector
% GlobalMatrix : Global Matrix
%
% Return:
% RHS : RHS vector with DBC applied
% GlobalMatrix : Global Matrix with DBC applied
%
% Francesco Berteau (fb552) - November 2023
    
    %upper Dirichlet BC
    GlobalMatrix(end, :) = 0;
    GlobalMatrix(end) = 1;
    RHS(end) = boundary.DirichletU;

    %lower Dirichlet BC
    GlobalMatrix(1, :) = 0;
    GlobalMatrix(1, 1) = 1;
    RHS(1) = boundary.DirichletL;
     
end
