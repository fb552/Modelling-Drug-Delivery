function [NB,NBnext] = NeumannBoundary(boundary,mesh)
%Computes the Neumann Boundary Conditions
% This is done according to the lower and upper boundary values contained
% within boundary
%
% Input:
%  Boundary : Data structure containing all the boundary conditions
%  mesh : Finite element mesh
% Return:
%  NB : current Neumann boundary conditions vector
%  NBnext : next Neumann boundary conditions vector
%
%Francesco Berteau (fb552) - November 2023

    %initialise boundary conditions vector with zeros
    NB = zeros(mesh.ngn,1);

    %Neumann boundaries implementation
    if isnumeric(boundary.NeumannL)
        %lower Neumann BC
        NB(1) = -boundary.NeumannL;
    elseif isnumeric(boundary.NeumannU)
        %upper Neumann BC
        NB(mesh.ngn) = boundary.NeumannU;
    end
    %store the next BC
    NBnext = NB;
end