function [NB,NBnext] = NeumannBoundary(boundary,mesh)
%This function modifies the Neumann Boundary conditions vector according to
%the values of Boundary.NeumannL and Boundary.NeumannR stored in the data
%structure Boundary.
%
%Input arguments:
%Boundary - Data structure containing all the boundary conditions
%mesh - Finite element mesh
%
%Return arguments:
%NBc - Current time point Neumann boundary conditions vector
%NBcnext - Next time point Neumann boundary conditions vector

    %initialise boundary conditions vector with zeros
    NB = zeros(mesh.ngn,1);

    %Neumann boundaries implementation
    if isnumeric(boundary.NeumannL)
        %lower Neumann BC
        NB(1) = -boundary.NeumannL;
    elseif isnumeric(boundary.NeumannR)
        %upper Neumann BC
        NB(mesh.ngn) = boundary.NeumannR;
    end
    %store the next BC
    NBnext = NB;
end