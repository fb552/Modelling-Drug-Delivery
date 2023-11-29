function [RHS,GlobalMatrix] = DirichletBoundary(boundary,mesh,GlobalMatrix,RHS)
%This function enforces Dirichlet boundary conditions by setting matrix
%rows to zero that correspond to Dirichlet boundary nodes, setting those
%rows of the RHS vector equal to the boundary condition value and setting
%the diagonal entry of the BC matrix rows to 1.
%
%Input arguments:
%Boundary - Data structure containing all the boundary conditions
%mesh - Finite element mesh
%GlobalMatrix - Global Matrix without Dirichlet boundary conditions
%RHS - RHS vector without boundary conditions
%
%Return arguments:
%RHS - New RHS vector with Dirichlet boundary conditions
%GlobalMatrix -New Global Matrix with Dirichlet boundary conditions
    
    %initialise boundary conditions vector with zeros
    DB = zeros(mesh.ngn,1);
    %Dirichlet boundary conditions
    DB{1} = boundary.DirichletL;        %lower BC
    DB{mesh.ngn} = boundary.DirichletR; %upper BC

    %upper Dirichlet BC
    GlobalMatrix(end, :) = 0;
    GlobalMatrix(end) = 1;
    RHS(end) = DB(mesh.ngn);

    %lower Dirichlet BC
    GlobalMatrix(1, :) = 0;
    GlobalMatrix(1, 1) = 1;
    RHS(1) = DB(1);

%     %Loop to apply upper and lower boundary conditions
%     for i = [1 mesh.ngn]
%         % upper and lower DBC
%         if isnumeric(DB{i})
%             GlobalMatrix(i,:) = 0;
%             GlobalMatrix(i,i) = 1;
%             RHS(i) = DB{i};
%         end
%     end      
end