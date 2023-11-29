function [mesh,c,GQ,time,GM] = TransientFEM(boundary,Xmin,Xmax,Ne,GQ,order,theta,time,parameters)
%This function solves the full transient form of the diffusion reaction
%equation.
%
%Input arguments:
%Boundary - Data structure containing all the boundary conditions
%xmin - Lower spatial boundary position
%xmax - Upper spatial boundary position
%Ne - Number of elements
%GQ - Switch between manual integration and Gaussian quadrature rule
%order - Switch between linear and quadratic basis functions
%theta - Selection of numerical scheme
%time - Data structure containing all the time-related data
%matparameters - Selection of material parameters
%
%Return arguments:
%mesh - Finite element mesh
%c - Numerical solution of the problem
%time - Time interval
%GM - Final Global Matrix with Dirichlet boundary conditions

    %finite element mesh between Xmin and Xmax with Ne number of elements
    mesh = OneDimLinearMeshGen(Xmin,Xmax,Ne,order,parameters);
    
    %get Gaussian points and weights
    [GQ] = GQscheme(GQ);

    %delta t
    dt = time.dt;
    
    %get global mass matrix
    [GMmass] = GlobalMassMatrix(Ne,mesh,GQ,order);
    
    %get global stiffness matrix
    [GMstiffness] = GlobalStiffnessMatrix(Ne,mesh,GQ,order);

    %current global matrix [M + θΔtK]
    GM = GMmass + theta*dt*GMstiffness;
    %previous global matrix [M - (1-θ)ΔtK]
    prevGM = GMmass - ((1-theta)*dt*GMstiffness);

    %initialise list of vector c with zeros 
    c = zeros(mesh.ngn,length(time.t));
    %insert given time value at t = 0
    c(:,1) = time.IC;

    %get global source vector
    [GVsource] = GlobalSourceVector(Ne,mesh,GQ,order);
    
    %Neumann boundary conditions, current and following
    [NB,NBnext] = NeumannBoundary(boundary,mesh);

    %compute RHS therefore c at each time point
    for n = 1:time.N
    
        %RHS vector of equation
        RHS = (prevGM*c(:,n)) + dt*theta*(GVsource+NBnext) + dt*(1-theta)*(GVsource + NB);
        
        %Dirichlet boundary conditionsfor RHS and GM
        [RHS,GM] = DirichletBoundary(boundary,mesh,GM,RHS);
        
        %next numerical according to equation
        c(:,n+1) = GM\RHS;
    
    end
end