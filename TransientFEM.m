function [c,mesh,GQ,time] = TransientFEM(Xmin,Xmax,Ne,order,theta,time,GQ,boundary,parameters)
%Solves the full transient form of the diffusion reaction equation.
% [M + θΔtK]cn+1 = [M  (1 - θ)ΔtK]cn + Δtθ[Fn+1 + NBcn+1] + Δt(1 - θ)[Fn + NBcn]
%
% Input:
%  xmin : Lower spatial boundary
%  xmax : Upper spatial boundary
%  Ne : Number of elements
%  order : weather the basis functions is linear or quadratic
%  theta : Method selection
%  time : All time related values combined
%  GQ : Gaussian Quadrature parameters
%  boundary : Dirichlet and Neumann conditions combined
%  parameters : Parameters for current material
% Return:
%  c : Solution to the full transient FEM
%  mesh : Finite element mesh
%  GQ : Gaussian Quadrature parameters
%  time : All time related values combined
%
%Francesco Berteau (fb552) - November 2023

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
    %initial condition c(x,0)=initial condition
    c(:,1) = time.ic;

    %get global source vector
    [GVsource] = GlobalSourceVector(Ne,mesh,GQ,order);
    
    %Neumann boundary conditions, current and following
    [NB,NBnext] = NeumannBoundary(boundary,mesh);

    %compute RHS hence c at each time point
    for n = 1:time.N
    
        %RHS vector of equation
        RHS = (prevGM*c(:,n)) + dt*theta*(GVsource+NBnext) + dt*(1-theta)*(GVsource + NB);
        
        %Dirichlet boundary conditionsfor RHS and GM
        [RHS,GM] = DirichletBoundary(boundary,RHS,GM);
        
        %next numerical according to equation
        c(:,n+1) = GM\RHS;
    
    end
end
