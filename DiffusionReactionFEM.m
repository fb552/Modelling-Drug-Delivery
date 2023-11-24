function [c] = DiffusionReactionFEM(D, lambda, Ne)
% Function with Finite Element Method code to solve diffusion-reaction equation
    Xmin = 0;   %domain lower range
    Xmax = 1;   %domain upper range
    f = 0;      %source vector
    mesh = OneDimLinearMeshGen(Xmin,Xmax,Ne);
    
    %solve the matrix for
    matrixTerm = GlobalMatrix(D, lambda, mesh); %diffusion & reaction
    sourceTerm = GlobalVector(f, mesh);         %source term

    %Boundary Conditions for
    LowerDBC = 0; %lower Dirichlet
    UpperDBC = 1; %upper Dirichlet

    %upper Dirichlet BC
    matrixTerm(end, :) = 0;
    matrixTerm(end) = 1;
    sourceTerm(end) = UpperDBC;

    %lower Dirichlet BC
    matrixTerm(1, :) = 0;
    matrixTerm(1, 1) = 1;
    sourceTerm(1, 1) = LowerDBC;

    %laplace solution
    c = matrixTerm \ sourceTerm;

    %analytical solution
    x = 0:1/Ne:1;           % (0, 0.25, 0.5, 0.75, 1)
    analyticalC = transpose(((exp(3))/((exp(6))-(1)))*((exp((3)*(x)))-(exp((-3)*(x)))));

    %plot
    figure(1) = plot(x, c,'Color', [0 0 1],'Marker','o', 'DisplayName', 'Diffusion-Reaction FEM Solution');
    hold on
    figure(2) = plot(x, analyticalC,'Color', [1 0 0], 'DisplayName', 'Analytical Solution');
    hold off
    xlabel('X')
    ylabel('c')
    legend(figure, 'Location', 'northwest');
    % save plot as picture
    saveas(gcf,'d-rFEM16','png')
end