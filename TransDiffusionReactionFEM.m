function [c] = TransDiffusionReactionFEM(D, lambda, Ne, T, dt)
% Function with Finite Element Method code to solve transient diffusion-reaction equation
    Xmin = 0;   % domain lower range
    Xmax = 1;   % domain upper range
    f = 0;      % source vector
    mesh = OneDimLinearMeshGen(Xmin,Xmax,Ne);
    
    % solve the matrix for diffusion & reaction
    matrixTerm = GlobalMatrix(D, lambda, mesh);
    
    % Boundary Conditions
    LowerDBC = 0; % lower Dirichlet
    UpperDBC = 1; % upper Dirichlet

    % Apply Dirichlet BC
    matrixTerm(1, :) = 0;
    matrixTerm(1, 1) = 1;
    
    matrixTerm(end, :) = 0;
    matrixTerm(end, end) = 1;

    % Initial condition
    c0 = zeros(Ne + 1, 1);
    
    % Time-stepping loop (Forward Euler)
    numSteps = round(T / dt);
    for n = 1:numSteps
        % Compute the source term at the current time step
        sourceTerm = GlobalVector(f, mesh) + lambda * c0;
        
        % Apply Dirichlet BC to the source term
        sourceTerm(1) = LowerDBC;
        sourceTerm(end) = UpperDBC;

        % Update the solution using Forward Euler method
        c = c0 + dt * (matrixTerm \ sourceTerm);

        % Set the current solution as the initial condition for the next time step
        c0 = c;
    end

    % Analytical solution
    x = 0:1/Ne:1;
    analyticalC = transpose(((exp(3))/((exp(6))-(1)))*((exp((3)*(x)))-(exp((-3)*(x)))));

    % Plot
    figure;
    plot(x, c, 'Color', [0 0 1], 'Marker', 'o', 'DisplayName', 'Diffusion-Reaction FEM Solution');
    hold on;
    plot(x, analyticalC, 'Color', [1 0 0], 'DisplayName', 'Analytical Solution');
    hold off;
    xlabel('X');
    ylabel('c');
    legend('Location', 'northwest');
    % Save plot as picture
    saveas(gcf, 'd-rFEM16', 'png');
end
