%This script runs the transient FEM solver and plots the L2 norm used to 
% test the convergence rate of the finite element method function.
%
%Francesco Berteau (fb552) - November 2023
close all

%------- Imput parameters for L2 norm ------------------------------------%

% Parameters of current material
parameters.selection = '1';  %material for Part 1
parameters.D = 1;            %Diffusion coefficient
parameters.lambda = 0;       %Reaction coefficient
parameters.f = 0;            %Source term

Xmin = 0;   %lower spatial boundary
Xmax = 1;   %upper spatial boundary
Ne = 10;    %number of elements
elements = [2 5 10 15 20 25 30 35 40 45]; %choosen elements

% Time related values
time.tmin = 0;  %start time
time.tmax = 1;  %end time
time.ic = 0;    %initial conditiopn time
time.dt = 0.0001; %delta t
time.t = time.tmin:time.dt:time.tmax;   %entire time vector
time.N = (time.tmax-time.tmin)/time.dt; %number of time steps
tpoints = [0.05 0.1 0.3 1.0];   %given time points

GQ.switch = '1';    %Gaussian Quadrature yes or no
GQ.npts = 3;        %number of gauss points (2N-1)
[GQ] = GQscheme(GQ);%create Gaussian Quadrature

% All boundaries combined
boundary.DirichletL = 0;    %Lower Dirichlet Boundary Condition
boundary.DirichletU = 1;    %Upper Dirichlet Boundary Condition
boundary.NeumannL = 'Na';   %Lower Neumann Boundary Condition
boundary.NeumannU = 'Na';   %Upper Neumann Boundary Condition


%------- Linear Vs Quadratic Basis Function ------------------------------%
%Loops through selected elements range switching between linear and quadratic basis functions
for order = [1 2]
    figure('Name',strcat('L2 Norm test, basis order = ',num2str(order)))
    %predefined line colours [rgb]
    linecolor = {[0.00 0.45 0.75],[0.85 0.35 0.10], [0.95 0.70 0.15], [0.50 0.20 0.55]};
    %Crank-Nicolson or Backward Euler methods
    for theta = [0.5 1]

        for i = 1:length(tpoints)
            %get L2 Norm error 
            time.range = tpoints(i);
            [L2N, h] = L2Norm(Xmin,Xmax,elements,order,theta,time,GQ,boundary,parameters);
            %Crank-Nicolson or Backward Euler methods
            switch theta
                case 0.5
                    %plot L2 Norm error against characteristic lenght using log-log plot
                    plot(log(h),log(L2N),'DisplayName',strcat('t=',num2str(time.range)),'LineStyle','-.','LineWidth',1.3,'color',linecolor{i});
                case 1
                    %plot L2 Norm error against characteristic lenght using log-log plot
                    plot(log(h),log(L2N),'DisplayName',strcat('t=',num2str(time.range)),'LineStyle','--','LineWidth',1.3,'color',linecolor{i});
            end
            hold on
        end
    end
    set(gca, 'XScale', 'log', 'YScale', 'log'); %use log axis
    grid on %use grid lines
    title(strcat('L2 Norm with basis order=',num2str(order)),'FontSize',14)
    xlabel('ln(h)','FontSize',12);
    ylabel('ln(E)','FontSize',12);
    legend('Location','SouthEast','FontSize',10,'NumColumns',2)
    % save plot as picture
    saveas(gcf,strcat('L2Norm_order=',num2str(order)),'png')
end


%------- Manual integration Vs. Gaussian Quadrature ----------------------%
order = 1;
time.range = time.tmax;
figure('Name','Numerical Methods Convergence')
%manual or gaussian quadrature integration methods
for j = [0 1]
    GQ.Switch = num2str(j);
    %Crank-Nicolson or Backward Euler methods
    for theta = [0.5 1]
        %get L2 Norm error
        [L2N, h] = L2Norm(Xmin,Xmax,elements,order,theta,time,GQ,boundary,parameters);
        %two line colours
        if theta == 0.5
            linecolor = [0.00 0.45 0.75];
        elseif theta == 1
            linecolor = [0.85 0.35 0.10];
        end
        %manual or gaussian quadrature
        switch j
            case 0
                %plot L2 Norm error against characteristic lenght using ln plot
                plot(log(h),log(L2N),'LineWidth',1.3,'color',linecolor);
            case 1
                %plot L2 Norm error against characteristic lenght using ln plot
                plot(log(h),log(L2N),'--','LineWidth',1.3,'color',linecolor);
        end
        hold on
    end
end
set(gca, 'XScale', 'log', 'YScale', 'log'); %use log axis
grid on %use grid lines
title('Numerical Methods Convergence','FontSize',14)
xlabel('ln(h)','FontSize',12);
ylabel('ln(E)','FontSize',12)
legend('Crank-Nicolson (Manual Integration)','Backward Euler (Manual Integration)',...
'Crank-Nicolson (Gaussian Quadrature)','Backward Euler (Gaussian Quadrature)','Location','SouthEast','FontSize',10)
% save plot as picture
saveas(gcf,'NumericalMethodsConvergence','png')