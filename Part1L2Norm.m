%This script is used to plot all the required graphs to complete the L2
%Norm Extension of the Part 1 of the Transient MATLAB-Based FEM Modelling.
%This extension analyses the convergence rate of the finite element method.
close all


%--------------------Input parameters for L2 Norm-------------------------%
xmin = 0; %Lower spatial boundary position
xmax = 1; %Upper spatial boundary position
Ne = 10; %Number of elements
GQ.Switch = '1'; %Use of Gaussian Quadrature
GQ.N = 3; %Order of solver value
[GQ] = GQscheme(GQ); %Generates Gaussian quadrature
time.tmin = 0; %Start value of time interval
time.tmax = 1; %Last value of time interval
time.IC = 0; %Initial time condition
time.dt = 0.0001; %Time step
time.t = time.tmin:time.dt:time.tmax; %Creates the time interval vector
time.N = (time.tmax-time.tmin)/time.dt; %Number of time steps
Boundary.DirichletL = 0; %Defines Lower Dirichlet boundary
Boundary.DirichletU = 1; %Defines Upper Dirichlet boundary
Boundary.NeumannL = 'Na'; %Defines Lower Neumann boundary
Boundary.NeumannU = 'Na'; %Defines Upper Neumann boundary
matparameters.selection = '1'; %Selects material parameters
matparameters.D = 1; %Diffusion coefficient
matparameters.lambda = 0; %Reaction coefficient
matparameters.f = 0; %Source term
elementsrange = [2 5 10 15 20 25 30 35 40 50 60]; %Selected element size
timepoints = [0.05 0.1 0.3 1.0]; %Selected time points


%------------------Linear Vs Quadratic Basis Functions -------------------%
%Loops through selected elements range switching between linear and
%quadratic basis functions
for orderbase = [1 2]
    figure('Name',strcat('L2Norm Testing with basis order = ',num2str(orderbase)))
    %%Defines line colours using MATLAB RGB triplet
    linecolor = {[0.00 0.45 0.74],[0.85 0.33 0.10], [0.93 0.69 0.13], [0.49...
    0.18 0.56]};
    %Implements Crank-Nicolson and Backward Euler Numerical method
    for theta = [0.5 1]
        %Loops through all the selected time points
        for i = 1:length(timepoints)
            %Calculates L2Norm error for the selected time points
            time.range = timepoints(i);
            [L2Norm, h, gradient] = L2Normapproximation(xmin,xmax,elementsrange,Boundary,GQ,...
            orderbase,theta,time,matparameters);
            %Switch between Crank-Nicolson and Backward Euler method
            switch theta
                case 0.5
                    %Plots L2Norm error against characteristic lenght using a Log-log
                    %scale plot
                    loglog(h,L2Norm,'DisplayName',strcat('t = ',num2str(time.range)),...
                    'LineStyle','--','LineWidth',1.3,'color',linecolor{i});
                case 1
                    %Plots L2Norm error against characteristic lenght using a Log-log
                    %scale plot
                    loglog(h,L2Norm,'DisplayName',strcat('t =',num2str(time.range)),'LineWidth',1.3,'color',linecolor{i});
            end
            hold on
            grid on %Sets grid lines
            xlabel('log(h)','FontSize',12); %Creates xlabel
            ylabel('log(L2Norm)','FontSize',12); %Creates ylabel
            legend('Location','SouthEast','FontSize',10,'NumColumns',2) %Creates
            legend
        end
    end
end


%-----------------Manual integration Vs. Gaussian quadrature -------------------%
orderbase = 1;
time.range = time.tmax;
figure('Name','Numerical Methods Convergence ')
%Loops through selected elements range switching manual integration and
%gaussian quadrature
for j = [0 1]
    GQ.Switch = num2str(j);
    %Implements Crank-Nicolson and Backward Euler Numerical method
    for theta = [0.5 1]
        %Calculates L2Norm error
        [L2Norm, h, gradient] = L2Normapproximation(xmin,xmax,elementsrange,Boundary,GQ,...
        orderbase,theta,time,matparameters);
        %Defines line colours according to the numerical method
        %implemented
        if theta == 0.5
            linecolor = [0.00 0.45 0.74];
        elseif theta == 1
            linecolor = [0.85 0.33 0.10];
        end
        %Switch between manual and gaussian quadrature integration method
        switch j
            case 0
                %Plots L2Norm error against characteristic lenght using a Log-log
                %scale plot
                loglog(h,L2Norm,'LineWidth',1.3,'color',linecolor);
            case 1
                %Plots L2Norm error against characteristic lenght using a Log-log
                %scale plot
                loglog(h,L2Norm,'-x','LineWidth',1.3,'color',linecolor);
        end
        hold on
        grid on %Sets grid lines
        xlabel('log(h)','FontSize',12); %Creates xlabel
        ylabel('log(L2Norm)','FontSize',12) %Creates ylabel
        %Creates legend
        legend('Crank-Nicolson (Manual)','Backward Euler (Manual)',...
        'Crank-Nicolson (GQ Scheme)','Backward Euler (GQScheme)','Location','SouthEast','FontSize',10)
    end
end