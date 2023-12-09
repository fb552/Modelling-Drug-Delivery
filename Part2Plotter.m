%This script solves the drug delivery through uman skin problem.
% Various graphs are used to show the results and to do some analysis for
% Part 2 of ME40064, Transient MATLAB-Based FEM Modelling, CourseWork.
%
%Francesco Berteau (fb552) - November 2023
close all

%% ------- Imput parameters ----------------------------------------------%

% Parameters of current material
parameters.selection = '2';  %material for Part 2
parameters.De = 4e-6;        %Diffusion coefficient epidermis
parameters.Dd = 5e-6;        %Diffusion coefficient dermis
parameters.Db = 2e-6;        %Diffusion coefficient sub-cutaneous
parameters.betaE = 0;        %Reaction coefficient blood epidermis
parameters.betaD = 0.01;     %Reaction coefficient blood dermis
parameters.betaB = 0.01;     %Reaction coefficient blood sub-cutaneous
parameters.gammaE = 0.02;    %Reaction coefficient degradation epidermis
parameters.gammaD = 0.02;    %Reaction coefficient degradation dermis
parameters.gammaB = 0.02;    %Reaction coefficient degradation sub-cutaneous
parameters.Xe = 0.00166667;  %Epidermis x coordinate
parameters.Xd = 0.005;       %Dermis x coordinate
parameters.Xb = 0.01;        %Sub-cutaneous x coordinate

% Space related values
Xmin = 0;   %lower spatial boundary
Xmax = 0.01;   %upper spatial boundary
Ne = 40;    %number of elements
order = 2;  %order (linear(1) or quadratic(2)) of basis function
theta = 0.5;  %0:Forward Euler, 1:Backward Euler, 0.5:Crank-Nicolson

% Time related values
time.tmin = 0;  %start time
time.tmax = 30; %end time
time.ic = 0;    %initial condition time
time.dt = 0.01; %delta t
time.t = time.tmin:time.dt:time.tmax;   %entire time vector
time.N = (time.tmax-time.tmin)/time.dt; %number of time steps
tpoints = [0.05 0.1 0.3 1 2 3 5 7 10 20 30];   %given time points

GQ.switch = '1';    %Gaussian Quadrature yes or no
GQ.npts = 4;        %number of gauss points (2N-1)
[GQ] = GQscheme(GQ);%create Gaussian Quadrature

% All boundaries combined
boundary.DirichletL = 30;   %Lower Dirichlet Boundary Condition
boundary.DirichletU = 0;    %Upper Dirichlet Boundary Condition
boundary.NeumannL = 'Na';   %Lower Neumann Boundary Condition
boundary.NeumannU = 'Na';   %Upper Neumann Boundary Condition

%colours for plotting (parula colormap)
colours = parula(length(tpoints));

%Calculates numerical solution for the transient problem
[c,mesh,GQ,time] = TransientFEM(Xmin,Xmax,Ne,order,theta,time,GQ,boundary,parameters);

%% ------- Drug concentration --------------------------------------------%
figure('Name', 'Concentration')

%Loop through all the selected time points
for i = 1:length(tpoints)
    
    %Plots temperature distribution through the tissue for a range of
    %time points
    element = time.t == tpoints(i);
    plot(mesh.nvec,c(:,element),'DisplayName',(strcat('t=',num2str(tpoints(i)))),'LineWidth',1.3,'color',colours(i,:));
    hold on
end   

%show where the epidermis and dermis end
xline(parameters.Xe,'-.k','Epidermis','FontSize',12,'LineWidth',0.75,'HandleVisibility','off');
xline(parameters.Xd,'-.k','Dermis','FontSize',12,'LineWidth',0.75,'HandleVisibility','off');

grid on %use grid lines
title('Drug concentration through the skin','FontSize',14)
xlabel('Skin Depth (m)','FontSize',12);
ylabel('Concentration  c(x,t)','FontSize',12);
legend('Location','NorthEast','FontSize',10,'NumColumns',2)
ylim([0,35])
% save plot as picture
saveas(gcf,'DrugConcentration','png')


%% ------- Drug effectiveness --------------------------------------------%
figure('Name', 'Effectiveness')

% concentration for effectiveness at given position
ceff = 40;
Xpos = 0.005;

% increase boundary dose untill integral of c is less that 1000
K = 0;
while K < 1000
    [c,mesh,GQ,time] = TransientFEM(Xmin,Xmax,Ne,order,theta,time,GQ,boundary,parameters);
    
    [K,teff,position,element] = MinEffectiveDose(mesh,time,order,c,Xpos,ceff);
    
    % increase dose applied to skin
    boundary.DirichletL = boundary.DirichletL + 1;
end

plot(time.t,c(element,:),'DisplayName','Concentration @ x=0.005','LineWidth',1.3);
hold on
plot(teff,c(element,position),'ro', 'MarkerSize', 10,'DisplayName','Minimum effective dose')
% Add label with coordinates
label = sprintf('[%.2f, %.1f]', teff, c(element,position));
text(teff + 1, c(element,position) - 1, label, 'FontSize', 10);

grid on %use grid lines
title('Drug effectiveness at Dermis end','FontSize',14)
xlabel('Time t (s)','FontSize',12);
ylabel('Concentration c(0.005,t)','FontSize',12);
legend('Location','SouthEast','FontSize',10)
% save plot as picture
saveas(gcf,'DrugEffect','png')

fprintf("Minimum initial dose: %d \n",boundary.DirichletL);

%% ------- Effective drug concentration ----------------------------------%
figure('Name', 'Concentration')

%Loop through all the selected time points
for i = 1:length(tpoints)
    
    %Plots temperature distribution through the tissue for a range of
    %time points
    element = time.t == tpoints(i);
    plot(mesh.nvec,c(:,element),'DisplayName',(strcat('t=',num2str(tpoints(i)))),'LineWidth',1.3,'color',colours(i,:));
    hold on
end   

%show where the epidermis and dermis end
xline(parameters.Xe,'-.k','Epidermis','FontSize',12,'LineWidth',0.75,'HandleVisibility','off');
xline(parameters.Xd,'-.k','Dermis','FontSize',12,'LineWidth',0.75,'HandleVisibility','off');
%show the minimum concentration for effectiveness
yline(ceff,'--r','Effectiveness','FontSize',12,'LineWidth',1,'HandleVisibility','off')

grid on %use grid lines
title('Effective drug concentration through the skin','FontSize',14)
xlabel('Skin Depth (m)','FontSize',12);
ylabel('Concentration  c(x,t)','FontSize',12);
legend('Location','NorthEast','FontSize',10,'NumColumns',2)
ylim([0,75])
% save plot as picture
saveas(gcf,'EffectiveDrugConcentration','png')

%% ------- Parameters influence ------------------------------------------%
figure('Name', 'Parameters')

for i = [1 2 3 4 5 6 7 8]
    parameters = CombinationsDR(i,parameters);
    boundary.DirichletL = 30;
    % increase boundary dose untill integral of c is less that 1000
    K = 0;
    while K < 1000
        [c,mesh,GQ,time] = TransientFEM(Xmin,Xmax,Ne,order,theta,time,GQ,boundary,parameters);
        
        [K,teff,position,element] = MinEffectiveDose(mesh,time,order,c,Xpos,ceff);
        
        % increase dose applied to skin
        boundary.DirichletL = boundary.DirichletL + 1;
    end
    
    fprintf("Minimum initial dose: %d, with combination %d \n",boundary.DirichletL,i);
    plot(time.t,c(element,:),'DisplayName',strcat('Solution for combination: ',num2str(i)),'LineWidth',1.3);
    hold on
end

yline(ceff,'--r','LineWidth',0.75,'DisplayName','Minimum effective dose')
grid on %use grid lines
title('Drug effectiveness with different combinations','FontSize',14)
xlabel('Time t (s)','FontSize',12);
ylabel('Concentration c(0.005,t)','FontSize',12);
legend('Location','SouthEast','FontSize',10)
ylim([0,45])
% save plot as picture
saveas(gcf,'ParameterComparison','png')