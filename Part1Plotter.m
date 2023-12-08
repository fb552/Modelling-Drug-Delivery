%This script runs the transient FEM solver and plots some graphs used to 
% prove the code functioning for Part 1 of ME40064, Transient MATLAB-Based 
% FEM Modelling, CourseWork.
%
%Francesco Berteau (fb552) - November 2023
close all

%% ------- Imput parameters ----------------------------------------------%

% Parameters of current material
parameters.selection = '1';  %material for Part 1
parameters.D = 1;            %Diffusion coefficient
parameters.lambda = 0;       %Reaction coefficient
parameters.f = 0;            %Source term

Xmin = 0;   %lower spatial boundary
Xmax = 1;   %upper spatial boundary
Ne = 10;    %number of elements
order = 1;  %order (linear(1) or quadratic(2)) of basis function
theta = 1;  %0:Forward Euler, 1:Backward Euler, 0.5:Crank-Nicolson
Xpoints = Xmin:0.01:Xmax; %vector with x points for analytical plotting

% Time related values
time.tmin = 0;  %start time
time.tmax = 1;  %end time
time.ic = 0;    %initial condition time
time.dt = 0.01; %delta t
time.t = time.tmin:time.dt:time.tmax;   %entire time vector
time.N = (time.tmax-time.tmin)/time.dt; %number of time steps
tpoints = [0.05 0.1 0.3 1.0];   %given time points

GQ.switch = '1';    %Gaussian Quadrature yes or no
GQ.npts = 2;        %number of gauss points (2N-1)
[GQ] = GQscheme(GQ);%create Gaussian Quadrature

% All boundaries combined
boundary.DirichletL = 0;    %Lower Dirichlet Boundary Condition
boundary.DirichletU = 1;    %Upper Dirichlet Boundary Condition
boundary.NeumannL = 'Na';   %Lower Neumann Boundary Condition
boundary.NeumannU = 'Na';   %Upper Neumann Boundary Condition

%get numerical solution for the transient FEM
[Cnum,mesh,GQ,time,GM] = TransientFEM(Xmin,Xmax,Ne,order,theta,time,GQ,boundary,parameters);


%% ------- Graphs for 1.1 ------------------------------------------------%
%first figure for computed solution alone
figure('Name','Part1 numerical')
for i = 1:length(tpoints)

    %Plot c(x) against x
    element = time.t == tpoints(i);
    plot(mesh.nvec,Cnum(:,element),'DisplayName',(strcat('t=',num2str(tpoints(i)))),'LineWidth',1.3);
    hold on
end
grid on %use grid lines
title('Transient FEM numerical solutions','FontSize',14)
xlabel('Spatial distance x','FontSize',12)
ylabel('Concentration c(x,t)','FontSize',12)
legend('Location','NorthWest','FontSize',10)
% save plot as picture
saveas(gcf,'TransientFEM','png')

%second figure for analytical vs numerical comparison
figure('Name', 'Part1 analytical')
%Defines line colours using MATLAB RGB triplet
colours = {[0.00 0.45 0.74],[0.85 0.33 0.10], [0.93 0.69 0.13], [0.49 0.18 0.56]};
Canalytical = zeros(length(Xpoints),length(tpoints));

for i = 1:length(tpoints)
    for m = 1:length(Xpoints)
        %Calculates analytical solution
        Canalytical(m,i) = TransientAnalyticSoln(Xpoints(m),tpoints(i));
    end

    %plot analytical solution against x
    plot(Xpoints, Canalytical(:,i),'DisplayName',strcat('Analytical Solution @ t=',num2str(tpoints(i))),'LineStyle','--','LineWidth',1.3,'color', colours{i});
    hold on

    %plot c(x) against x
    element = time.t == tpoints(i);
    plot(mesh.nvec,Cnum(:,element),'DisplayName',(strcat('Numerical Solution @ t=',num2str(tpoints(i)))),'LineWidth',1.3,'color',colours{i});
    hold on
end
grid on %use grid lines
title('Numerical vs Analytical solutions','FontSize',14)
xlabel('Spatial distance x','FontSize',12)
ylabel('Concentration c(x,t)','FontSize',12)
legend('Location','NorthWest','FontSize',10)
% save plot as picture
saveas(gcf,'TransientFEM-analytical','png')

%% ------- Graphs for 1.2 ------------------------------------------------%
%figure for analytical vs numerical comparison @ x = 0.8
figure('Name','Part1 analytical @ 0.8')

element = mesh.nvec == 0.8;
plot(time.t,Cnum(element,:),'DisplayName','Numerical Solution @ x=0.8','LineWidth',1.3);
hold on

Canalytical = TransientAnalyticSoln(0.8,time.t);
plot(time.t,Canalytical,'DisplayName','Analytical Solution @ x=0.8','LineWidth',1.3);
hold on

grid on %use grid lines
title('Numerical vs Analytical solution at x=0.8','FontSize',14)
xlabel('Time t (s)','FontSize',12)
ylabel('Concentration c(0.8,t)','FontSize',12)
legend('Location','SouthEast','FontSize',10)
% save plot as picture
saveas(gcf,'TransientFEM-analytical@08','png')

%% ------- Backward Euler vs Crank-Nicolson ------------------------------%
%figure for different theta scheme comparison @ x = 0.8
figure('Name','Crank-Nicolson, Backward and Forward @ x = 0.8')

for theta = [0 0.5 1]
    %get numerical solution for given theta
    [Cnum,mesh,GQ,time,GM] = TransientFEM(Xmin,Xmax,Ne,order,theta,time,GQ,boundary,parameters);
    
    element = mesh.nvec == 0.8;
    plot(time.t,Cnum(element,:),'-x','DisplayName',strcat('Numerical Solution θ=',num2str(theta)),'LineWidth',1.3);
    hold on
end

%Plots numerical solution at x = 0.8 for the established time interval
Canalytical = TransientAnalyticSoln(0.8,time.t);
plot(time.t,Canalytical,'DisplayName','Analytical Solution','LineWidth',1.3);

grid on %use grid lines
title('Crank-Nicolson vs Euler methods','FontSize',14)
xlabel('Time t (s)','FontSize',12)
ylabel('Concentration c(0.8,t)','FontSize',12)
legend('Location','SouthEast','FontSize',10)
ylim([0,0.8])
xlim([0,0.2])
% save plot as picture
saveas(gcf,'Crank-Nicolson_vs_Euler','png')

%% ------- Linear vs Quadratic basis function ----------------------------%
%figure for different basis function order
figure('Name','Basis function order')

Canalytical = zeros(length(Xpoints),1);
for m = 1:length(Xpoints)
    %Calculates analytical solution
    Canalytical(m,1) = TransientAnalyticSoln(Xpoints(m),0.05);
end
%plot analytical solution against x
plot(Xpoints, Canalytical(:,1),'DisplayName','Analytical Solution','LineWidth',1.3);
hold on

for order = [1 2]
    [Cnum,mesh,GQ,time,GM] = TransientFEM(Xmin,Xmax,Ne,order,theta,time,GQ,boundary,parameters);
    %plot c(x) against x
    element = time.t == 0.05;
    plot(mesh.nvec,Cnum(:,element),'-x','DisplayName',(strcat('Numerical Solution order=',num2str(order))),'LineWidth',1.3);
    hold on
end

grid on %use grid lines
title('Basis Function order','FontSize',14)
xlabel('Spatial distance x','FontSize',12)
ylabel('Concentration c(x,0.05)','FontSize',12)
legend('Location','NorthWest','FontSize',10)
xlim([0.4,0.8])
% save plot as picture
saveas(gcf,'BasisFunctionOrder','png')