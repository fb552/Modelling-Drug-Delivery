%This script runs the transient FEM solver and plots some graphs used to 
% prove the code functioning for Part 1 of ME40064, Transient MATLAB-Based 
% FEM Modelling, CourseWork.
close all

%------- Input parameters ------------------------------------------------%

parameters.selection = '1';  %material for Part 1
parameters.D = 1;            %Diffusion coefficient
parameters.lambda = 0;       %Reaction coefficient
parameters.f = 0;            %Source term

Xmin = 0;   %lower spatial boundary
Xmax = 1;   %upper spatial boundary
Ne = 10;    %number of elements
order = 1;  %Selects linear basis functions
theta = 1;  %0:Forward Euler, 1:Backward Euler, 0.5:Crank-Nicolson
Xpoints = Xmin:0.01:Xmax; %vector with all the needed x points

time.tmin = 0;  %start time
time.tmax = 1;  %end time
time.IC = 0;    %initial time
time.dt = 0.01; %delta t
time.t = time.tmin:time.dt:time.tmax;   %entire time vector
time.N = (time.tmax-time.tmin)/time.dt; %number of time steps
tpoints = [0.05 0.1 0.3 1.0];   %given time points

GQ.Switch = '1';    %Gaussian Quadrature yes or no
GQ.N = 2;           %order of solver
[GQ] = GQscheme(GQ);%create Gaussian Quadrature

boundary.DirichletL = 0;    %Lower Dirichlet Boundary Condition
boundary.DirichletR = 1;    %Upper Dirichlet Boundary Condition
boundary.NeumannL = 'Na';   %Lower Neumann Boundary Condition
boundary.NeumannR = 'Na';   %Upper Neumann Boundary Condition

%get numerical solution for the transient FEM
[mesh,Cnum,GQ,time,GlobalMatrix] = TransientFEM(boundary,Xmin,Xmax,Ne,GQ,order,theta,time,parameters);


%------- Graphs for 1.1 --------------------------------------------------%
%first figure for computed solution alone
figure('Name','Part1numerical')
for i = 1:length(tpoints)

    %Plot c(x) against x
    findplace = time.t == tpoints(i);
    plot(mesh.nvec,Cnum(:,findplace),'DisplayName',(strcat('t =',num2str(tpoints(i)))),'LineWidth',1.3);
    hold on

    grid on %use grid lines
    xlabel('Spatial distance x','FontSize',12)
    ylabel('Concentration c(x,t)','FontSize',12)
    legend('Location','NorthWest','FontSize',10)

end
%second figure for analytical vs numerical comparison
figure('Name', 'Part1analytical')
%Defines line colours using MATLAB RGB triplet
linecolor = {[0.00 0.45 0.74],[0.85 0.33 0.10], [0.93 0.69 0.13], [0.49 0.18 0.56]};
Canalytical = zeros(length(Xpoints),length(tpoints));

for i = 1:length(tpoints)
    for m = 1:length(Xpoints)
        %Calculates analytical solution
        Canalytical(m,i) = TransientAnalyticSoln(Xpoints(m),tpoints(i));
    end

    %plot analytical solution against x
    plot(Xpoints, Canalytical(:,i),'DisplayName',strcat('Analytical Solution at t = ',num2str(tpoints(i))),'LineStyle','--','LineWidth',1.3,'color', linecolor{i});
    hold on

    %plot c(x) against x
    findplace = time.t == tpoints(i);
    plot(mesh.nvec,Cnum(:,findplace),'DisplayName',(strcat('Numerical Solution at t=',num2str(tpoints(i)))),'LineWidth',1.3,'color',linecolor{i});
    hold on

    grid on %use grid lines
    xlabel('Spatial distance x','FontSize',12)
    ylabel('Concentration c(x,t)','FontSize',12)
    legend('Location','NorthWest','FontSize',10)
end

%------- Graphs for 1.2 --------------------------------------------------%
%figure for analytical vs numerical comparison @ x = 0.8
figure('Name','Part1analytical@0.8')

findplace = mesh.nvec == 0.8;
plot(time.t,Cnum(findplace,:),'DisplayName','Numerical Solution x=0.8','LineWidth',1.3);
hold on

Canalytical = TransientAnalyticSoln(0.8,time.t);
plot(time.t,Canalytical,'DisplayName','Analytical Solution x=0.8','LineWidth',1.3);
hold on

grid on %use grid lines
xlabel('Time t (s)','FontSize',12)
ylabel('Concentration c(0.8,t)','FontSize',12)
legend('Location','SouthEast','FontSize',10)