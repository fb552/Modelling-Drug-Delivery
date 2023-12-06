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
order = 1;  %order (linear(1) or quadratic(2)) of basis function
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
GQ.npts = 3;        %number of gauss points (2N-1)
[GQ] = GQscheme(GQ);%create Gaussian Quadrature

% All boundaries combined
boundary.DirichletL = 30;   %Lower Dirichlet Boundary Condition
boundary.DirichletU = 0;    %Upper Dirichlet Boundary Condition
boundary.NeumannL = 'Na';   %Lower Neumann Boundary Condition
boundary.NeumannU = 'Na';   %Upper Neumann Boundary Condition

%colours for plotting (parula colormap)
colours = parula(length(tpoints));

%Calculates numerical solution for the transient problem
[c,mesh,GQ,time,GM] = TransientFEM(Xmin,Xmax,Ne,order,theta,time,GQ,boundary,parameters);

%% ------- Drug concentration --------------------------------------------%
figure('Name', 'Concentration')

%Loop through all the selected time points
for i = 1:length(tpoints)
    
    %Plots temperature distribution through the tissue for a range of
    %time points
    findplace = time.t == tpoints(i);
    plot(mesh.nvec,c(:,findplace),'DisplayName',(strcat('t=',num2str(tpoints(i)))),'LineWidth',1.3,'color',colours(i,:));
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
[c,mesh,GQ,time,GM] = TransientFEM(Xmin,Xmax,Ne,order,theta,time,GQ,boundary,parameters);

%time step
dt = time.dt;

element = round(order*0.005/(mesh.nvec(end)/mesh.ne));

position = find(c(element,:) > 40,1,'first');
teff = time.t(position);
K = trapz(c(element,position:end))*dt;

plot(time.t,c(element,:));
hold on
plot([1 1]*time.t(position),[0 c(element,position)], '--r')



fprintf(strcat(" Effectiveness: ",num2str(DrugEffect)),'\n');
% 













% %----------------------Minimum temperature reduction ---------------------%
% 
% G = [0 0.0375]; %Defines blood flow vector
% 
% %Implements linear and quadratic basis functions
% for orderbase = [1 2]
%     
%     %Implements Crank-Nicolson and Backward Euler numerical method
%     for theta = [0.5 1]
%         
%         %Implements blood flow values
%         for j = 1:length(G)
%             
%             matparameters.G = G(j); %Defines blood flow value
%             
%             %Calculates numerical solution for the transient problem
%             [mesh,Cnum,GQ,time,GlobalMatrix] = FiniteElemSolver(Boundary,xmin,xmax,Ne,GQ,orderbase,theta,time,matparameters);
%             
%             %Implements dermis and epidermis boundaries
%             for xtissue = [0.005 0.0016667]
%                 
%                 %Calculates temperature reduction
%                 [TempReduction] = TempReductionCalc(Boundary,xmin,xmax,Ne,GQ,orderbase,theta,time,matparameters,xtissue,Tburn);
%                 %Write formatted blood flow and temperature reduction data to text file
%                 fprintf(strcat('\n'," G = ", num2str(G(j))," Temperature reduction required: ",num2str(TempReduction),"(K)"))
%             end
%         end
%     end
% end
% 
% %-------------------------Convergence Analysis----------------------------%
% 
% G = [0 0.0375]; %Defines blood flow vector
% 
% %Loop to calculate temperature values for each blood flow case
% for j = 1:length(G)
%     
%     matparameters.G = G(j); %Defines blood flow value
%     
%     %Calculates numerical solution for the transient problem
%     [mesh,Cnum,GQ,time,GlobalMatrix] = FiniteElemSolver(Boundary,xmin,xmax,Ne,GQ,orderbase,theta,time,matparameters);
%     
%     figure(5)
%     
%     %Plots temperature against time at the dermis layer
%     findplace = mesh.nvec == 0.005;
%     plot(time.t,Cnum(findplace,:),'DisplayName',(strcat('G = ',num2str(G(j)))),...
%     'LineWidth',1.3,'color',linecolor{j});
%     hold on
%     
%     grid on %Sets grid lines
%     xlabel('Time t (s)','FontSize',12); %Creates xlabel
%     ylabel('Temperature T(0.005,t)','FontSize',12); %Creates ylabel
%     legend('Location','SouthEast','FontSize',10,'NumColumns',1) %Creates legend
% end