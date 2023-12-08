function [L2Norm,h] = L2Norm(Xmin,Xmax,elements,order,theta,time,GQ,boundary,parameters)
%Computes the L2 Norm of the transient diffusion reaction equation.
% The L2 Norm is the error between the numerical and the analytical
% solution. This is done to test convergence of the TransientFEM calculator
% over a range of elements.
%
% Input:
%  Xmin : Lower spatial boundary
%  Xmax : Upper spatial boundary
%  elements : given elements
%  order : weather the basis functions is linear or quadratic
%  theta : Method selection
%  time : All time related values combined
%  GQ : Gaussian Quadrature parameters
%  boundary : Dirichlet and Neumann conditions combined
%  parameters : Parameters for current material
% Return:
%  L2Norm : L2 Norm error
%  h : Characteristic length
%
%Francesco Berteau (fb552) - November 2023

    %number of gauss points
    N = GQ.npts;
    %preallocations for fast processing
    h = zeros(1,length(elements));
    L2Norm = zeros(1,length(elements));
    E = zeros(1,N);

    for b = 1:length(elements)
        %Defines number of elements
        Ne = elements(b);
        Etot = zeros(1,Ne); %preallocation for fast processing
        %length for each element
        h(b) = ((Xmax-Xmin)/Ne);
        %solve transient FEM
        [Cnum,mesh,GQ,time] = TransientFEM(Xmin,Xmax,Ne,order,theta,time,GQ,boundary,parameters);
        %Returns the solution column at the selected time
        cnumerical = Cnum(:,time.t == time.range);
        for eN = 1:Ne
            %Jacobian for a given element
            J = mesh.elem(eN).J;
            %rows of the local nodes at given eN
            Cnoderows = cnumerical(mesh.elem(eN).n(1):mesh.elem(eN).n(end),1);
            for i = 1:N
                w = GQ.gw(i);          %Gauss weights
                xipts = GQ.xipts(i);   %Gauss points
                %basis function at given point
                [psi,~] = EvalBasis(order,xipts);
                %interpolation between nodes
                CNum = psi'*Cnoderows;
                %Finds x position
                x = mesh.elem(eN).x*psi;
                %analytical solution
                CAnalyt = TransientAnalyticSoln(x,time.range);
                %error squared
                E(i) = w*J*(CAnalyt-CNum)^2;
            end
            %total error
            Etot(eN) = sum(E);
        end
        %L2 Norm is RMS of total error
        L2Norm(b) = sqrt(sum(Etot));
    end
    %coefficients of polynomial to fit the data
    y = polyfit(log(h),log(L2Norm),1);
    %error gradient
    gradient = y(1);
end