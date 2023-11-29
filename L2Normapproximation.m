function [L2Norm, h,gradient] = L2Normapproximation(xmin,xmax,elementsrange,Boundary,GQ,orderbase,theta,time,matparameters)
%This function calculates the L2 Norm error between the analytical and
%numerical solution of the transient form of the diffusion-reaction
%equation. This is carried out for a range of elements in order to test
%time and spatial convergence.
%
%Input arguments:
%xmin - Lower spatial boundary position
%xmax - Upper spatial boundary position
%elementsrange - Range of elements
%Boundary - Data structure containing all the boundary conditions
%GQ - Switch between manual integration and Gaussian quadrature rule
%orderbase - Switch between linear and quadratic basis functions
%theta - Selection of numerical scheme
%time - Data structure containing all the time data
%matparameters - Selection of material parameters
%
%Return arguments:
%L2Norm - L2 Norm error
%h - Characteristic length
%gradient - Output line gradient
    %Defines the order of the quadrature rule
    N = GQ.N;
    %Loops through all the elements defined in the elementsrange vector
    for b = 1:length(elementsrange)
        %Defines number of elements
        Ne = elementsrange(b);
        %Calculates the characteristic length for each element
        h(b) = ((xmax-xmin)/Ne);
        %Solves the transient form of the diffusion-reaction equation
        [Cnum,mesh,GQ,time,~] = TransientFEM(xmin,xmax,Ne,orderbase,theta,time,GQ,Boundary,matparameters);
        %Returns the solution column at the selected time
        cnumerical = Cnum(:,time.t == time.range);
        for eID = 1:Ne
            %Access Jacobian elements
            J = mesh.elem(eID).J;
            %Returns the rows of the local nodes at the given eID
            Cnoderows = cnumerical(mesh.elem(eID).n(1):mesh.elem(eID).n(end),1);
            %Loops over number of Gauss points (xpoint)
            for i = 1:N
                weight = GQ.gw(i); %Access Gauss weights
                xpoint = GQ.xipts(i); %Access Gauss points
                %Returns the corresponding value of the basis function at
                %the given Gauss point
                [psi,~] = EvalBasis(orderbase,xpoint);
                %Interpolation to solve between nodes multiplying by the
                %transpose of psi
                CNum = psi'*Cnoderows;
                %Finds x position
                x = mesh.elem(eID).x*psi;
                %Determines analytical solution
                CAnalyt = TransientAnalyticSoln(x,time.range);
                %Calculates the squared error
                E(i) = weight*J*(CAnalyt-CNum)^2;
            end
            %Calculates total error by adding together the error at each
            %element
            Etot(eID) = sum(E);
        end
        %Takes the square root of the error to apply L2Norm concept (RMS)
        L2Norm(b) = sqrt(sum(Etot));
    end
    %Finds the coefficients of a polynomial of degree N that fits the data best
    %in a least-squares sense
    y = polyfit(log(h),log(L2Norm),1);
    %Defines the error gradient
    gradient = y(1);
end