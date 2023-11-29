function [gq] = GQscheme(gq)
%Creates a Gauss-Legendre Quadrature scheme data structure of order N
% The scheme stores both the quadrature weights and the Legendre points.
% Values found @
% https://en.wikipedia.org/wiki/Gauss-Legendre_quadrature#Definition
%
% Input:
%  GQ : Gaussian Quadrature initial parameters
% Return:
%  GQ : Gaussian Quadrature with weights and points
%
%Francesco Berteau (fb552) - November 2023

    %order of quadrature rule
    N = gq.N;
    
    if (N>0) && (N<6)
        %initialise zero array of size equal gauss quadrature order
        gq.gw = zeros(N,1);     %Gauss weights
        gq.xipts = zeros(N,1);  %Gauss points
    
        switch N
            case 1  
              gq.gw(1) = 2; 
              gq.xipts(1) = 0;
            case 2
              gq.gw(1) = 1;  
              gq.gw(2) = 1;  
              gq.xipts(1) = -sqrt(1/3);
              gq.xipts(2) = sqrt(1/3);
            case 3
              gq.gw(1) = 5/9;  
              gq.gw(2) = 8/9; 
              gq.gw(3) = 5/9; 
              gq.xipts(1) = -sqrt(3/5);
              gq.xipts(2) = 0;
              gq.xipts(3) = sqrt(3/5);
            case 4
              gq.gw(1) = (18+sqrt(30))/36;  
              gq.gw(2) = (18+sqrt(30))/36;                
              gq.gw(3) = (18-sqrt(30))/36;
              gq.gw(4) = (18-sqrt(30))/36;  
              gq.xipts(1) = -sqrt((3/7)-(2/7)*sqrt(6/5));
              gq.xipts(2) = sqrt((3/7)-(2/7)*sqrt(6/5));
              gq.xipts(3) = -sqrt((3/7)+(2/7)*sqrt(6/5));  
              gq.xipts(4) = sqrt((3/7)+(2/7)*sqrt(6/5));
            case 5
              gq.gw(1) = 128/225; 
              gq.gw(2) = (322+13*sqrt(70))/900;
              gq.gw(3) = (322+13*sqrt(70))/900;              
              gq.gw(4) = (322-13*sqrt(70))/900;  
              gq.gw(5) = (322-13*sqrt(70))/900;  
              gq.xipts(1) = 0.0;
              gq.xipts(2) = -(1/3)*sqrt(5-2*sqrt(10/7));
              gq.xipts(3) = (1/3)*sqrt(5-2*sqrt(10/7));
              gq.xipts(4) = -(1/3)*sqrt(5+2*sqrt(10/7));
              gq.xipts(5) = (1/3)*sqrt(5+2*sqrt(10/7));
        end
    else
      fprintf('Invalid number of Gauss points (N)');
    end 
end