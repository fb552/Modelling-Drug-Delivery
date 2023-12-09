function [K,teff,position,element] = MinEffectiveDose(mesh,time,order,c,Xpos,ceff)
%Computes the Minimum effective drug dose
% Given a position in the mesh and a value for the concentration for
% effectivensess this returns the integral of the curve between the
% effectiveness point and the end.
%
% Input:
%  mesh : Finite element mesh
%  time : All time related values combined
%  order : weather the basis functions is linear or quadratic
%  c : Solution to the full transient FEM
%  Xpos : Position along the mesh
%  ceff : effectiveness concentration
% Return:
%  K : current Neumann boundary conditions vector
%  teff : time of effectiveness
%  position : index of given Xpos
%  element : index of effective time
%
%Francesco Berteau (fb552) - November 2023
    
    % find in c(Xpos,:) where the effectiveness concentration is reached
    element = round(order*Xpos/(mesh.nvec(end)/mesh.ne));
    position = find(c(element,:) > ceff,1,'first');
    % use the index to return a time point
    teff = time.t(position);
    
    % compute integral between current time and final of c in dt
    K = trapz(c(element,position:end))*time.dt;
end