function [K,teff,position,element] = MinEffectiveDose(mesh,time,order,c,Xpos,ceff)
    
    % find in c(Xpos,:) where the effectiveness concentration is reached
    element = round(order*Xpos/(mesh.nvec(end)/mesh.ne));
    position = find(c(element,:) > ceff,1,'first');
    % use the index to return a time point
    teff = time.t(position);
    
    % compute integral between current time and final of c in dt
    K = trapz(c(element,position:end))*time.dt;


end