function [DrugEffect] = MinEffectiveDose(mesh,time,order,c,xtissue,ceff)

    %position of the corresponding layer in the finite element mesh
    findposition = find((mesh.nvec<xtissue),1,'last');

            
    %Establish local nodes
    x0 = mesh.nvec(findposition);
    x1 = mesh.nvec(findposition + 1);
    
    %Calculates xpoints by interpolating the skin tissue layer
    %within the established x domain
    xpoint = (2*((xtissue - x0)/(x1-x0))-1);
    
    %Returns the corresponding value of the basis function at the given
    %Gauss point
    [psi, ~] = EvalBasis(order,xpoint);
    
    %Extracts temperature values for the input tissue layer position
    Tempvector = psi' * c(findposition:findposition + 1,:);
 
    
    %Extracts all the temperature values above the established temperature limit
    %(Tburn)
    effect = Tempvector(Tempvector >= ceff);
    
    %time step
    dt = time.dt;

    %Evaluate the integral numerically using the trapezium rule
    DrugEffect = trapz(effect)*dt;
end