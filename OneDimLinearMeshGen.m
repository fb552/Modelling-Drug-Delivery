function [mesh] = OneDimLinearMeshGen(xmin,xmax,Ne,order,parameters)
%%This function generates a one dimensional, equispaced, linear finite
%%element mesh, with Ne number of elements, between the points at x
%%position xmin and xmax.

    mesh.ne = Ne; %set number of elements
    mesh.ngn = Ne+1 + (Ne*(order -1)); %set number of global nodes
    mesh.nvec = zeros(mesh.ngn,1); %allocate vector to store global node values
    dx = (xmax - xmin)/Ne; %calculate element size

    mesh.nvec = xmin:dx-(0.5*dx*(order-1)):xmax;
    
    %loop over elements & set the element properties
    for i=1:Ne
        %set element Jacobian based on mapping to standard element
        mesh.elem(i).J = 0.5*dx; %this is assuming standard element of -1 to 1

        switch order
            case 1 %linear
                %set spatial positions of nodes
                mesh.elem(i).x(1) = xmin + (i-1)*dx;
                mesh.elem(i).x(2) = xmin + i*dx ;
                %set global IDs of the nodes
                mesh.elem(i).n(1) = i;
                mesh.elem(i).n(2) = i+1;

            case 2 %quadratic
                %set spatial positions of nodes
                mesh.elem(i).x(1) = xmin + (i-1)*dx;
                mesh.elem(i).x(2) = xmin + i*dx - dx/2;
                mesh.elem(i).x(3) = xmin + i*dx ;
                %set global IDs of the nodes
                mesh.elem(i).n(1) = (i*2)-1;
                mesh.elem(i).n(2) = (i*2);
                mesh.elem(i).n(3) = (i*2)+1;
        end

        %check weather using material 1 or 2
        if parameters.selection == '1'

            %diffusion coefficient for each element
            mesh.elem(i).D = parameters.D;

            %reaction coefficient for each element
            mesh.elem(i).lambda = parameters.lambda;

            %source term for each element
            mesh.elem(i).f = parameters.f;

        elseif parameters.selection == '2'
            %human sking parameters
            mesh.elem(i).rho = parameters.rho;  %tissue density
            mesh.elem(i).c = parameters.c;      %tissue specific heat capacity
            mesh.elem(i).G = parameters.G;      %blood flow
            mesh.elem(i).Tb = parameters.Tb;    %blood temperature
            mesh.elem(i).rhob = parameters.rhob;%blood density
            mesh.elem(i).cb = parameters.cb;    %blood Specific heat capacity

            Xparameters = mesh.elem(i).x(1);    %x-position of element

            %within epidermis layer
            if(Xparameters >= 0 && Xparameters < parameters.xEpidermis)
                %set thermal conductivity and assume no blood flow
                mesh.elem(i).k = parameters.kEpidermis;
                mesh.elem(i).G = 0;

            %within dermis layer
            elseif (Xparameters >= parameters.xEpidermis && Xparameters <= parameters.xDermis)
                %set thermal conductivity
                mesh.elem(i).k = parameters.kDermis; 

            %within subcutaneous layer
            elseif (Xparameters > parameters.xDermis && Xparameters <= parameters.xBody)
                %set thermal conductivity
                mesh.elem(i).k = parameters.ksubcut; 
            end

            %diffusion coefficient with set parameters
            mesh.elem(i).D = (mesh.elem(i).k)/((mesh.elem(i).rho * mesh.elem(i).c));

            %reaction coefficient with set parameters
            mesh.elem(i).lambda = -((mesh.elem(i).G * mesh.elem(i).rhob * mesh.elem(i).cb)/(mesh.elem(i).rho * mesh.elem(i).c));
            
            %source term with set parameters
            mesh.elem(i).f = ((mesh.elem(i).G * mesh.elem(i).rhob * mesh.elem(i).cb * mesh.elem(i).Tb)/(mesh.elem(i).rho * mesh.elem(i).c));
        end
    end
end


    
            
        
