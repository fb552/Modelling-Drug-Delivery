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

            %diffusion coefficient of current element
            mesh.elem(i).D = parameters.D;
            %reaction coefficient of current element
            mesh.elem(i).lambda = parameters.lambda;
            %source term of current element
            mesh.elem(i).f = parameters.f;

        elseif parameters.selection == '2'
            %position on the mesh of current element
            position = mesh.elem(i).x(1);    

            % Set parameters values based on current position
            %epidermis layer
            if(position >= 0 && position < parameters.Xe)
                %diffusion coefficient and no blood flow
                mesh.elem(i).D = parameters.De;
                mesh.elem(i).beta = parameters.betaE;
                mesh.elem(i).gamma = parameters.gammaE;

            %dermis layer
            elseif (position >= parameters.Xe && position < parameters.Xd)
                %diffusion coefficient
                mesh.elem(i).D = parameters.Dd;
                mesh.elem(i).beta = parameters.betaD;
                mesh.elem(i).gamma = parameters.gammaD;

            %sub-cutaneous layer
            elseif (position >= parameters.Xd && position <= parameters.Xb)
                %diffusion coefficient
                mesh.elem(i).D = parameters.Db;
                mesh.elem(i).beta = parameters.betaB;
                mesh.elem(i).gamma = parameters.gammaB;
            end

            %sum of the reaction terms
            mesh.elem(i).lambda = -(mesh.elem(i).beta + mesh.elem(i).gamma);
            %no source term
            mesh.elem(i).f = 0;
        end
    end
end