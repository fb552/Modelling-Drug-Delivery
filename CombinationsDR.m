function [parameters] = CombinationsDR(I,parameters)

    % this function should only be called in material 2 testing
    parameters.selection = '2';  %material for Part 2

    switch I
        case 1
            parameters.De = 4e-5;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-5;        %Diffusion coefficient dermis
            parameters.Db = 2e-5;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.01;     %Reaction coefficient blood dermis
            parameters.betaB = 0.01;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.02;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.02;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.02;    %Reaction coefficient degradation sub-cutaneous
        case 2
            parameters.De = 4e-7;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-7;        %Diffusion coefficient dermis
            parameters.Db = 2e-7;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.01;     %Reaction coefficient blood dermis
            parameters.betaB = 0.01;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.02;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.02;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.02;    %Reaction coefficient degradation sub-cutaneous
        case 3
            parameters.De = 4e-6;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-6;        %Diffusion coefficient dermis
            parameters.Db = 2e-6;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.01;     %Reaction coefficient blood dermis
            parameters.betaB = 0.01;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.2;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.2;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.2;    %Reaction coefficient degradation sub-cutaneous
        case 4
            parameters.De = 4e-6;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-6;        %Diffusion coefficient dermis
            parameters.Db = 2e-6;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.001;     %Reaction coefficient blood dermis
            parameters.betaB = 0.001;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.2;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.2;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.2;    %Reaction coefficient degradation sub-cutaneous
        case 5
            parameters.De = 4e-6;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-6;        %Diffusion coefficient dermis
            parameters.Db = 2e-6;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.1;     %Reaction coefficient blood dermis
            parameters.betaB = 0.1;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.002;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.002;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.002;    %Reaction coefficient degradation sub-cutaneous
        case 6
            parameters.De = 4e-6;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-6;        %Diffusion coefficient dermis
            parameters.Db = 2e-6;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.1;     %Reaction coefficient blood dermis
            parameters.betaB = 0.1;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.02;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.02;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.02;    %Reaction coefficient degradation sub-cutaneous
        case 7
            parameters.De = 4e-6;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-6;        %Diffusion coefficient dermis
            parameters.Db = 2e-6;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.001;     %Reaction coefficient blood dermis
            parameters.betaB = 0.001;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.02;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.02;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.02;    %Reaction coefficient degradation sub-cutaneous
        case 8
            parameters.De = 4e-6;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-6;        %Diffusion coefficient dermis
            parameters.Db = 2e-6;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.01;     %Reaction coefficient blood dermis
            parameters.betaB = 0.01;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.002;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.002;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.002;    %Reaction coefficient degradation sub-cutaneous

    end
end