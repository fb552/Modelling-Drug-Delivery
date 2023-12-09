function [parameters] = CombinationsDR(I,parameters)
%Stores a series of parameters combinations for diffusion and reaction coefficients.
% Changes the parameters data structure to have different coefficients for
% diffusion, beta and gamma. Entera a combination value in range [1:1:8]
%
% Input:
%  I : combination number
%  parameters : current parameters data structure
% Return:
%  parameters : new parameters data structure
%
%Francesco Berteau (fb552) - November 2023

    % this function should only be called in material 2 testing
    parameters.selection = '2';  %material for Part 2

    switch I
        case 1  % High diffusion
            parameters.De = 4e-5;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-5;        %Diffusion coefficient dermis
            parameters.Db = 2e-5;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.01;     %Reaction coefficient blood dermis
            parameters.betaB = 0.01;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.02;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.02;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.02;    %Reaction coefficient degradation sub-cutaneous
        case 2  % Low diffusion
            parameters.De = 4e-7;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-7;        %Diffusion coefficient dermis
            parameters.Db = 2e-7;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.01;     %Reaction coefficient blood dermis
            parameters.betaB = 0.01;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.02;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.02;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.02;    %Reaction coefficient degradation sub-cutaneous
        case 3  % high degradation
            parameters.De = 4e-6;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-6;        %Diffusion coefficient dermis
            parameters.Db = 2e-6;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.01;     %Reaction coefficient blood dermis
            parameters.betaB = 0.01;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.2;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.2;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.2;    %Reaction coefficient degradation sub-cutaneous
        case 4  % high degradation, low blood
            parameters.De = 4e-6;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-6;        %Diffusion coefficient dermis
            parameters.Db = 2e-6;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.001;     %Reaction coefficient blood dermis
            parameters.betaB = 0.001;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.2;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.2;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.2;    %Reaction coefficient degradation sub-cutaneous
        case 5  % low degradation, high blood
            parameters.De = 4e-6;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-6;        %Diffusion coefficient dermis
            parameters.Db = 2e-6;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.1;     %Reaction coefficient blood dermis
            parameters.betaB = 0.1;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.002;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.002;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.002;    %Reaction coefficient degradation sub-cutaneous
        case 6  % high blood
            parameters.De = 4e-6;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-6;        %Diffusion coefficient dermis
            parameters.Db = 2e-6;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.1;     %Reaction coefficient blood dermis
            parameters.betaB = 0.1;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.02;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.02;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.02;    %Reaction coefficient degradation sub-cutaneous
        case 7  % low blood
            parameters.De = 4e-6;        %Diffusion coefficient epidermis
            parameters.Dd = 5e-6;        %Diffusion coefficient dermis
            parameters.Db = 2e-6;        %Diffusion coefficient sub-cutaneous
            parameters.betaE = 0;        %Reaction coefficient blood epidermis
            parameters.betaD = 0.001;     %Reaction coefficient blood dermis
            parameters.betaB = 0.001;     %Reaction coefficient blood sub-cutaneous
            parameters.gammaE = 0.02;    %Reaction coefficient degradation epidermis
            parameters.gammaD = 0.02;    %Reaction coefficient degradation dermis
            parameters.gammaB = 0.02;    %Reaction coefficient degradation sub-cutaneous
        case 8  % low degradation
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