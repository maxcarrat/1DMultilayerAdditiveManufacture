function [ projectedTemperature ] = L2projectionIGA( previousTemperature, problem, integrationOrder,...
    integrationModalOrder, layerLength, initialTemperature, eXtendedFlag, previousProblem)
%L2PROJECTIONLINEARDISTRIBUTION project the previous solution onto the updated mesh at the
%new time step
%Input:
%previousTemperature = temperature coefficients of the last converged
%solution
%problem = poisson Problem struct
%integrationOrder = number of quadrature points for Gauss integration
%integrationModalOrder = number of quadrature points for the enriched
%elements
%layerLength = length of the layer
%initialTemperature = initial temperature of the powder
%eXtendedFlag = true if ROM-phase, false if training-phase
%previousProblem = poisson problem struct of the previous layer
%Output:
%projectedTemperature = L2 projection of the last converged solution onto
%the new domain

if strcmp(eXtendedFlag,'true')
    [ M, f ] = assemblyL2ProjectionXIGAMatrixAndVector( previousTemperature, problem, integrationOrder,...
        integrationModalOrder, layerLength, initialTemperature, previousProblem );
elseif strcmp(eXtendedFlag,'transition')
    [ M, f ] = assemblyL2ProjectionIGAOntoXIGAMatrixAndVector( previousTemperature, problem, integrationOrder,...
        integrationModalOrder, layerLength, initialTemperature, previousProblem );
else
    [ M, f ] = assemblyL2ProjectionIGAMatrixAndVector( previousTemperature, problem, integrationOrder,...
        layerLength, initialTemperature );
end

projectedTemperature = M\f;

end
