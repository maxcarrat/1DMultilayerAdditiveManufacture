function [ projectedBaseTemperature, projectedOverlayTemperature ] = ...
    L2MultiscaleProjection( previousBaseTemperature, previousOverlayTemperature,...
    baseProblem, overlayProblem, integrationOrder, integrationModalOrder,...
    layerLength, initialTemperature, eXtendedFlag, previousBaseProblem)
%L2MULTISCALEPROJECTION project the previous solution onto the updated mesh 
%at the new time step.
%Input:
%previousBaseTemperature = temperature coefficients of the last converged
%solution on the base mesh
%previousOverlayTemperature = temperature coefficients of the last converged
%solution on the overlay mesh
%baseProblem = poisson Problem struct of the base mesh
%overlayProblem = poisson Problem struct of the overlay mesh
%integrationOrder = number of quadrature points for Gauss integration
%integrationModalOrder = number of quadrature points for the enriched
%elements
%layerLength = length of the layer
%initialTemperature = initial temperature of the powder
%eXtendedFlag = true if ROM-phase, false if training-phase
%previousBaseProblem = poisson problem struct of the previous layer base
%mesh
%Output:
%projectedTemperature = L2 projection of the last converged solution onto
%the new domain

if strcmp(eXtendedFlag,'true')
    
    [ M, f ] = assemblyL2ProjectionMultiscalePODXIGA( previousBaseTemperature,...
        previousOverlayTemperature, baseProblem,...
        overlayProblem, integrationOrder, integrationModalOrder, layerLength,...
        initialTemperature, previousBaseProblem );
    
elseif strcmp(eXtendedFlag,'transition')
    
    [ M, f ] = assemblyL2ProjectionMutiscalehIGAOntoPODXIGA( previousBaseTemperature,...
        previousOverlayTemperature, baseProblem,...
        overlayProblem, integrationOrder, integrationModalOrder, layerLength,...
        initialTemperature, previousBaseProblem );
else
    
    [ M, f ] = assemblyL2ProjectionMultiscalehIGA( previousBaseTemperature,...
        previousOverlayTemperature, baseProblem,...
        overlayProblem, previousBaseProblem, integrationOrder, initialTemperature );
end

projectedTemperature = M\f;

projectedBaseTemperature = projectedTemperature(1:baseProblem.gdof);
projectedOverlayTemperature = projectedTemperature(1+baseProblem.gdof:end);

end

