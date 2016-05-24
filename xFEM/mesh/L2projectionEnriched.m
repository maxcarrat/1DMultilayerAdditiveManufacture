function [ projectedTemperature ] = L2projectionEnriched(previousTemperature, updatedMesh, previousMesh, modes)
%L2PROJECTIONENRICHED project the previous solution onto the updated mesh at the
%new time step considering the enrichement modes from POD.
%   previousTemperature = temeprature distribution of the previous mesh
%   updatedMesh = actual mesh configuration
%   previousmesh = previous mesh configuration
%   modes = number of enrichment modes

projectedTemperature = zeros(size(updatedMesh,2), 1);

for i=1:size(previousMesh,2)
    X = previousMesh(i);
    for j=1:size(updatedMesh,2)
        x = updatedMesh(j);
        if(X == x)
            projectedTemperature(j) = previousTemperature(i);
        end
    end
end

end

