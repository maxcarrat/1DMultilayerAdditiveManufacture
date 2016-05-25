function [ projectedTemperature ] = L2projection(previousTemperature, updatedMesh, previousMesh)
%L2PROJECTION project the previous solution onto the updated mesh at the
%new time step
%   previousTemperature = temeprature distribution of the previous mesh
%   updatedMesh = actual mesh configuration
%   previousmesh = previous mesh configuration

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

