function [ rhs ] = externalSource( x, t, xEnd, numberOfLayers,...
    tEnd, numberOfTimeStepsPerLayer, refinementDepth, sourceValue )
%EXTERNALSOURCE returns the extrnal heat source on the doamin
% for a given set of cordinates x at the
% given time t, the heat sorce is applied only at the first time step for
% each layer.

deltaT_layer = tEnd /numberOfLayers;
deltaT = tEnd /(numberOfLayers*numberOfTimeStepsPerLayer);

timeStep_layer = ceil(t/deltaT/numberOfTimeStepsPerLayer);

deltaX = xEnd/numberOfLayers;

boundaryConditionPoint = deltaX * timeStep_layer;
rhs = zeros(numel(x), 1);

% deltaRefined = deltaX / 2^(refinementDepth);

timeToApplyTheSource = (timeStep_layer-1)*deltaT_layer + deltaT;

for i=1:numel(x)
   if (x <=  boundaryConditionPoint && x >= (boundaryConditionPoint-deltaX)) ...
           && t <= timeToApplyTheSource
       rhs(i) = sourceValue;
   else
       rhs(i) = 0.0;
   end
end

end

