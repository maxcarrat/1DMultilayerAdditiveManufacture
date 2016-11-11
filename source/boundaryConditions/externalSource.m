function [ rhs ] = externalSource( x, t, xEnd, numberOfLayers,...
    tEnd, numberOfTimeStepsPerLayer, numberOfHeatingTimeSteps, sourceValue )
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

heatingTime = (timeStep_layer-1)*deltaT_layer +...
    deltaT * numberOfHeatingTimeSteps;

for i=1:numel(x)
   if (x ==  boundaryConditionPoint) ...
           && t <= heatingTime
       rhs(i) = sourceValue;
   else
       rhs(i) = 0.0;
   end
end

end

