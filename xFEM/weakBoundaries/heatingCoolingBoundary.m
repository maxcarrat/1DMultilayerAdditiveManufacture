function [ T ] = heatingCoolingBoundary( t, numberOfLayers,...
    tEnd, numberOfTimeStepsPerLayer, numberOfHeatingTimeSteps,...
     boundaryValue )
%HEATINGCOOLINGBOUNDARY Returns a boundary condition with a given
%temperature in the first time interval and no boundaries in the second
%time interval

deltaT_layer = tEnd /numberOfLayers;
deltaT = tEnd /(numberOfLayers * numberOfTimeStepsPerLayer);

timeStep_layer = ceil(t/deltaT/numberOfTimeStepsPerLayer);

T = [];

heatingTime = (timeStep_layer-1)*deltaT_layer +...
    deltaT * numberOfHeatingTimeSteps;

if  t <= heatingTime
    T = [T boundaryValue];
end

end

