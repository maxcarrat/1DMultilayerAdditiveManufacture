function [ k ] = steelThermalConductivity( x, t, xEnd, numberOfLayers,...
    tEnd, numberOfTimeStepsPerLayer, T)
%STEELTHERMALCONDUCTIVITY returns the steel thermal conductivity for the
%powder and the solid part.

deltaT = tEnd /(numberOfLayers*numberOfTimeStepsPerLayer);
timeStep_layer = ceil(t/deltaT/numberOfTimeStepsPerLayer);
deltaX = xEnd/numberOfLayers;
boundaryConditionPoint = deltaX * timeStep_layer;

k = zeros(numel(T), 1);

for i=1:numel(T)
     k(i) = 54.0 - 26.7 * tanh( ( T(i) - 800.0 ) / 800.0 + 1.0 ); % [W/(m K)]
end

end


