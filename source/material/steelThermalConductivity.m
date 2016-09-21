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
%    if (x <=  boundaryConditionPoint && x >= (boundaryConditionPoint-deltaX)) 
%        k(i) = 54.0 - 26.7 * tanh( ( T - 800.0 ) / 800.0 + 1.0 ); % [W/(m K)]
%    else
%        k(i) = 25;
%    end

%    if T(i) > 0 && T(i) <= 800
%        k(i) = 54.0 - 27.3 * (T(i) / 800); % [W/(m K)]
%    else
%        k(i) = 26.7;
%    end

     k(i) = 54.0 - 26.7 * tanh( ( T(i) - 800.0 ) / 800.0 + 1.0 ); % [W/(m K)]

   
end

end


