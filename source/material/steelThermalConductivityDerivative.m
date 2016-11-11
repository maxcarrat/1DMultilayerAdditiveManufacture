function [ kDerivative ] = steelThermalConductivityDerivative( x, t, xEnd, numberOfLayers,...
    tEnd, numberOfTimeStepsPerLayer, T )
%STEELTHERMALCONDUCTIVITYDERIVATIVE returns the steel thermal conductivity derivative for the
%powder and the solid part.

kDerivative = zeros(numel(T), 1);

for i=1:numel(T)
     kDerivative(i) = - 0.033375 * sech(0.00125*T(i))^2; % [W/(m K ^2)]
end

end

