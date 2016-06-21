function [ activeTemp ] = getActiveTemperatureSolution( temperatureSolutions, initialTemperature,...
time, timeSteps )
%[activetemp] = GETACTIVETEMPERATURESOLUTION: get the active temperature
%solutions at a given time Step

offset = size(temperatureSolutions,1) - timeSteps;
activeTemp = zeros(time+offset, 1);
activeTemp(:) = initialTemperature;

for i=1:time+offset-1
    activeTemp(i) = temperatureSolutions(i,time-1);
end

end

