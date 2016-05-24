function [ activeTemp ] = getActiveTemperatureSolution( temperatureSolutions, time, timeSteps )
%[activetemp] = GETACTIVETEMPERATURESOLUTION: get the active temperature
%solutions at a given time Step

offset = size(temperatureSolutions,1) - timeSteps;
activeTemp = temperatureSolutions(1:(time+offset),time-1);

end

