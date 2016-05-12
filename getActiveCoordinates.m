function [ activeCoords ] = getActiveCoordinates( coords, time, timeSteps )
%[activeCoords] = GETACTIVECOORDINATES: get the active coords at a given
%time Step

offset = size(coords,2) - timeSteps;
activeCoords = coords(1:(time+offset));

end

