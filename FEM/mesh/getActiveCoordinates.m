function [ activeCoords ] = getActiveCoordinates( coords, time, timeSteps )
%[activeCoords] = GETACTIVECOORDINATES: get the active coords at a given
%time Step
%   coords = 1D mesh
%   time = actual time step
%   timeSteps = number of time steps

offset = size(coords,2) - timeSteps;
activeCoords = coords(1:((time)+offset));

end

