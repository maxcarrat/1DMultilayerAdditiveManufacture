%% 1D Multi-layer Linear Additive Manufacturing Overkilled
% evaluate a 1D poisson transient linear problem using h-FEM

clear all;
clc;

writerObj = VideoWriter('AdditiveManufacturingOverkilled.avi');
writerObj.Quality = 100;
writerObj.FrameRate = 5;
open(writerObj);

%% Problem SetUp
% Define the problem parameter, the boundary conditions
% and the discretization.

rho = 7500.0;                                               % density [kg/m^3]
c = 600.0;                                                  % specific heat [J/(kg°C)]
k = 56.0;                                                   % thermal conductivity [W/(m°C)]
T0 = 20.0;                                                  % Initial temperature [°C]
heatCapacity= rho*c;                                        % heat capacity [kJ / kg °C]
Tsource = 2000.0;                                           % source temperature [°C]

tEnd = 0.1;                                                 % total time [sec]
xEnd = 0.001;                                               % length of the bar 0.001 [m]

dirichletLeftBC = @(t) T0;
dirichletRightBC = @(t) T0 + Tsource;
nuemannRightBC = 0.0e+08;
bodySource = 0.0e+12;

% maxTrainingTimeSteps = timeSteps*0.5;
% relErrorEnergyNorm = zeros(maxTrainingTimeSteps-3,1);
% modes = zeros(maxTrainingTimeSteps-3,1);
% tainingVector = linspace(3,maxTrainingTimeSteps, maxTrainingTimeSteps-2);

numberOfLayers = 20;
trainingTimeSteps = numberOfLayers;
numberOfTimeStepsPerLayer = 5;

integrationOrder = 1;

% for trainingTimeSteps = 3:maxTrainingTimeSteps

numberOfElementsInX = numberOfLayers;
timeSteps = 4;
refinementDepth = 6;
numberOfPODModes = 10;

bodyLoad = @(x, t) externalSource(x, t, xEnd, numberOfLayers, tEnd,...
    numberOfTimeStepsPerLayer, refinementDepth, bodySource);


t = linspace(0, tEnd, numberOfTimeStepsPerLayer*numberOfLayers + 1);        % time discretization
x = linspace(0.0, xEnd, numberOfElementsInX + 1);                           % spatial discretization X
layerCoords = linspace(0.0, xEnd, numberOfElementsInX + 1);                 % layer spatial discretization X

x_ref = linspace(0.0, 1.0, 2^refinementDepth + 1);                          % refinement discretization X_ref
x_PostProcess = linspace(0.0, xEnd, (numberOfElementsInX +1) * 2.^8);          % post-processing coordinates

[X, T] = meshgrid(x_PostProcess, t);
[X_ref, T_ref] = meshgrid(x_ref, t);


%% Analysis
[temperatureSolution, heatFlux, internalEnergy] = backwardEulerNoCoarseningGaussIntegrationSolver(x, x_PostProcess, bodyLoad, T0,...
    dirichletLeftBC, dirichletRightBC, nuemannRightBC, k, heatCapacity, t, refinementDepth, numberOfTimeStepsPerLayer, numberOfLayers, integrationOrder);

%% Post-Process

figure(0002)
% Create axes
axes1 = axes;

for i=1:size(temperatureSolution, 2)
    
    % Create multiple lines using matrix input to plot
    plot(x_PostProcess,temperatureSolution(:,i)');
    
    % Title
    title('Temperature distibution at each time step', 'Interpreter','latex');
    
    % Axis labels
    ylabel('\fontname{Latin Modern Math} Temperature [°C]');
    xlabel('\fontname{Latin Modern Math} length [m]');
    
    % Axis limit
    xlim(axes1,[0 0.001]);
    ylim(axes1,[-100 3000]);
    
    box(axes1,'on');
    % Set the remaining axes properties
    set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
    
    hold on
    
end

hold off

% Write results to a file
filename = 'myReferenceResultsFile.txt';            % Name of the file
resultFile = fopen(filename, 'wt');                 % Open for writing
for i=1:size(temperatureSolution, 1)
    for j=1:numberOfTimeStepsPerLayer
        fprintf(resultFile, '%d, %\t', temperatureSolution(i,end-(j-1)));
    end
    fprintf(resultFile, '\n');
end
fclose(resultFile);

figure(0004)
% Create axes
axes1 = axes;

F(size(t,2)) = struct('cdata',[],'colormap',[]);

for i=1:size(t,2)
    plot(x_PostProcess,temperatureSolution(:,i)')
    
    % Title
    title('Temperature evolution of the bar', 'Interpreter','latex');
    
    % Create ylabel
    ylabel('\fontname{Latin Modern Math} Temperature [°C]');
    xlabel('\fontname{Latin Modern Math} length [m]');
    
    % Axis limit
    xlim(axes1,[0 0.001]);
    ylim(axes1,[-100 3000]);
    
    box(axes1,'on');
    % Set the remaining axes properties
    set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
    
    grid on
    drawnow
    F(i) = getframe;
    writeVideo(writerObj, getframe(gcf, [ 0 0 560 420 ]));
end

close(writerObj);