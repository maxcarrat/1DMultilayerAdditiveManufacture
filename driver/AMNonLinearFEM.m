%% 1D Multi-layer Non-Linear Additive Manufacturing
% evaluate a 1D poisson transient non-linear problem using h-FEM

clear all;
clc;

writerObj = VideoWriter('AMNonLinearHFEM.avi');
writerObj.Quality = 100;
writerObj.FrameRate = 5;
open(writerObj);

%% Problem SetUp
% Define the problem parameter, the boundary conditions
% and the discretization.

rho = 7500.0;                                               % density [kg/m^3]
c = 600.0;                                                  % specific heat [J/(kg°C)]
% k = @(T) 54.0 - 26.7 * tanh( ( T - 800.0 ) / 800.0 + 1.0 ); % thermal conductivity [W/(m°C)]
T0 = 20.0;                                                  % Initial temperature [°C]
heatCapacity= rho*c;                                        % heat capacity [kJ / kg °C]
Tsource = 2000.0;                                           % source temperature [°C]

tEnd = 45;                                                 % total time [sec]
xEnd = 0.01;                                               % length of the bar 0.001 [m]

dirichletLeftBC = @(t) T0;
dirichletRightBCValue = T0 + Tsource;
nuemannRightBC = 0.0e+00;
bodySource = 0.0e+13;

tolerance = 1.0e-03;
maxIteration = 20;

for depth = 1:5
    % maxTrainingTimeSteps = timeSteps*0.5;
    % relErrorEnergyNorm = zeros(maxTrainingTimeSteps-3,1);
    % modes = zeros(maxTrainingTimeSteps-3,1);
    % tainingVector = linspace(3,maxTrainingTimeSteps, maxTrainingTimeSteps-2);
    
    numberOfLayers = 200;
    trainingTimeSteps = numberOfLayers;
    numberOfTimeStepsPerLayer = 10;     % total time per layer 0.225 sec
    numberOfHeatingTimeSteps = 4;       % heating laser time per layer 0.090 sec

    integrationOrder = 2;
    
    % for trainingTimeSteps = 3:maxTrainingTimeSteps
    
    numberOfElementsInX = numberOfLayers;
    timeSteps = 4;
    refinementDepth = depth;
    numberOfPODModes = 0;
    
    numberOfRefinedElementsToBeKept = 1;
 
    dirichletRightBC = @(t) heatingCoolingBoundary(t, numberOfLayers, tEnd,...
        numberOfTimeStepsPerLayer, numberOfHeatingTimeSteps, dirichletRightBCValue);
    
    bodyLoad = @(x, t) externalSource(x, t, xEnd, numberOfLayers, tEnd,...
        numberOfTimeStepsPerLayer, refinementDepth, bodySource);
    
    k = @(x, t, T) steelThermalConductivity(x, t, xEnd, numberOfLayers, tEnd,...
        numberOfTimeStepsPerLayer, T);  % thermal conductivity [W/(m K)]
    
    
    t = linspace(0, tEnd, numberOfTimeStepsPerLayer*numberOfLayers + 1);        % time discretization
    x = linspace(0.0, xEnd, numberOfElementsInX + 1);                           % spatial discretization X
    layerCoords = linspace(0.0, xEnd, numberOfElementsInX + 1);                 % layer spatial discretization X
    
    x_ref = linspace(0.0, 1.0, 2^refinementDepth + 1);                          % refinement discretization X_ref
    x_PostProcess = linspace(0.0, xEnd, (numberOfElementsInX +1) * 2^8);          % post-processing coordinates
    
    [X, T] = meshgrid(x_PostProcess, t);
    [X_ref, T_ref] = meshgrid(x_ref, t);
    
    
    %% Analysis
        [temperatureSolution, heatFlux, internalEnergy, CPUTime] = nonLinearBackwardEulerGaussIntegrationSolver(x, x_PostProcess, bodyLoad, T0,...
            dirichletLeftBC, dirichletRightBC, nuemannRightBC, k, heatCapacity, t, tolerance, maxIteration,numberOfRefinedElementsToBeKept,...
            refinementDepth, trainingTimeSteps, numberOfTimeStepsPerLayer,...
            numberOfLayers, numberOfPODModes, integrationOrder);
%     [temperatureSolution, heatFlux, internalEnergy] = nonLinearBackwardEulerNoCoarseningGaussIntegrationSolver(x, x_PostProcess, bodyLoad, T0,...
%         dirichletLeftBC, dirichletRightBC, nuemannRightBC, k, heatCapacity, t, tolerance, maxIteration, ...
%         refinementDepth, numberOfTimeStepsPerLayer, numberOfLayers, integrationOrder);
    
    %% Post-Process
    
    figure(depth+6)
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
        xlim(axes1,[0 xEnd]);
        ylim(axes1,[-100 3000]);
        
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    
    % Create ylabel
    
    
    hold off
    
    % Write results to a file
    formatSpec = 'myFEMNonLinearResultsFile_Coarsening_%d.txt';
    filename = sprintf(formatSpec,depth);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:size(temperatureSolution, 1)
        for j=1:numberOfTimeStepsPerLayer
            fprintf(resultFile, '%d, %\t', temperatureSolution(i,end-(j-1)));
        end
        fprintf(resultFile, '\n');
    end
    fclose(resultFile);
    
    % Write CPU time to a file
    formatSpec = 'myFEMNonLinearTimeFile_Coarsening_%d.txt';
    filename = sprintf(formatSpec,depth);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:numel(CPUTime)
        fprintf(resultFile, '%d, %\t', CPUTime(i));
    end
    
    fclose(resultFile);
    
    
    
end % depth loop


figure(10101)

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
    xlim(axes1,[0 xEnd]);
    ylim(axes1,[-100 3000]);
    
    box(axes1,'on');
    % Set the remaining axes properties
    set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
    
    grid on
    drawnow
    F(i) = getframe;
    writeVideo(writerObj, getframe(gcf, [ 0 0 500 400 ]));
end

close(writerObj);