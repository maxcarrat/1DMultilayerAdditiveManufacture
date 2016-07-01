%% 1D Multi-layer Non-Linear Additive Manufacturing with enriched local functions
% evaluate a 1D poisson transient non-linear problem using X-FEM

clear all;
clc;

%% Problem SetUp
% Define the problem parameter, the boundary conditions
% and the discretization.

rho = 7500.0;                                               % density [kg/m^3]
c = 600.0;                                                  % specific heat [J/(kg°C)]
T0 = 20.0;                                                  % Initial temperature [°C]
heatCapacity= rho*c;                                        % heat capacity [kJ / kg °C]
Tsource = 2000.0;                                           % source temperature [°C]

tEnd = 4.5;                                                 % total time [sec]
xEnd = 0.01;                                                % length of the bar [m]

dirichletLeftBC = @(t) T0;
dirichletRightBCValue =  T0 + Tsource;
nuemannRightBC = 0.0e+00;
bodySource = 0.0e+13;

tolerance = 1.0e-03;
maxIteration = 20;

for modes = 1:5
    % maxTrainingTimeSteps = timeSteps*0.5;
    % relErrorEnergyNorm = zeros(maxTrainingTimeSteps-3,1);
    % modes = zeros(maxTrainingTimeSteps-3,1);
    % tainingVector = linspace(3,maxTrainingTimeSteps, maxTrainingTimeSteps-2);
    
    myVideoSpec = 'AdditiveManufacturingNonLinearPODBasis_%d.avi';
    myVideo = sprintf(myVideoSpec,modes);
    writerObj = VideoWriter(myVideo);
    writerObj.Quality = 100;
    writerObj.FrameRate = 10;
    open(writerObj);
    
    numberOfLayers = 20;
    trainingTimeSteps = 10;
    numberOfTimeStepsPerLayer = 10;
    numberOfHeatingTimeSteps = 2;
    refinementDepth = 4;
    
    integrationOrder = 2;
    integrationModesOrder = (modes - 1) ^ 2  + 1 + modes;
    
    numberOfRefinedElementsToBeKept = 1;
    
    % for trainingTimeSteps = 3:maxTrainingTimeSteps
    
    numberOfElementsInX = numberOfLayers;
    timeSteps = 4;
    numberOfPODModes = modes;
    
    dirichletRightBC = @(t) heatingCoolingBoundary(t, numberOfLayers, tEnd,...
        numberOfTimeStepsPerLayer, numberOfHeatingTimeSteps, dirichletRightBCValue);
    
    bodyLoad = @(x, t) externalSource(x, t, xEnd, numberOfLayers, tEnd,...
        numberOfTimeStepsPerLayer, refinementDepth, bodySource);
    
    k = @(x, t, T) steelThermalConductivity(x, t, xEnd, numberOfLayers, tEnd,...
        numberOfTimeStepsPerLayer, T);  % thermal conductivity [W/(m K)]
    
    t = linspace(0, tEnd, numberOfTimeStepsPerLayer*numberOfLayers + 1);                        % time discretization
    x = linspace(0.0, xEnd, numberOfElementsInX + 1);                                           % spatial discretization X
    layerCoords = linspace(0.0, xEnd, numberOfElementsInX + 1);                                 % layer spatial discretization X
    
    x_PostProcess = linspace(0.0, xEnd, (numberOfElementsInX +1) * 2.^8);                       % post-processing coordinates
    
    [X, T] = meshgrid(x_PostProcess, t);
    
    
    %% Analysis
    [temperatureSolution, heatFlux, internalEnergy, CPUTime] = nonLinearBackwardEulerGaussIntegrationSolver(x, x_PostProcess, bodyLoad, T0,...
        dirichletLeftBC, dirichletRightBC, nuemannRightBC, k, heatCapacity, t, tolerance, maxIteration, numberOfRefinedElementsToBeKept,...
        refinementDepth, trainingTimeSteps, numberOfTimeStepsPerLayer,...
        numberOfLayers, numberOfPODModes, integrationOrder, integrationModesOrder);
    
    %% Post-Process
    figure(modes + 111111)
    
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
        ylim(axes1,[-100 2500]);
        
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    
    hold off
    
    % Write results to a file
    formatSpec = 'myPODXFEMNonLinearResultsFile_%d.txt';
    filename = sprintf(formatSpec,modes);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:size(temperatureSolution, 1)
        for j=1:numberOfTimeStepsPerLayer
            fprintf(resultFile, '%d, %\t', temperatureSolution(i,end-(j-1)));
        end
        fprintf(resultFile, '\n');
    end
    fclose(resultFile);
    
    % Write CPU time to a file
    formatSpec = 'myPODXFEMNonLinearTimeFile_%d.txt';
    filename = sprintf(formatSpec,modes);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:numel(CPUTime)
        fprintf(resultFile, '%d, %\t', CPUTime(i));
    end
    
    fclose(resultFile);
    
    figure(modes + 11111111)
    % Create axes
    axes1 = axes;
    
    F(size(t,2)) = struct('cdata',[],'colormap',[]);
    
    for i=1:size(t,2)
        plot(x_PostProcess,temperatureSolution(:,i)')
        
        % Title
        myTitleSpec = 'Temperature evolution of the bar using %d. modes';
        myTitle = sprintf(myTitleSpec,modes);
        title( myTitle, 'Interpreter','latex');
        
        % Create ylabel
        ylabel('\fontname{Latin Modern Math} Temperature [°C]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes1,[0 xEnd]);
        ylim(axes1,[-100 2500]);
        
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        grid on
        drawnow
        F(i) = getframe;
        writeVideo(writerObj, getframe(gcf, [ 0 0 500 400 ]));
    end
    
    close(writerObj);
    
    
    
end % modes loop


