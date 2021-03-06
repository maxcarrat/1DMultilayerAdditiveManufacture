%% 1D Multi-layer Non-Linear Multi-phase Additive Manufacturing with enriched local functions
% evaluate a 1D poisson transient non-linear problem using X-FEM

clear all;
clc;

mkdir('XFEMTest/Temperature');
mkdir('XFEMTest/Fluxes');
mkdir('XFEMTest/Time');
mkdir('XFEMTest/Dofs');
mkdir('XFEMTest/Files4Pics');


%% Problem SetUp
% Define the problem parameter, the boundary conditions
% and the discretization.

rho = 7500.0;                                               % density [kg/m^3]
c = 600.0;                                                  % specific heat [J/(kg°C)]
T0 = 20.0;                                                  % Initial temperature [°C]
entalphyJump = 261e+03;                                     % entalphy [J/Kg]
Tsource = 1500.0;                                           % source temperature [°C]
Tmelt = 1450.0;                                             % melting temperature [°C]

tEnd = 4.5;                                                 % total time [sec]
xEnd = 0.001;                                               % length of the bar [m]

dirichletLeftBC = @(t) T0;
dirichletRightBCValue =  T0 + Tsource;
nuemannRightBC = 0.0e+04;
bodySource = 0.0e+13;

tolerance = 1.0e-05;
maxIteration = 50;

initialNumberOfModes = 3;
maxNumberOfModes = 3;
DOFs = zeros(maxNumberOfModes, 1);
% modes = 1;

for modes = initialNumberOfModes:maxNumberOfModes
    % for refinementDepth = 1:8
    
    numberOfLayers = 20;
    trainingTimeSteps = 5;
%     trainingTimeSteps = 20;
    
    numberOfTimeStepsPerLayer = 10;     % total time per layer 0.225 sec
    numberOfHeatingTimeSteps = 4;       % heating laser time per layer 0.090 sec
    
    refinementDepth = 8;
    PODRefinementDepth = 0;
    
    integrationOrder = 2;
    integrationModesOrder = modes + 30;  %(modes - 1) ^ 2  + 1 + modes;
    
    numberOfRefinedElementsToBeKept = 1;
    
    numberOfElementsInX = numberOfLayers;
    timeSteps = 4;
    numberOfPODModes = modes;
    
    dirichletRightBC = @(t) heatingCoolingBoundary(t, numberOfLayers, tEnd,...
        numberOfTimeStepsPerLayer, numberOfHeatingTimeSteps, dirichletRightBCValue);
    
    bodyLoad = @(x, t) externalSource(x, t, xEnd, numberOfLayers, tEnd,...
        numberOfTimeStepsPerLayer, refinementDepth, bodySource);
    
    %     k = @(x, t, T) 27.0;
    k = @(x, t, T)  steelThermalConductivity(x, t, xEnd, numberOfLayers, tEnd,...
        numberOfTimeStepsPerLayer, T);                                                          % thermal conductivity [W/(m K)]
    
    %     heatCapacity= @(x, t, T) 7500 * 600;
    heatCapacity= @(x, T, T_last)  capacityPhaseTransition(x, T, T_last, c, rho,...
        entalphyJump, Tsource, Tmelt);                              % heat capacity [kJ / kg °C]
    
    t = linspace(0, tEnd, numberOfTimeStepsPerLayer*numberOfLayers + 1);                        % time discretization
    x = linspace(0.0, xEnd, numberOfElementsInX + 1);                                           % spatial discretization X
    layerCoords = linspace(0.0, xEnd, numberOfElementsInX + 1);                                 % layer spatial discretization X
    
    %     x_PostProcess = linspace(0.0, xEnd, ( numberOfElementsInX + 1 ) * 2.^8);                    % post-processing coordinates
    x_PostProcess = linspace(0.0, xEnd, 6000);                    % post-processing coordinates
    [X, T] = meshgrid(x_PostProcess, t);
    
    %% Analysis
    [temperatureSolution, heatFlux, internalEnergy, CPUTime, DOFs(modes), ...
        localRefinedTemperatureSolutions, solutionReductionOperator] = multiPhaseBackwardEulerSolver(x, x_PostProcess, bodyLoad, T0,...
        dirichletLeftBC, dirichletRightBC, nuemannRightBC, k, heatCapacity, t, tolerance, maxIteration, numberOfRefinedElementsToBeKept,...
        refinementDepth, PODRefinementDepth, trainingTimeSteps, numberOfTimeStepsPerLayer,...
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
        ylim(axes1,[-100 2000]);
        
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        grid on
        hold on
        
    end
    
    hold off
    
    figure(modes+9999)
    % Create axes
    axes2 = axes;
    
    for i=1:size(heatFlux, 2)
        
        % Create multiple lines using matrix input to plot
        plot(x_PostProcess,heatFlux(:,i)');
        
        % Title
        title('Heat flux distibution at each time step', 'Interpreter','latex');
        
        % Axis labels
        ylabel('\fontname{Latin Modern Math} Heat flux [W/m^2]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes2,[0 xEnd]);
        ylim(axes2,[0.0e+01 4.0e+08]);
        
        box(axes2,'on');
        % Set the remaining axes properties
        set(axes2,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    
    % Create ylabel
    
    
    hold off
    
    % Write results to a file
    
    
    formatSpec = 'XFEMTest/Temperature/myXFEMMultiphaseResultsFile_%d.txt';
    filename = sprintf(formatSpec,refinementDepth);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:size(temperatureSolution, 1)
        fprintf(resultFile, '%d, %\t', x_PostProcess(i));
        for j=1:numberOfTimeStepsPerLayer
            fprintf(resultFile, '%d, %\t', temperatureSolution(i,end-(j-1)));
        end
        fprintf(resultFile, '\n');
    end
    fclose(resultFile);
    
    % Write fluxes to a file
    formatSpec = 'XFEMTest/Fluxes/myXFEMMultiphaseFluxesFile_%d.txt';
    filename = sprintf(formatSpec,refinementDepth);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:size(heatFlux, 1)
        fprintf(resultFile, '%d, %\t', x_PostProcess(i));
        for j=1:numberOfTimeStepsPerLayer
            fprintf(resultFile, '%d, %\t', heatFlux(i,end-(j-1)));
        end
        fprintf(resultFile, '\n');
    end
    fclose(resultFile);
    
    % Write CPU time to a file
    formatSpec = 'XFEMTest/Time/myXFEMMultiphaseTimeFile_%d.txt';
    filename = sprintf(formatSpec,refinementDepth);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:numel(CPUTime)
        fprintf(resultFile, '%d, %\t', i);
        fprintf(resultFile, '%d, %\n', CPUTime(i));
        fprintf(resultFile, '\n');
    end
    fclose(resultFile);
    
    %% Write files for presentation
    %Temperature Overlay mesh training phase
    formatSpec = 'XFEMTest/Files4Pics/TrainingPhaseOverlayTemperature_%d.txt';
    filename = sprintf(formatSpec,modes);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:size(localRefinedTemperatureSolutions, 1)
        fprintf(resultFile, '%d, %\t', i);
        for k = 1:trainingTimeSteps*numberOfTimeStepsPerLayer+1
            fprintf(resultFile, '%d, %\t', localRefinedTemperatureSolutions(i,k));
        end
        fprintf(resultFile, '\n');
    end
    fclose(resultFile);
    
    %Temperature training phase
    formatSpec = 'XFEMTest/Files4Pics/TrainingPhaseTemperature_%d.txt';
    filename = sprintf(formatSpec,modes);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:size(heatFlux, 1)
        fprintf(resultFile, '%d, %\t', x_PostProcess(i));
        for k = 1:trainingTimeSteps*numberOfTimeStepsPerLayer
            fprintf(resultFile, '%d, %\t', temperatureSolution(i,k));
        end
        fprintf(resultFile, '\n');
    end
    fclose(resultFile);
    
    
    %POD modes
    formatSpec = 'XFEMTest/Files4Pics/PODModes_%d.txt';
    filename = sprintf(formatSpec,modes);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:size(solutionReductionOperator, 1)
        fprintf(resultFile, '%d, %\t', i);
        for k = 1:modes
            fprintf(resultFile, '%d, %\t', solutionReductionOperator(i,k));
        end
        fprintf(resultFile, '\n');
    end
    
    fclose(resultFile);
    
    
    
end % modes loop


