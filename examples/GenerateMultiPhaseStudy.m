%% Generate Convergence Study
% Read data from .txt file and plot the convergnce study of h-IGA and
% XPODIGA w.r.t. Reference overkilled solution.

figure(10000000)
% Create axes
axes1 = axes;

referenceFile = fopen('myMultiPhaseReferenceResultsFile.txt','r');
formatSpec = '%f;\n';
% referenceTemperatureSolution = fscanf(referenceFile, formatSpec);
% referenceTemperatureSolution = dlmread('XIGAResults/Temperature/myIGAMultiPhaseResultsFile_8.txt');
referenceTemperatureSolution = dlmread('PODXIGA/Temperature/myPODXIGAResultsFileComparison_6.txt');
referenceFluxSolution = dlmread('PODXIGA/Fluxes/myPODXIGAFluxesFileComparison_6.txt');
% referenceFluxSolution = dlmread('XIGAResults/Fluxes/myIGAMultiPhaseFluxesFile_8.txt');
% referenceFluxSolution = dlmread('XFEMResults/Fluxes/myFEMNonLinearFluxesFile_8.txt');


%% hIGA plot
% plot the rlative error on the temperature L2 norm using diffrent level of
% refinement for classical IGA discretization

dofsVector=[22, 23, 25, 29, 37, 53, 85, 149, 277]; % refDepth = 1,..,8
relError = [];
relFluxError = [];

for depth=0:4
    formatSpec = 'PODXIGA/Temperature/myPODXIGAResultsFileComparison_%d.txt';
    IGAFileID = sprintf(formatSpec,depth);
    IGATemperatureSolution = dlmread(IGAFileID);
    
    formatSpec = 'PODXIGA/Fluxes/myPODXIGAFluxesFileComparison_%d.txt';
    IGAFileID = sprintf(formatSpec,depth);
    IGAFile = fopen(IGAFileID, 'r');
    IGAFluxSolution = dlmread(IGAFileID);
    
    err = 0.0;
    fluxErr = 0.0;
    
    for i=2:size(IGATemperatureSolution,1)-1
        
        for j=2:size(IGATemperatureSolution,2)-1
            
            err = err + sqrt(( referenceTemperatureSolution(i,j) - IGATemperatureSolution(i,j))^2 ...
                / referenceTemperatureSolution(i,j)^2 ) * 100;
            fluxErr = fluxErr + sqrt(( referenceFluxSolution(i,j) - IGAFluxSolution(i,j))^2 ...
                / referenceFluxSolution(i,j)^2 ) * 100;
        end
    end
    
    relError = [relError err/(size(IGATemperatureSolution, 1)-2)/(size(IGATemperatureSolution,2)-2)];
    relFluxError = [relFluxError fluxErr/(size(IGATemperatureSolution, 1)-2)/(size(IGATemperatureSolution,2)-2)];

end
% Create multiple lines using matrix input to plot
loglog(dofsVector(1:5),relError,'-kd');
hold on
% loglog(dofsVector(1:6),relFluxError,'-.ko'); 
% hold on
grid on


%% X-PODIGA plot 
% plot the rlative error on the temperature L2 norm using diffrent level of
% refinement for XPODIGA discretization

dofsVectorPOD=[23, 24, 25, 26, 27, 28, 29, 30]; % PODmodes = 1,..,8
relErrorPOD = [];
relFluxErrorPOD = [];
pointwiseRelError = [];
pointwiseFluxRelError = [];

for modes=1:3
    formatSpec = 'PODXIGA/Temperature/myeXtendedPODXIGAResultsFileComparison_%d.txt';
    XIGAFileID = sprintf(formatSpec,modes);
    XIGAFile = fopen(XIGAFileID, 'r');
    XIGATemperatureSolution = dlmread(XIGAFileID);
    
    formatSpec = 'PODXIGA/Fluxes/myeXtendedPODXIGAFluxesFileComparison_%d.txt';
    XIGAFileID = sprintf(formatSpec,modes);
    XIGAFile = fopen(XIGAFileID, 'r');
    XIGAFluxSolution = dlmread(XIGAFileID);
    
    err = 0.0;
    fluxErr = 0.0;
    for i=2:size(XIGATemperatureSolution,1)-10
        
        timeStepErrorAtPoint = 0.0;
        timeStepFluxErrorAtPoint = 0.0;

        for j= 2:size(XIGATemperatureSolution,2)-1
            pointError = sqrt(( referenceTemperatureSolution(i,j) - XIGATemperatureSolution(i,j))^2 ...
                / referenceTemperatureSolution(i,j)^2 );
            
            timeStepErrorAtPoint = timeStepErrorAtPoint + pointError;
            
            pointFluxError = sqrt(( referenceFluxSolution(i,j) - XIGAFluxSolution(i,j))^2 ...
                / referenceFluxSolution(i,j)^2 );            
            timeStepFluxErrorAtPoint = timeStepFluxErrorAtPoint + pointFluxError;

        end
        pointwiseRelError = [pointwiseRelError timeStepErrorAtPoint/(size(XIGATemperatureSolution,2)-2)];
        pointwiseFluxRelError = [pointwiseFluxRelError timeStepFluxErrorAtPoint/(size(XIGATemperatureSolution,2)-2)];

        err = err + timeStepErrorAtPoint * 100;
        fluxErr = fluxErr + timeStepFluxErrorAtPoint * 100;
    end
    
    relErrorPOD = [relErrorPOD err/(size(XIGATemperatureSolution, 1)-11)/(size(XIGATemperatureSolution,2)-2)];
    relFluxErrorPOD = [relFluxErrorPOD fluxErr/(size(XIGATemperatureSolution, 1)-11)/(size(XIGATemperatureSolution,2)-2)];
end

% Create multiple lines using matrix input to plot
loglog(dofsVectorPOD(1:3),relErrorPOD, '--b*');
hold on
% loglog(dofsVectorPOD(1:3),relFluxErrorPOD, '-.r+');
% hold on
grid on

% %% X-PODFEM plot
% % plot the rlative error on the temperature L2 norm using diffrent level of
% % refinement for XPODFEM discretization
% 
% dofsVectorPOD=[22, 23, 24, 25, 26, 27, 28, 29]; % PODmodes = 1,..,8
% relErrorPOD = [];
% relFluxErrorPOD = [];
% pointwiseRelError = [];
% pointwiseFluxRelError = [];
% 
% for modes=1:3
%     formatSpec = 'XFEMResults/Temperature/myXFEMMultiPhaseResultsFile_%d.txt';
%     XFEMFileID = sprintf(formatSpec,modes);
%     XFEMFile = fopen(XFEMFileID, 'r');
%     XFEMTemperatureSolution = dlmread(XFEMFileID);
%     
%     formatSpec = 'XFEMResults/Fluxes/myXFEMMultiPhaseFluxesFile_%d.txt';
%     XFEMFileID = sprintf(formatSpec,modes);
%     XFEMFile = fopen(XFEMFileID, 'r');
%     XFEMFluxSolution = dlmread(XFEMFileID);
%     
%     err = 0.0;
%     fluxErr = 0.0;
%     for i=5750:size(XFEMTemperatureSolution,1)-1
%         
%         timeStepErrorAtPoint = 0.0;
%         timeStepFluxErrorAtPoint = 0.0;
% 
%         for j= 1:size(XFEMTemperatureSolution,2)-1
%             pointError = sqrt(( referenceTemperatureSolution(i,j) - XFEMTemperatureSolution(i,j))^2 ...
%                 / referenceTemperatureSolution(i,j)^2 );
%             timeStepErrorAtPoint = timeStepErrorAtPoint + pointError;
%             pointFluxError = sqrt(( referenceFluxSolution(i,j) - XFEMFluxSolution(i,j))^2 ...
%                 / referenceFluxSolution(i,j)^2 );            
%             timeStepFluxErrorAtPoint = timeStepFluxErrorAtPoint + pointFluxError;
% 
%         end
%         pointwiseRelError = [pointwiseRelError timeStepErrorAtPoint/(size(XFEMTemperatureSolution,2)-1)];
%         pointwiseFluxRelError = [pointwiseFluxRelError timeStepFluxErrorAtPoint/(size(XFEMTemperatureSolution,2)-1)];
% 
%         err = err + timeStepErrorAtPoint * 100;
%         fluxErr = fluxErr + timeStepFluxErrorAtPoint * 100;
%     end
%     
%     relErrorPOD = [relErrorPOD err/(size(XFEMTemperatureSolution, 1)-1)/(size(XFEMTemperatureSolution,2)-1)];
%     relFluxErrorPOD = [relFluxErrorPOD fluxErr/(size(XFEMTemperatureSolution, 1)-1)/(size(XFEMTemperatureSolution,2)-1)];
% end
% 
% % Create multiple lines using matrix input to plot
% % loglog(dofsVectorPOD(1:4),relErrorPOD, '--b*');
% % hold on
% loglog(dofsVectorPOD(1:3),relFluxErrorPOD, '-.r+');
% hold on
% grid on

% Title
title({'\centering{\quad \quad \quad Non-linear X-PODIGA}'; '\centering{using refined hIGA (depth 8) as reference}'}...
    , 'Interpreter','latex');

% Legend
legend('IGA depth=1..7', 'X-PODIGA m=1..5','Location', 'southeast');

% Axis labels
ylabel('\fontname{Latin Modern Math} Relative error on the energy norm [%]');
xlabel('\fontname{Latin Modern Math} DOFs');

box(axes1,'on');
grid on;
% Set the remaining axes properties
set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');


hold off

figure(123456789)
plot(linspace(0.0, 0.001, numel(pointwiseRelError)), pointwiseRelError,'-r');
hold on
plot(linspace(0.0, 0.001, numel(pointwiseFluxRelError)), pointwiseFluxRelError, '-.b');
hold on
% Title
title({'Pointwise relative error'}...
    , 'Interpreter','latex');

% % Legend
% legend('IGA', 'POD', 'Location', 'southeast');

% Axis labels
ylabel('\fontname{Latin Modern Math} Relative error on average temperature (last layer)');
xlabel('\fontname{Latin Modern Math} bar length');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');



figure(99999)
formatSpec = 'XIGAResults/Time/myXIGAMultiPhaseTimeFile_1.txt';
XIGATimeFileID = sprintf(formatSpec,modes);
XIGATime = dlmread(XIGATimeFileID);
t = linspace(0, numel(XIGATime), numel(XIGATime)-1);
semilogy(t, XIGATime(1:end-1),'-bx')
hold on

formatSpec = 'XIGAResults/Time/myXIGAMultiPhaseTimeFile_2.txt';
XIGATimeFileID = sprintf(formatSpec,modes);
XIGATime = dlmread(XIGATimeFileID);
t = linspace(0, numel(XIGATime), numel(XIGATime)-1);
semilogy(t, XIGATime(1:end-1),'-rx')
hold on

formatSpec = 'XIGAResults/Time/myXIGAMultiPhaseTimeFile_3.txt';
XIGATimeFileID = sprintf(formatSpec,modes);
XIGATime = dlmread(XIGATimeFileID);
t = linspace(0, numel(XIGATime), numel(XIGATime)-1);
semilogy(t, XIGATime(1:end-1),'-mx')
hold on

formatSpec = 'XIGAResults/Time/myIGAMultiPhaseTimeFile_5.txt';
IGATimeFileID = sprintf(formatSpec,modes);
IGATime = dlmread(IGATimeFileID);
t = linspace(0, numel(IGATime), numel(IGATime)-1);
semilogy(t, IGATime(1:end-1),'-kd')
hold on

formatSpec = 'XIGAResults/Time/myIGAMultiPhaseTimeFile_6.txt';
IGATimeFileID = sprintf(formatSpec,modes);
IGATime = dlmread(IGATimeFileID);
t = linspace(0, numel(IGATime), numel(IGATime)-1);
semilogy(t, IGATime(1:end-1),'-yd')
hold on

formatSpec = 'XIGAResults/Time/myIGAMultiPhaseTimeFile_7.txt';
IGATimeFileID = sprintf(formatSpec,modes);
IGATime = dlmread(IGATimeFileID);
t = linspace(0, numel(IGATime), numel(IGATime)-1);
semilogy(t, IGATime(1:end-1),'-gd')
hold on

formatSpec = 'XIGAResults/Time/myIGAMultiPhaseTimeFile_8.txt';
IGATimeFileID = sprintf(formatSpec,modes);
IGATime = dlmread(IGATimeFileID);
t = linspace(0, numel(IGATime), numel(IGATime)-1);
semilogy(t, IGATime(1:end-1),'-rd')
hold on

% Title
title({'CPU Time at each time step'}...
    , 'Interpreter','latex');

% Legend
legend('X-PODIGA m = 1', 'X-PODIGA m = 2', 'X-PODIGA m = 3', 'IGA depth 5', 'IGA depth 6', 'IGA depth 7', 'Reference', 'Location', 'northeast');

% Axis labels
ylabel('\fontname{Latin Modern Math} CPU Time [msec]');
xlabel('\fontname{Latin Modern Math} Time Step');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
hold off

