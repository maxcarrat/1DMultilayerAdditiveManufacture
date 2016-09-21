%% Generate Convergence Study
% Read data from .txt file and plot the convergnce study of h-FEM and
% XPODFEM w.r.t. Reference overkilled solution.

figure(10000000)
% Create axes
axes1 = axes;

referenceFile = fopen('myMultiPhaseReferenceResultsFile.txt','r');
formatSpec = '%f;\n';
% referenceTemperatureSolution = fscanf(referenceFile, formatSpec);
referenceTemperatureSolution = dlmread('XFEMResults/Temperature/myFEMMultiPhaseResultsFile_8.txt');
referenceFluxSolution = dlmread('XFEMResults/Fluxes/myFEMMultiPhaseFluxesFile_8.txt');


%% hFEM plot
% plot the rlative error on the temperature L2 norm using diffrent level of
% refinement for classical FEM discretization

dofsVector=[41, 81, 161, 321, 641, 1281, 2561, 5121]; % refDepth = 1,..,8
relError = [];
relFluxError = [];
for depth=1:7
    formatSpec = 'XFEMResults/Temperature/myFEMMultiPhaseResultsFile_%d.txt';
    FEMFileID = sprintf(formatSpec,depth);
    FEMFile = fopen(FEMFileID, 'r');
    FEMTemperatureSolution = dlmread(FEMFileID);
    
    formatSpec = 'XFEMResults/Fluxes/myFEMMultiPhaseFluxesFile_%d.txt';
    FEMFileID = sprintf(formatSpec,depth);
    FEMFile = fopen(FEMFileID, 'r');
    FEMFluxSolution = dlmread(FEMFileID);
    
    err = 0.0;
    fluxErr = 0.0;
    for i=2:size(FEMTemperatureSolution,1)-1
        for j=1:size(FEMTemperatureSolution,2)-1
            err = err + sqrt(( referenceTemperatureSolution(i,j) - FEMTemperatureSolution(i,j))^2 ...
                / referenceTemperatureSolution(i,j)^2 ) * 100;
            fluxErr = fluxErr + sqrt(( referenceFluxSolution(i,j) - FEMFluxSolution(i,j))^2 ...
                / referenceFluxSolution(i,j)^2 ) * 100;
        end
    end
    
    relError = [relError err/size(FEMTemperatureSolution,1)/(size(FEMTemperatureSolution,2)-1)];
    relFluxError = [relFluxError fluxErr/size(FEMTemperatureSolution,1)/(size(FEMTemperatureSolution,2)-1)];

end
% Create multiple lines using matrix input to plot
% loglog(dofsVector(1:5),relError,'-kd');
% hold on
loglog(dofsVector(1:7),relFluxError,'-.ko'); 
hold on

% %% X-PODFEM plot M+1
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
%     formatSpec = 'myMultiPhaseXFEMResultsFileShort_IntegrationM+1_%d_ref0.txt';
%     XFEMFileID = sprintf(formatSpec,modes);
%     XFEMFile = fopen(XFEMFileID, 'r');
%     XFEMTemperatureSolution = dlmread(XFEMFileID);
%     
%     formatSpec = 'myMultiPhaseXFEMFluxesFileShort_IntegrationM+1_%d_ref0.txt';
%     XFEMFileID = sprintf(formatSpec,modes);
%     XFEMFile = fopen(XFEMFileID, 'r');
%     XFEMFluxSolution = dlmread(XFEMFileID);
%     
%     err = 0.0;
%     fluxErr = 0.0;
%     for i=2:size(XFEMTemperatureSolution,1)-1
%         
%         timeStepErrorAtPoint = 0.0;
%         timeStepFluxErrorAtPoint = 0.0;
% 
%         for j=1:size(XFEMTemperatureSolution,2)-1
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
% 
% % Create multiple lines using matrix input to plot
% % loglog(dofsVectorPOD(1:4),relErrorPOD, '--b*');
% % hold on
% loglog(dofsVectorPOD(1:3),relFluxErrorPOD, '-.b*');
% hold on
% 
% %% X-PODFEM plot M+5 ref depth 0
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
%     formatSpec = 'myMultiPhaseXFEMResultsFileShort_IntegrationM+5_%d_ref0.txt';
%     XFEMFileID = sprintf(formatSpec,modes);
%     XFEMFile = fopen(XFEMFileID, 'r');
%     XFEMTemperatureSolution = dlmread(XFEMFileID);
%     
%     formatSpec = 'myMultiPhaseXFEMFluxesFileShort_IntegrationM+5_%d_ref0.txt';
%     XFEMFileID = sprintf(formatSpec,modes);
%     XFEMFile = fopen(XFEMFileID, 'r');
%     XFEMFluxSolution = dlmread(XFEMFileID);
%     
%     err = 0.0;
%     fluxErr = 0.0;
%     for i=2:size(XFEMTemperatureSolution,1)-1
%         
%         timeStepErrorAtPoint = 0.0;
%         timeStepFluxErrorAtPoint = 0.0;
% 
%         for j=1:size(XFEMTemperatureSolution,2)-1
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
% % loglog(dofsVectorPOD(1:4),relErrorPOD, '-.md');
% % hold on
% loglog(dofsVectorPOD(1:3),relFluxErrorPOD, '-.md');
% hold on
% 
% % %% X-PODFEM plot M+5 ref depth 1
% % % plot the rlative error on the temperature L2 norm using diffrent level of
% % % refinement for XPODFEM discretization
% % 
% % dofsVectorPOD=[22, 23, 24, 25, 26, 27, 28, 29]; % PODmodes = 1,..,8
% % relErrorPOD = [];
% % relFluxErrorPOD = [];
% % pointwiseRelError = [];
% % pointwiseFluxRelError = [];
% % 
% % for modes=1:3
% %     formatSpec = 'myXFEMMultiPhaseResultsFile_IntegrationM+5_%d_ref1.txt';
% %     XFEMFileID = sprintf(formatSpec,modes);
% %     XFEMFile = fopen(XFEMFileID, 'r');
% %     XFEMTemperatureSolution = dlmread(XFEMFileID);
% %     
% %     formatSpec = 'myXFEMMultiPhaseFluxesFile_IntegrationM+5_%d_ref1.txt';
% %     XFEMFileID = sprintf(formatSpec,modes);
% %     XFEMFile = fopen(XFEMFileID, 'r');
% %     XFEMFluxSolution = dlmread(XFEMFileID);
% %     
% %     err = 0.0;
% %     fluxErr = 0.0;
% %     for i=1320:size(XFEMTemperatureSolution,1)-1
% %         
% %         timeStepErrorAtPoint = 0.0;
% %         timeStepFluxErrorAtPoint = 0.0;
% % 
% %         for j=1:size(XFEMTemperatureSolution,2)-1
% %             pointError = sqrt(( referenceTemperatureSolution(i,j) - XFEMTemperatureSolution(i,j))^2 ...
% %                 / referenceTemperatureSolution(i,j)^2 );
% %             timeStepErrorAtPoint = timeStepErrorAtPoint + pointError;
% %             pointFluxError = sqrt(( referenceFluxSolution(i,j) - XFEMFluxSolution(i,j))^2 ...
% %                 / referenceFluxSolution(i,j)^2 );            
% %             timeStepFluxErrorAtPoint = timeStepFluxErrorAtPoint + pointFluxError;
% % 
% %         end
% %         pointwiseRelError = [pointwiseRelError timeStepErrorAtPoint/(size(XFEMTemperatureSolution,2)-1)];
% %         pointwiseFluxRelError = [pointwiseFluxRelError timeStepFluxErrorAtPoint/(size(XFEMTemperatureSolution,2)-1)];
% % 
% %         err = err + timeStepErrorAtPoint * 100;
% %         fluxErr = fluxErr + timeStepFluxErrorAtPoint * 100;
% %     end
% %     
% %     relErrorPOD = [relErrorPOD err/(size(XFEMTemperatureSolution, 1)-1)/(size(XFEMTemperatureSolution,2)-1)];
% %     relFluxErrorPOD = [relFluxErrorPOD fluxErr/(size(XFEMTemperatureSolution, 1)-1)/(size(XFEMTemperatureSolution,2)-1)];
% % end
% 
% % % Create multiple lines using matrix input to plot
% % % loglog(dofsVectorPOD(1:4),relErrorPOD, '-.bo');
% % % hold on
% % % loglog(dofsVectorPOD(1:1),relFluxErrorPOD, '-.bo');
% % % hold on
% % 
% % % X-PODFEM plot M+5 ref depth 2
% % % plot the rlative error on the temperature L2 norm using diffrent level of
% % % refinement for XPODFEM discretization
% % 
% % dofsVectorPOD=[22, 23, 24, 25, 26, 27, 28, 29]; % PODmodes = 1,..,8
% % relErrorPOD = [];
% % relFluxErrorPOD = [];
% % pointwiseRelError = [];
% % pointwiseFluxRelError = [];
% % 
% % for modes=1:4
% %     formatSpec = 'myFEMMultiPhaseResultsFile_p%d.txt';
% %     XFEMFileID = sprintf(formatSpec,modes);
% %     XFEMFile = fopen(XFEMFileID, 'r');
% %     XFEMTemperatureSolution = dlmread(XFEMFileID);
% %     
% %     formatSpec = 'myFEMMultiPhaseFluxesFile_p%d.txt';
% %     XFEMFileID = sprintf(formatSpec,modes);
% %     XFEMFile = fopen(XFEMFileID, 'r');
% %     XFEMFluxSolution = dlmread(XFEMFileID);
% %     
% %     err = 0.0;
% %     fluxErr = 0.0;
% %     for i=1320:size(XFEMTemperatureSolution,1)-1
% %         
% %         timeStepErrorAtPoint = 0.0;
% %         timeStepFluxErrorAtPoint = 0.0;
% % 
% %         for j=1:size(XFEMTemperatureSolution,2)-1
% %             pointError = sqrt(( referenceTemperatureSolution(i,j) - XFEMTemperatureSolution(i,j))^2 ...
% %                 / referenceTemperatureSolution(i,j)^2 );
% %             timeStepErrorAtPoint = timeStepErrorAtPoint + pointError;
% %             pointFluxError = sqrt(( referenceFluxSolution(i,j) - XFEMFluxSolution(i,j))^2 ...
% %                 / referenceFluxSolution(i,j)^2 );            
% %             timeStepFluxErrorAtPoint = timeStepFluxErrorAtPoint + pointFluxError;
% % 
% %         end
% %         pointwiseRelError = [pointwiseRelError timeStepErrorAtPoint/(size(XFEMTemperatureSolution,2)-1)];
% %         pointwiseFluxRelError = [pointwiseFluxRelError timeStepFluxErrorAtPoint/(size(XFEMTemperatureSolution,2)-1)];
% % 
% %         err = err + timeStepErrorAtPoint * 100;
% %         fluxErr = fluxErr + timeStepFluxErrorAtPoint * 100;
% %     end
% %     
% %     relErrorPOD = [relErrorPOD err/(size(XFEMTemperatureSolution, 1)-1)/(size(XFEMTemperatureSolution,2)-1)];
% %     relFluxErrorPOD = [relFluxErrorPOD fluxErr/(size(XFEMTemperatureSolution, 1)-1)/(size(XFEMTemperatureSolution,2)-1)];
% % end
% % 
% % % Create multiple lines using matrix input to plot
% % loglog(dofsVectorPOD(1:4),relErrorPOD, '-.bo');
% % hold on
% % loglog(dofsVectorPOD(1:1),relFluxErrorPOD, '-.r*');
% % hold on

%% X-PODFEM plot M+10
% plot the rlative error on the temperature L2 norm using diffrent level of
% refinement for XPODFEM discretization

dofsVectorPOD=[22, 23, 24, 25, 26, 27, 28, 29]; % PODmodes = 1,..,8
relErrorPOD = [];
relFluxErrorPOD = [];
pointwiseRelError = [];
pointwiseFluxRelError = [];

for modes=1:3
    formatSpec = 'XFEMResults/Temperature/myXFEMMultiPhaseResultsFile_%d.txt';
    XFEMFileID = sprintf(formatSpec,modes);
    XFEMFile = fopen(XFEMFileID, 'r');
    XFEMTemperatureSolution = dlmread(XFEMFileID);
    
    formatSpec = 'XFEMResults/Fluxes/myXFEMMultiPhaseFluxesFile_%d.txt';
    XFEMFileID = sprintf(formatSpec,modes);
    XFEMFile = fopen(XFEMFileID, 'r');
    XFEMFluxSolution = dlmread(XFEMFileID);
    
    err = 0.0;
    fluxErr = 0.0;
    for i=2:size(XFEMTemperatureSolution,1)-1
        
        timeStepErrorAtPoint = 0.0;
        timeStepFluxErrorAtPoint = 0.0;

        for j= 1:size(XFEMTemperatureSolution,2)-1
            pointError = sqrt(( referenceTemperatureSolution(i,j) - XFEMTemperatureSolution(i,j))^2 ...
                / referenceTemperatureSolution(i,j)^2 );
            timeStepErrorAtPoint = timeStepErrorAtPoint + pointError;
            pointFluxError = sqrt(( referenceFluxSolution(i,j) - XFEMFluxSolution(i,j))^2 ...
                / referenceFluxSolution(i,j)^2 );            
            timeStepFluxErrorAtPoint = timeStepFluxErrorAtPoint + pointFluxError;

        end
        pointwiseRelError = [pointwiseRelError timeStepErrorAtPoint/(size(XFEMTemperatureSolution,2)-1)];
        pointwiseFluxRelError = [pointwiseFluxRelError timeStepFluxErrorAtPoint/(size(XFEMTemperatureSolution,2)-1)];

        err = err + timeStepErrorAtPoint * 100;
        fluxErr = fluxErr + timeStepFluxErrorAtPoint * 100;
    end
    
    relErrorPOD = [relErrorPOD err/(size(XFEMTemperatureSolution, 1)-1)/(size(XFEMTemperatureSolution,2)-1)];
    relFluxErrorPOD = [relFluxErrorPOD fluxErr/(size(XFEMTemperatureSolution, 1)-1)/(size(XFEMTemperatureSolution,2)-1)];
end

% Create multiple lines using matrix input to plot
% loglog(dofsVectorPOD(1:4),relErrorPOD, '--b*');
% hold on
loglog(dofsVectorPOD(1:3),relFluxErrorPOD, '-.r+');
hold on

% %% X-PODFEM plot M+20
% % plot the rlative error on the temperature L2 norm using diffrent level of
% % refinement for XPODFEM discretization
% 
% dofsVectorPOD=[22, 23, 24, 25, 26, 27, 28, 29]; % PODmodes = 1,..,8
% relErrorPOD = [];
% relFluxErrorPOD = [];
% pointwiseRelError = [];
% pointwiseFluxRelError = [];
% 
% for modes=1:4
%     formatSpec = 'myXFEMMultiPhaseResultsFile_IntegrationM+20_%d.txt';
%     XFEMFileID = sprintf(formatSpec,modes);
%     XFEMFile = fopen(XFEMFileID, 'r');
%     XFEMTemperatureSolution = dlmread(XFEMFileID);
%     
%     formatSpec = 'myXFEMMultiPhaseFluxesFile_IntegrationM+20_%d.txt';
%     XFEMFileID = sprintf(formatSpec,modes);
%     XFEMFile = fopen(XFEMFileID, 'r');
%     XFEMFluxSolution = dlmread(XFEMFileID);
%     
%     err = 0.0;
%     fluxErr = 0.0;
%     for i=1:size(XFEMTemperatureSolution,1)-1
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
% loglog(dofsVectorPOD(1:4),relFluxErrorPOD, '-.co');
% hold on
% 
% Title
title({'\centering{\quad \quad \quad Non-linear X-PODFEM}'; '\centering{using refined hFEM (depth 8) as reference}'}...
    , 'Interpreter','latex');

% Legend
legend('FEM', 'POD (M+5) 1 Xelement', 'POD (M+5) 2 Xelements', 'POD (M+5) 4 Xelements','Location', 'southeast');

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
% legend('FEM', 'POD', 'Location', 'southeast');

% Axis labels
ylabel('\fontname{Latin Modern Math} Relative error on average temperature (last layer)');
xlabel('\fontname{Latin Modern Math} bar length');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');



% figure(99999)
% formatSpec = 'myXFEMMultiPhaseTimeFile_IntegrationM+5_1_ref0.txt';
% XFEMTimeFileID = sprintf(formatSpec,modes);
% XFEMTime = dlmread(XFEMTimeFileID);
% t = linspace(0, numel(XFEMTime), numel(XFEMTime)-1);
% plot(t, XFEMTime(1:end-1),'-bo')
% hold on
% 
% formatSpec = 'myXFEMMultiPhaseTimeFile_IntegrationM+5_1_ref1.txt';
% XFEMTimeFileID = sprintf(formatSpec,modes);
% XFEMTime = dlmread(XFEMTimeFileID);
% t = linspace(0, numel(XFEMTime), numel(XFEMTime)-1);
% plot(t, XFEMTime(1:end-1),'-rx')
% hold on
% 
% formatSpec = 'myXFEMMultiPhaseTimeFile_IntegrationM+5_1_ref2.txt';
% XFEMTimeFileID = sprintf(formatSpec,modes);
% XFEMTime = dlmread(XFEMTimeFileID);
% t = linspace(0, numel(XFEMTime), numel(XFEMTime)-1);
% plot(t, XFEMTime(1:end-1),'-m*')
% hold on
% 
% formatSpec = 'myFEMMultiPhaseTimeFile_6.txt';
% FEMTimeFileID = sprintf(formatSpec,modes);
% FEMTime = dlmread(FEMTimeFileID);
% t = linspace(0, numel(FEMTime), numel(FEMTime)-1);
% plot(t, FEMTime(1:end-1),'-kd')
% hold on
% 
% % Title
% title({'CPU Time at each time step'}...
%     , 'Interpreter','latex');
% 
% % Legend
% legend('POD-1 Mode, 1 Xelement', 'POD-1 Modes, 2 Xelements', 'POD-1 Modes, 4 Xelements', 'FEM 64 elements', 'Location', 'southeast');
% 
% % Axis labels
% ylabel('\fontname{Latin Modern Math} CPU Time [msec]');
% xlabel('\fontname{Latin Modern Math} Time Step');
% 
% box(axes1,'on');
% % Set the remaining axes properties
% set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
% hold off

