%---------------------------------------------------------------------------------------------------
% For Paper
% "Distributed Model Predictive Flocking with Obstacle Avoidance and Asymmetric Interaction Forces"
% by P. Hastedt and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Philipp Hastedt
%---------------------------------------------------------------------------------------------------

clear all;
close all;
addpath(genpath('evaluation'),genpath('simulation'), genpath('data'))
dataPath = "data/";

%% Select Data to Evaluate
% videos are available in the "root/video" directory

% Available Data
simData = [
    "hastedt_scenario_1"        % 1
    "hastedt_scenario_2"        % 2
    "huang_scenario_1"          % 3
    "huang_scenario_2"          % 4
    "olfati-saber_scenario_1"   % 5
    "olfati-saber_scenario_2"   % 6
    ];

% Select Data to compare by adding the indices from the simData array
dataSelection = [1,3,5]; % scenario 1 comparison
% dataSelection = [2,4,6]; % scenario 2 comparison

%% Comparison/Evaluation
% Minimum distance variable initialization
t_distN = {};
distN = {};
t_distO = {};
distO = {};

% calculate performance indicators
for i=1:length(dataSelection)
    load(dataPath+simData(dataSelection(i)));
    t_distN{i} = out.t;
    t_distO{i} = out.t;
    [distN{i},distO{i}] = calculateMinimumDistances(out,param,0,7,6);
    distO{i}=changem(distO{i},nan);
end

%% Agent trajectories with markers for initial and final position
for j=1:length(dataSelection)
    figure()
    load(dataPath+simData(dataSelection(j)));
    viscircles(param.obstacles(1:2,:)',param.obstacles(3,:),'Color','black', 'LineWidth', 1); hold on;
    for i = 1:size(out.data.position,3)
        plot(out.data.position(:,1,i),out.data.position(:,2,i),'b'); hold on;
    end
    plot(squeeze(out.data.position(1,1,:)),squeeze(out.data.position(1,2,:)),'kx','MarkerSize',10,'LineWidth',1); hold on;
    plot(squeeze(out.data.position(end,1,:)),squeeze(out.data.position(end,2,:)),'kx','MarkerSize',10,'LineWidth',1); hold on;
    title("Agent Trajectories "+ replace(erase(simData(dataSelection(j)),".mat"),"_","\_"));
    xlabel('x');
    ylabel('y');
    xlim([-5 100]);
    ylim([-5 100]);
end

%% Minimum distances
figure()
for i = 1:size(distN,2)
    plot(t_distN{i},distN{i}, 'DisplayName',replace(erase(simData(dataSelection(i)),".mat"),"_","\_")+" q_{ij}");
    hold on;
    plot(t_distO{i},distO{i}, '--','DisplayName',replace(erase(simData(dataSelection(i)),".mat"),"_","\_")+" q_{io}");
    hold on;
end
xlim([0,400])
ylim([0,8.4])
grid on;
title('Minimum inter-agent and agent-obstacle distances')
xlabel('time in s');
ylabel('distance');
legend show;
legend('Location','southeast');

%% Input Comparison
names = {};
for i = 1:length(dataSelection)
    load(dataPath+simData(dataSelection(i)));
    u_rms(i) = calculateTotalInputUsed(out);
    names{i} = replace(erase(simData(dataSelection(i)),".mat"),"_","\_");
end
figure()
bar(1:i,u_rms);
set(gca, 'XTick', 1:length(names),'XTickLabel',names);
title('RMS Input Values')

%% Minimum Obstacle Distance Comparison
names = {};
for i = 1:length(dataSelection)
    load(dataPath+simData(dataSelection(i)));
    d_min(i) = min(distO{i});
    names{i} = replace(erase(simData(dataSelection(i)),".mat"),"_","\_");
end
figure()
bar(1:i,d_min);
set(gca, 'XTick', 1:length(names),'XTickLabel',names);
title('Minimum Obstacle Distance')

%% Minimum Agent Distance Comparison
names = {};
for i = 1:length(dataSelection)
    load(dataPath+simData(dataSelection(i)));
    d_min(i) = min(distN{i}(101:end));
    names{i} = replace(erase(simData(dataSelection(i)),".mat"),"_","\_");
end
figure()
bar(1:i,d_min);
set(gca, 'XTick', 1:length(names),'XTickLabel',names);
title('Minimum Agent Distance')

%% Obstacle Clearance Time Comparison
names = {};
for i = 1:length(dataSelection)
    load(dataPath+simData(dataSelection(i)));
    indices = find((~isnan(distO{i}))>0.5);
    clearanceTime(i) = t_distO{i}(indices(end))-t_distO{i}(indices(1));
    names{i} = replace(erase(simData(dataSelection(i)),".mat"),"_","\_");
end
figure()
bar(1:i,clearanceTime);
set(gca, 'XTick', 1:length(names),'XTickLabel',names);
title('Obstacle Clearance Time')

%% Potential Comparison
plotPotentials()