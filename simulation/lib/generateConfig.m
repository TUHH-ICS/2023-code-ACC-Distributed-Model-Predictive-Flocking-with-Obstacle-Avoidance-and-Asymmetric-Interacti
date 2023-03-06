%---------------------------------------------------------------------------------------------------
% For Paper
% "Distributed Model Predictive Flocking with Obstacle Avoidance and Asymmetric Interaction Forces"
% by P. Hastedt and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Philipp Hastedt
%---------------------------------------------------------------------------------------------------

function cfg = generateConfig(algorithms, index)

% Algorithms
% 1: Hastedt MPC
% 2: Huang MPC
% 3: Olfati-Saber

% common parameters
cfg.d       = 7;    % desired inter-agent distance
cfg.do      = 6;    % desired obstacle separation
cfg.m       = 2;    % space dimension
cfg.Hp      = 5;    % prediction horizon
cfg.u_max   = 1;    % input constraint

switch(index)
    % 1: Hastedt MPC parameters
    case (1)
        cfg.lambda_u        = 0.05;     % control weight
        cfg.lambda_beta     = 1.5;      % obstacle weight
        cfg.lambda_a_plus   = 0.1;      % attractive force weight
        cfg.lambda_a_minus  = 1;        % repulsive force weight
        cfg.lookAhead       = 1;        % reference look ahead distance
        cfg.q_pos           = 0.07;     % reference position weight
        cfg.q_vel           = 0.01;     % reference velocity weight
        
    % 2: Huang MPC parameters
    case (2)
        cfg.lambda  = 0.01;         % control weight
        cfg.c       = 0.1;          % velocity consensus weight
        cfg.C1      = 0.2;          % reference position gain
        cfg.C2      = 0.35;         % reference velocity gain
        
    % 3: Olfati-Saber parameters
    case (3)
        cfg.ha = 0.2;               % alpha flocking bump function parameter
        cfg.hb = 0.9;               % beta flocking bump function parameter
        cfg.epsilon_sigma = 0.1;    % sigma norm parameter
        cfg.c1a = 0.9;              % P gain alpha flocking
        cfg.c2a = 2*sqrt(cfg.c1a);  % D gain alpha flocking
        cfg.c1b = 0.9;              % P gain beta flocking
        cfg.c2b = 2*sqrt(cfg.c1b);  % D gain beta flocking
        cfg.c1g = 0.9;              % P gain gamma flocking
        cfg.c2g = 2*sqrt(cfg.c1g);  % D gain gamma flocking
        cfg.pot_a = 4;              % potential field parameter a
        cfg.pot_b = 4;              % potential field parameter b
        cfg.isSaturated = true;     % saturation flag
end
save(strcat(algorithms(index),'/cfg/config.mat'),'cfg');
end

