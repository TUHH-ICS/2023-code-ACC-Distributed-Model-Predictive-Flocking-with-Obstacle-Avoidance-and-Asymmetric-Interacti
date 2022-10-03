%---------------------------------------------------------------------------------------------------
% For Paper
% "Distributed Model Predictive Flocking with Obstacle Avoidance and Asymmetric Interaction Forces"
% by P. Hastedt and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Philipp Hastedt
%---------------------------------------------------------------------------------------------------

% This class is a implementation of the flocking algorithm derived by
% Olfati-Saber in the paper "Flocking for multi-agent dynamic systems: 
% algorithms and theory"
% (https://doi.org/10.1109/TAC.2005.864190)
classdef OlfatiSaberAgent < DoubleIntegratorAgent
    
    % Define constants of the flocking protocol
    properties
        da;             % alpha flocking desired distance
        db;             % beta flocking desired distance
        ra;             % alpha flocking interaction range
        rb;             % beta flocking interaction range
        ha;             % alpha flocking bump function parameter
        hb;             % beta flocking bump function paramter
        epsilon_sigma;  % sigma norm parameter
        da_sigma;       % scaled da
        db_sigma;       % scaled db
        ra_sigma;       % scaled ra
        rb_sigma;       % scaled rb
        c1a;            % P gain alpha flocking
        c1b;            % P gain beta flocking
        c1g;            % P gain gamma flocking 
        c2a;            % D gain alpha flocking
        c2b;            % D gain beta flocking 
        c2g;            % D gain gamma flocking
        pot_a;          % potential field parameter a
        pot_b;          % potential field parameter b
        u_max;          % maximum input
        isSaturated;    % saturation flag
    end
    
    % data to be transmitted in addition to position, velocity, and u
    properties
        num_N;      % number of neighbors
        num_O;      % number of obstacles
        neighbors;  % recognized neighbors
        reference;  % gamma agent
        obstacles;  % set of obstacles
    end
    
    % member variables
    properties(GetAccess = private, SetAccess = private)
        T;              % sampling time
        m;              % space dimension  
        numberAgents;   % number of agents in simulation scenario
    end
    
    methods
        %OlfatiSaberAgent Initialize agent
        % Inputs:
        %   id          : agent id
        %   dT          : sampling time
        %   initialPos  : initial position
        %   initialVel  : initial velocity
        %   cfg         : configuration parameters
        function obj = OlfatiSaberAgent(id, param, initialPos, initialVel, cfg)
            % call parent constructor
            obj@DoubleIntegratorAgent(id, param.dT, initialPos, initialVel);
            
            % load parameters
            obj.da = cfg.d;
            obj.db = cfg.do;
            obj.ra = param.range;
            obj.rb = param.ro;
            obj.ha = cfg.ha;
            obj.hb = cfg.hb;
            obj.epsilon_sigma = cfg.epsilon_sigma;
            
            obj.c1a = cfg.c1a;
            obj.c2a = cfg.c2a;
            obj.c1b = cfg.c1a;
            obj.c2b = cfg.c2a;
            obj.c1g = cfg.c1a;
            obj.c2g = cfg.c2a;
            obj.pot_a = cfg.pot_a;
            obj.pot_b = cfg.pot_b;
            obj.u_max = cfg.u_max;
            
            obj.isSaturated = cfg.isSaturated;
            
            % calculate sigma values
            obj.da_sigma = sigma_norm(obj.da, obj.epsilon_sigma);
            obj.db_sigma = sigma_norm(obj.db, obj.epsilon_sigma);
            obj.ra_sigma = sigma_norm(obj.ra, obj.epsilon_sigma);
            obj.rb_sigma = sigma_norm(obj.rb, obj.epsilon_sigma);
            
            % initialize member variables
            obj.T = param.dT;
            obj.m = cfg.m;
            obj.num_N = 0;
            obj.numberAgents = param.agentCount;
            obj.obstacles = param.obstacles;
            
            % set reference
            if numel(param.reference) == 0
                obj.reference = [];
            elseif numel(param.reference) == 2*obj.m
                obj.reference = param.reference;
            else
                obj.reference = param.reference(:,id);
            end
        end
        
        function step(obj)
            u = zeros(obj.m, 1);
            obj.neighbors = zeros(obj.numberAgents,1);
            % Receive messages from the network
            messages = obj.receive();
            obj.num_N = length(messages);
            
            % alpha flocking component
            if obj.num_N ~= 0
                for message = messages
                    obj.neighbors(message.data.id) = 1;
                    % compute sigma distance and gradient
                    qij = message.data.position - obj.position;
                    pij = message.data.velocity - obj.velocity;
                    [qij_sigma, grad_qij_sigma] = sigma_norm(qij, obj.epsilon_sigma);
                    
                    % calculate u_alpha_q
                    u_alpha_q = obj.c1a * phi_alpha(qij_sigma, obj.ra_sigma, obj.da_sigma, obj.ha, obj.pot_a, obj.pot_b) * grad_qij_sigma;
                    
                    % calculate u_alpha_p
                    aij = rho_h(qij_sigma / obj.ra_sigma, obj.ha);
                    u_alpha_p = obj.c2a * aij * pij;
                    
                    % add to overall input
                    u = u + u_alpha_q + u_alpha_p;
                end
            end
            
            % beta flocking component
            if ~isempty(obj.obstacles)
                [obj.num_O, betaAgents] = getObstaclesInRange([obj.position; obj.velocity], obj.obstacles, obj.rb);
                if obj.num_O ~= 0
                    betaAgents = reshape(betaAgents,2*obj.m,[]);
                    betaAgents(obj.m+1:end,:) = [];
                    for i = 1:obj.num_O
                        % compute sigma distance and gradient
                        qio = betaAgents(:,i) - obj.position;
                        [qio_sigma, grad_qio_sigma] = sigma_norm(qio, obj.epsilon_sigma);
                        
                        % calculate u_beta_q
                        u_beta_q = obj.c1b * phi_beta(qio_sigma, obj.db_sigma, obj.hb) * grad_qio_sigma;
                        
                        % add to overall input
                        u = u + u_beta_q;
                    end
                end
            end
            
            % gamma flocking component
            if ~isempty(obj.reference)
                qir = obj.reference(1:obj.m) - obj.position;
                pir = obj.reference(obj.m+1:end) - obj.velocity;
                [~, grad_qir_sigma] = sigma_norm(qir, 1);
                
                % calculate u_gamma_q
                u_gamma_q = obj.c1g * grad_qir_sigma;
                
                % calculate u_gamma_p
                u_gamma_p = obj.c2g * pir;
                
                % add to overall input
                u = u + u_gamma_q + u_gamma_p;
            end
            
            % saturation
            if obj.isSaturated
                for i = 1:obj.m
                    u(i) = max(-obj.u_max, min(obj.u_max, u(i)));
                end
            end
            
            % Evaluate double integrator dynamics
            obj.move(u);
            
            % Send message to network, include position and velocity
            data = struct;
            data.position = obj.position;
            data.velocity = obj.velocity;
            data.u = u;
            data.num_N = obj.num_N;
            data.id = obj.id;
            data.neighbors = obj.neighbors;
            obj.send(data)
        end
    end
end

