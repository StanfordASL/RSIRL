%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Risk-sensitive Inverse Reinforcement Learning %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Choose participant
clear;clc;close all;
name_str = 'p2';
l = 3.476; % length of car
dt = 0.1; % time step discretization
K = 15; % number of time steps within a decision stage
overlap = 8; % overlap between robot's and human's decisions

%%%%%%%%%%%%%%%%%%%%%%%%%%%   Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 2; % number of decision stages

nb_clusters_1 = 15; % number of centroids used in K-Means to discretize the action space in the first decision stage
nb_clusters_2 = 5; % number of centroids used in K-Means to discretize the action space in the second decision stage

beta = 3.0; % inverse temperature of Boltzmann distribution
T = 10; % number of time steps in gradient ascent

k_delay = 3; % offset (tune it manually)

check_offset = false;
load_means = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process training data
if (~check_offset)
    load(['Data/Inference/Means/k_delay_train_',name_str,'.mat'], 'k_delay');
end
process_training_data;

%% Compute mean actions for each decision stage

if (~load_means)
    u_kmeans = {};
    for i=1:numExps
        u_kmeans{i} = u_expert(:,(i-1)*K+1:i*K);
    end
    
    [means_1, ~] = segment_actions(u_kmeans, K, nb_clusters_1);
    action_space_1 = {};
    for i=1:nb_clusters_1
        action_space_1{i} = [means_1(i,1:K); means_1(i,K+1:2*K)];
    end
    
    [means_2, ~] = segment_actions(u_kmeans, K, nb_clusters_2);
    action_space_2 = {};
    for i=1:nb_clusters_2
        action_space_2{i} = [means_2(i,1:K); means_2(i,K+1:2*K)];
    end
    
    u_disc{1} = action_space_1;
    u_disc{2} = action_space_2;
    
    save(['Data/Inference/Means/means_', num2str(N),'_',name_str,'.mat'],'means_1','means_2','action_space_1','action_space_2','u_disc');
else
    load(['Data/Inference/Means/means_', num2str(N),'_',name_str,'.mat']);
end

%{
for k=1:15
    plot(1:20, means_1(k,1:20));
    hold all;
end
%}


%% Load data
data.x = x_expert;
data.counter = counters;
data.dist_last = dist_last;
u = {};
for i=1:numExps
    u{i} = u_expert(:,(i-1)*K+1:i*K);
end
[data.u, data.opt_action] = map_actions(u, means_1);

%% Initialize polytope parameterization and dynamics
parameterized_polytope;
define_dynamics;
par(3) = L;
%% Parameters for gradient ascent
alpha = 0.7;
mask = [1,1];
w_initial = [0.970942973174348
   0.861734604362476
   0.184879392473038
   0.965328976338833
   0.020000000000000
   0.138265395637524
   0.806063580701310
   0.025613996835515
];
%w_initial = 0.15*ones(M,1);
w_initial = tighten(w_initial, A, b, L);
disp('w_initial inside P_proj:');
disp(P_proj.contains(w_initial));

wc_initial = [0.029449850218653
   0.350030083691099
   0.305069896997852
   0.275755130667104
   0.024519450053457
   0.015175588371834];
%wc_initial = [1;5;2;3;1;2]; wc_initial = wc_initial/sum(wc_initial);

%% Outliers

%compute log likelihood to find outliers
l_array = compute_log_like(1:numExps,w_initial, wc_initial, A, b, dynamics, disturbances, data, u_disc, beta, par, LP_sol, LP_prob);
close all;
plot(1:numExps,l_array);

%keyboard;
%Exps = [1,3:28,30:43];
Exps = 1:numExps;

%% Learn polytope and cost weights by maximizing log-likelihood using full-batch gradient ascent

delta_l = 1;
delta_l_eps = 0.01;
w_f = w_initial; wc_f = wc_initial;

while (delta_l > delta_l_eps)
    
    res = subgradientdescent(mask,Exps,w_initial, wc_initial, A, b, dynamics, disturbances, data, u_disc, beta, alpha, T, P_proj, par, LP_sol, LP_prob);
    
    %net change (use for convergence test)
    delta_l = res.log_likelihood(end)-res.log_likelihood(1);
    fprintf('Delta(log_like): %f\n',delta_l);
    disp('**********************************************');
        
    if (delta_l > 0) %successful sweep
        %update record
        l_array = res.log_like_array(:,end);
        l_w = res.log_likelihood;
        w_f = res.w_s(:,end-1); % offsets' parameters
        wc_f = res.wc_s(:,end-1); % cost weights
    
        %update starting point
        w_initial = w_f; wc_initial = wc_f;
    end
end
    
% Recovered polytope
P_simp = Polyhedron('V', eye(L));
P_w = Polyhedron('A', A, 'b', b-w_f) & P_simp;
P_w = P_w.minVRep();

%% Save results
% keyboard;
infres.P_w = P_w;
infres.l_w = l_w;
infres.w_f = w_f;
infres.wc_f = wc_f;

par_ascent.Exps = Exps;
par_ascent.w_initial = w_initial;
par_ascent.wc_initial = wc_initial;
par_ascent.T = T;
par_ascent.beta = beta; 
par_ascent.alpha = alpha;
save(['Data/Inference/results_N_',num2str(N),'_',name_str,'.mat'], 'beta', 'infres', 'disturbances', 'dynamics', 'data', 'x_expert', 'u_disc', 'par', 'par_ascent');









