function res = subgradientdescent(mask,Idx_do,w_0, w_c_0, A, b, dynamics, disturbances, data, u_disc, beta, alpha, T, P_proj, par, LP_sol, LP_prob)
%% parameters
n = par(1); L = par(3); N = par(4); K = par(5);
x_expert = data.x; 
u_expert = data.u;
a_expert = data.opt_action;
dist_last = data.dist_last;
counters = data.counter;
% numExps = length(u_expert);
M = length(w_0);

%% Outliers

numExps = numel(Idx_do);

%% initialize gradient ascent
w = w_0;
w_c = w_c_0;
log_likelihood = zeros(1,T);
log_like_array = zeros(numExps,T);
w_s = zeros(length(w),T);
wc_s = zeros(length(w_c),T);

for t=1:T
    disp(['Time step: ', num2str(t)]);   
    l_array = zeros(numExps,1);
    grad_w = zeros(length(w),numExps);   
    grad_wc = zeros(length(w_c),numExps);
    
    P_w = Polyhedron('A', A, 'b', b-w) & Polyhedron('V', eye(L));
    V_w = P_w.minVRep().V;
    P_w_data = struct('A',A,'b',b-w);
    
    %compute gradient of log-likelihood
    parfor j=1:numExps
        i = Idx_do(j);
        x_obs = x_expert(:,(i-1)*K+1); 
        counter = counters(i); % indicates if robot on left or right lane
        a_obs = a_expert(i,1);

        [~, tau, D_tau_w, ~, D_tau_wc, ~] = est_H(x_obs, 0, M, dynamics, disturbances, dist_last(i),...
                        counter, P_w_data, V_w, u_disc, w_c, beta, 1, par, LP_sol, LP_prob);

        B_weights = exp(-beta*tau)./sum(exp(-beta*tau));

        gradient_W = beta * ( D_tau_w * B_weights - D_tau_w(:,a_obs) );
        gradient_wc = beta * ( D_tau_wc * B_weights - D_tau_wc(:,a_obs) );
        grad_w(:,j) = gradient_W; % gradient of log likelihood of i-th observation wrt w
        grad_wc(:,j) = gradient_wc; % gradient of log likelihood of i-th observation wrt wc
        
        l_array(j) = log( B_weights(a_obs) ); % log likelihood of i-th observation
    end
    
    subgradient_w = sum(grad_w,2)/numExps;
    subgradient_wc = sum(grad_wc,2)/numExps;
    
    disp(subgradient_w)
    disp(subgradient_wc)
    
    log_likelihood(t) = sum(l_array)/ numExps;
    disp(log_likelihood(t))
    log_like_array(:,t) = l_array;

    %Backtracking search:
    step = 1; adv_w = -1; bt_count = 0;
    w_new = w;
    fprintf('adv_w:');
    while (mask(1) && adv_w < 0 && bt_count < 4)
        % update parameter
        w_new = w + (step) * subgradient_w;
        % project parameter on feasible space and tighten polytope (still convex set)
        w_proj = P_proj.project(w_new); w_new = w_proj.x; 
        w_new = tighten(w_new, A, b, L); 
        
        try
            l_array_new = compute_log_like(Idx_do,w_new, w_c, A, b, dynamics, disturbances, data, u_disc, beta, par, LP_sol, LP_prob);
            adv_w = mean(l_array_new) - (log_likelihood(t)+(step/2)*(subgradient_w'*(w_new-w)));
        catch
            adv_w = -1;
        end
        
        if (adv_w < 0)
            step = alpha*step;
        end
        fprintf('%f, ', adv_w);
        bt_count = bt_count + 1;
    end
    if (adv_w < 0)
        w_new = w;
    end
    fprintf('\n');
    
    step = 0.15; adv_wc = -1; bt_count = 0;
    w_c_new = w_c;
    fprintf('adv_wc:');
    while (mask(2) && adv_wc < 0 && bt_count < 4)
        [val,~] =  max(subgradient_wc);
        subgradient_wc = subgradient_wc - (val>0)*val;
        w_c_new = w_c;
        for j = 1:length(w_c)
            w_c_new(j) = w_c(j)*exp((step)*subgradient_wc(j))/...
                (w_c'*exp((step)*subgradient_wc));
        end
        
        l_array_new = compute_log_like(Idx_do,w, w_c_new, A, b, dynamics, disturbances, data, u_disc, beta, par, LP_sol, LP_prob);
        adv_wc = mean(l_array_new) - (log_likelihood(t)+(step/2)*(subgradient_wc'*(w_c_new-w_c)));
        
        if (adv_wc < 0)
            step = alpha*step;
        end
        fprintf('%f, ', adv_wc);
        bt_count = bt_count + 1;
    end
    if (adv_wc < 0)
        w_c_new = wc;
    end
    fprintf('\n'); 
    
    w = w_new;
    disp(w)
    w_s(:,t) = w;
    
    w_c = w_c_new;
    disp(w_c);
    wc_s(:,t) = w_c; 
    
    %early termination
    if (adv_w < 0 && adv_wc < 0)
       for tt = t+1:T
           w_s(:,tt) = w;
           wc_s(:,tt) = w_c;
           log_likelihood(tt) = log_likelihood(t);
           log_like_array(:,tt) = l_array;
       end
       break;
    end
end

res.w_s = w_s;
res.wc_s = wc_s;
res.log_likelihood = log_likelihood;
res.log_like_array = log_like_array;

end



