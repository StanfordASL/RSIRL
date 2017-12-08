function res = subgradientdescent_neutral(Idx_do,w_c_0, p_true, dynamics, disturbances, data, u_disc, beta, alpha, T, par)
%% parameters
n = par(1); L = par(3); N = par(4); K = par(5);
x_expert = data.x; 
u_expert = data.u;
a_expert = data.opt_action;
counters = data.counter;
dist_last = data.dist_last;
% numExps = length(u_expert);

%% Outliers

numExps = numel(Idx_do);

%% initialize gradient ascent
w_c = w_c_0;
log_likelihood = zeros(1,T);
log_like_array = zeros(numExps,T);
wc_s = zeros(length(w_c),T);

for t=1:T
    disp(['Time step: ', num2str(t)]);   
    l_array = zeros(numExps,1); 
    grad_wc = zeros(length(w_c),numExps);
    
    %compute gradient of log-likelihood
    parfor j=1:numExps
        i = Idx_do(j);
        x_obs = x_expert(:,(i-1)*K+1); 
        counter = counters(i); % indicates if robot on left or right lane
        a_obs = a_expert(i,1);
        
        %disp(['iteration: ', num2str(i)])
        %disp(['counter: ', num2str(counter)])
        %disp(['last disturbance: ', num2str(dist_last(i))])
        
        [~, tau, D_tau_wc, ~] = est_H_neutral(x_obs, 0, dynamics, disturbances, dist_last(i),...
                        counter, p_true, u_disc, w_c, beta, 1, par);

        B_weights = exp(-beta*tau)./sum(exp(-beta*tau));

        gradient_wc = beta * ( D_tau_wc * B_weights - D_tau_wc(:,a_obs) );
        grad_wc(:,j) = gradient_wc; % gradient of log likelihood of i-th observation wrt wc
        
        l_array(j) = log( B_weights(a_obs) ); % log likelihood of i-th observation
    end
    
    subgradient_wc = sum(grad_wc,2)/numExps;
    
    disp(subgradient_wc)
    
    log_likelihood(t) = sum(l_array)/ numExps;
    disp(log_likelihood(t))
    
    log_like_array(:,t) = l_array;

    step = 0.15; adv = -1; bt_count = 0;
    w_c_new = w_c;
    fprintf('adv_wc:');
    while (adv < 0 && bt_count < 4)
        [val,~] =  max(subgradient_wc);
        subgradient_wc = subgradient_wc - (val>0)*val;
        w_c_new = w_c;
        
        for j = 1:length(w_c)
            w_c_new(j) = w_c(j)*exp((step)*subgradient_wc(j))/...
                (w_c'*exp((step)*subgradient_wc));
        end

        l_array_new = compute_log_like_neutral(Idx_do,p_true,w_c_new,dynamics, disturbances, data, u_disc, beta, par);
        adv = mean(l_array_new) - (log_likelihood(t)+(step/2)*(subgradient_wc'*(w_c_new-w_c)));
        
        if (adv < 0)
            step = alpha*step;
        end
        fprintf('%f, ', adv);
        bt_count = bt_count + 1;
    end
    if (adv < 0)
        w_c_new = w_c;
    end
    fprintf('\n'); 
    
    w_c = w_c_new;
    disp(w_c);
    wc_s(:,t) = w_c;
    
    %early termination
    if (adv < 0)
       for tt = t+1:T
           wc_s(:,tt) = w_c;
           log_likelihood(tt) = log_likelihood(t);
           log_like_array(:,tt) = l_array;
       end
       break;
    end
    
end

res.wc_s = wc_s;
res.log_likelihood = log_likelihood;
res.log_like_array = log_like_array;

end





