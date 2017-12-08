load('results_N_2_ed_s.mat')
infres_RS = infres;
load('results_neutral_N_2_ed_s.mat')
infres_RN = infres;

wc_RS = infres_RS.wc_f;
wc_RN = infres_RN.wc_f;


plot([1,2,3,4], wc_RS(1:4), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red'); hold on;
plot([1 2 3 4], wc_RN(1:4), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'blue')
xlim([0, 5])

%% 
x_rel = 1.:0.1:30;
r_1 = 0.05;
r_2 = 1.;
features_1 =  (x_rel > 2.5).*(log(1+exp(r_1*(x_rel-2.5))) - log(2) ); % not too far from leader
features_2 = (x_rel < 2.5).*( log(1+exp(-r_2*(x_rel-2.5))) -log(2) );

RS_cost = wc_RS(1)*features_1 + wc_RS(2)*features_2;
RN_cost = wc_RN(1)*features_1 + wc_RN(2)*features_2;

figure()
plot(x_rel, RS_cost, 'r', 'LineWidth', 4); hold on;
plot(x_rel, RN_cost, 'b', 'LineWidth', 4);
xlabel('$x_{\textrm{rel}}$', 'Interpreter', 'latex');
%ylabel('Cost w.r.t. combination of $\phi_1$ and $\phi_2$', 'Interpreter', 'latex');
xlim([1., 30])
set(gca,'FontSize',30);
legend('RS model', 'RN model')








