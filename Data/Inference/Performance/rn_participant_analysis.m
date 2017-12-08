%% Load data (risk-averse)
load perf_N_2_p4.mat;

%% Full trajectory with no prediction
tt = 1:780; K = 15;
x_rel_actual = x_expert(6,:) - x_expert(1,:);
figure()
plot(0.1*tt, x_rel_actual(tt)-2.5, 'g', 'LineWidth', 5.);
hold on;
xlabel('time (s)')
ylabel('$x_{\textrm{rel}}$', 'Interpreter', 'latex', 'FontSize', 30)
ylim([0,45])
set(findall(gcf,'type','text'),'FontSize',28);set(gca,'FontSize',28);
%legend('Human', 'Location', 'northwest');

%% Plot x relative over total time
tt = 1:780; K = 15;
x_rel_actual = x_expert(6,:) - x_expert(1,:);
figure()
plot(tt, x_rel_actual(tt), 'r');
hold on;
for k=1:51
    t_k = (k-1)*K+1:k*K+1;
    %plot(t_k, x_rel_1(1,:,k), 'b');
    %hold on;
    plot(t_k, x_rel(1,:,k), 'k');
    hold on;
    plot(t_k, x_rel_n(1,:,k), 'g');
end

%% Sequence
l_car = 4.2;
figure()
ttt = 350:400; ttt = 0.1*ttt;
plot(ttt, (x_rel_actual(350:400)-2.5)/l_car, 'g', 'LineWidth', 5.);
hold on;
for k=25:26
    t_k = (k-1)*K+1:k*K+1;
    t_k = 0.1*t_k;
    plot(t_k, (x_rel(1,:,k)-2.5)/l_car, 'r', 'LineWidth', 5, 'MarkerSize', 3);
    hold on;
    plot(t_k, (x_rel_n(1,:,k)-2.5)/l_car, 'b', 'LineWidth', 5, 'MarkerSize', 3);
    plot(0.1*(K*(k-1)+1)*ones(1,36), 2.5:0.1:6, 'k', 'LineWidth', 4)
    plot(0.1*(K*(k-1)+9)*ones(1,36), 2.5:0.1:6, 'k:', 'LineWidth', 4)
    plot(0.1*(K*k+1)*ones(1,36), 2.5:0.1:6, 'k', 'LineWidth', 4)
    annotation('doublearrow', [0.31 0.41], [0.25,0.25], 'LineWidth', 3);
    text(36.7,2.6,['Prepare'],'HorizontalAlignment','right', 'Fontsize', 25);
    annotation('doublearrow', [0.44 0.52], [0.25 0.25], 'LineWidth', 3);
    text(37.45,2.6,['React'],'HorizontalAlignment','right', 'Fontsize', 25);
    annotation('doublearrow', [0.54 0.65], [0.25 0.25], 'LineWidth', 3);
    text(38.25,2.6,['Prepare'],'HorizontalAlignment','right', 'Fontsize', 25);
    annotation('doublearrow', [0.67 0.75], [0.25 0.25], 'LineWidth', 3);
    text(38.9,2.6,['React'],'HorizontalAlignment','right', 'Fontsize', 25);
end
xlabel('time (s)')
ylabel('$x_{\textrm{rel}}/l_{\textrm{car}}$', 'Interpreter', 'latex')
ylim([2,7])
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30);
legend('Human', 'RS-IRL, T=2', 'RN-IRL, T=2', 'Decision time', 'Disturbance time', 'Location', 'northwest');

%% Polytope
P_simp = Polyhedron('V',eye(4));
proj_idx = [1;2;3];
figure()
P_simp.projection(proj_idx).plot('alpha',0.); hold on;
P_w.projection(proj_idx).plot('color','red', 'alpha', 0.7); hold on;
xlabel('p(1)   (nothing)')
ylabel('p(2)   (acceleration)')
zlabel('p(3)   (deceleration)')
set(gca,'FontSize', 30);
axis equal;
grid off;

%% Polytope 2 dimensions
P_simp = Polyhedron('V',eye(4));
proj_idx = [1;3];
figure()
P_simp.projection(proj_idx).plot('alpha',0.); hold on;
P_w.projection(proj_idx).plot('color','red', 'alpha', 0.7); hold on;
xlabel('p(1)   (nothing)')
%ylabel('p(2)   (acceleration)')
ylabel('p(3)   (deceleration)')
set(gca,'FontSize', 30);
axis equal;
grid off;

%% Absolute bar plots
l_car = 4.2;
fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(1:numExps,err_x_n/l_car,0.7,'facecolor','b','edgecolor','b','facealpha',0.6); hold on
bar(1:numExps,err_x/l_car,0.5,'facecolor','r','edgecolor','r','facealpha',0.45);
xlabel('Demonstration Index in Test Trajectory'); 
legend('RN-IRL, T=2', 'RS-IRL, T=2', 'Location', 'northwest');
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
%saveas(fig,['Plots/Bar_expected/abs_bar_expected_',name_str],'epsc');

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(1:numExps,err_vx_n/l_car,0.7,'facecolor','b','edgecolor','b','facealpha',0.6); hold on
bar(1:numExps,err_vx/l_car,0.5,'facecolor','r','edgecolor','r','facealpha',0.45);
xlabel('Demonstration Index in Test Trajectory'); 
legend('RN-IRL, T=2', 'RS-IRL, T=2', 'Location', 'northwest');
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
%saveas(fig,['Plots/Bar_expected/abs_bar_expected_vx_',name_str],'epsc');

%% Full trajectory plot with boxes
l_car = 4.2;
tt = 1:780; K = 15;

x_rel_actual = x_expert(6,:) - x_expert(1,:);
figure()
plot(0.1*tt, (x_rel_actual(tt)-2.5)/l_car, 'g', 'LineWidth', 5);
hold on;

xlabel('time (s)')
ylabel('$x_{\textrm{rel}}/l_{\textrm{car}}$', 'Interpreter', 'latex')
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30);

%% Percentage plot
fig = figure('units','normalized','outerposition',[0 0 1 1]);
ind_RS = find(err_x_n(1:end-1) > err_x(1:end-1));
ind_RN = find(err_x_n(1:end-1) < err_x(1:end-1));
improvement_RS = (err_x_n(ind_RS)'-err_x(ind_RS)')./err_x_n(ind_RS)';
bar(ind_RS,100*improvement_RS, 'r', 'BarWidth', 0.6);
hold on;
improvement_RN = (err_x_n(ind_RN)'-err_x(ind_RN)')./err_x_n(ind_RN)';
bar(ind_RN,100*improvement_RN, 'b', 'BarWidth', 0.6);
xlabel('Demonstration Index in Test Trajectory');
lim_sup = 100*max(max(improvement_RN), max(improvement_RS));
ylim([-lim_sup-5,lim_sup+5])
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
%saveas(fig,['Plots/Bar_expected/bar_expected_percentage_',name_str],'epsc');


%% Percentage plot vx
fig = figure('units','normalized','outerposition',[0 0 1 1]);
ind_RS = find(err_vx_n(1:end-1) > err_vx(1:end-1));
ind_RN = find(err_vx_n(1:end-1) < err_vx(1:end-1));
improvement_RS = (err_vx_n(ind_RS)'-err_vx(ind_RS)')./err_vx_n(ind_RS)';
bar(ind_RS,100*improvement_RS, 'r', 'BarWidth', 0.6);
hold on;
improvement_RN = (err_vx_n(ind_RN)'-err_vx(ind_RN)')./err_vx_n(ind_RN)';
bar(ind_RN,100*improvement_RN, 'b', 'BarWidth', 0.6);
xlabel('Demonstration Index in Test Trajectory');
lim_sup = 100*max(max(improvement_RN), max(improvement_RS));
ylim([-lim_sup-5,lim_sup+5])
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
%saveas(fig,['Plots/Bar_expected/bar_expected_vx_percentage_',name_str],'epsc');













