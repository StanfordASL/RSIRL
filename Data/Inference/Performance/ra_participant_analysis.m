%% Load data (risk-averse)
load perf_N_2_p2.mat;

%% Full trajectory with no prediction
tt = 1:780; K = 15;
x_rel_actual = x_expert(6,:) - x_expert(1,:);
figure()
plot(0.1*tt, x_rel_actual(tt), 'g', 'LineWidth', 5.);
hold on;
xlabel('time (s)')
ylabel('$x_{\textrm{rel}}$', 'Interpreter', 'latex', 'FontSize', 30)
ylim([0,45])
set(findall(gcf,'type','text'),'FontSize',28);set(gca,'FontSize',28);
%legend('Human', 'Location', 'northwest');

%% Full trajectory with prediction
tt = 1:780; K = 15;
x_rel_actual = x_expert(6,:) - x_expert(1,:);
figure()
plot(tt, x_rel_actual(tt), 'g');
hold on;
for k=1:51
    t_k = (k-1)*K+1:k*K+1;
    %plot(t_k, x_rel_1(1,:,k), 'b');
    %hold on;
    plot(t_k, x_rel(1,:,k), 'b');
    hold on;
    plot(t_k, x_rel_n(1,:,k), 'r');
end

%% Sequence 1
%brown: [0.54 0.27 0.07]
l_car = 4.2; x_rel_actual = x_expert(6,:) - x_expert(1,:);
figure()
ttt = 710:760; ttt = 0.1*ttt;
plot(ttt, (x_rel_actual(710:760))/l_car, 'g', 'LineWidth', 5.);
hold on;
for k=49:50
    t_k = (k-1)*K+1:k*K+1;
    t_k = 0.1*t_k;
    plot(t_k, (x_rel(1,:,k))/l_car, 'r', 'LineWidth', 5, 'MarkerSize', 3);
    hold on;
    plot(t_k, (x_rel_n(1,:,k))/l_car, 'b', 'LineWidth', 5, 'MarkerSize', 3);
    plot(0.1*(K*(k-1)+1)*ones(1,41), 2.:0.1:6, 'k', 'LineWidth', 4)
    plot(0.1*(K*(k-1)+9)*ones(1,41), 2.:0.1:6, 'k:', 'LineWidth', 4)
    plot(0.1*(K*k+1)*ones(1,41), 2.:0.1:6, 'k', 'LineWidth', 4)
    annotation('doublearrow', [0.26 0.39], [0.6,0.6], 'LineWidth', 3);
    text(72.75,5.15,['Prepare'],'HorizontalAlignment','right', 'Fontsize', 25);
    annotation('doublearrow', [0.42 0.52], [0.6 0.6], 'LineWidth', 3);
    text(73.4,5.15,['React'],'HorizontalAlignment','right', 'Fontsize', 25);
    annotation('doublearrow', [0.55 0.68], [0.6 0.6], 'LineWidth', 3);
    text(74.25,5.15,['Prepare'],'HorizontalAlignment','right', 'Fontsize', 25);
    annotation('doublearrow', [0.71 0.81], [0.6 0.6], 'LineWidth', 3);
    text(74.9,5.15,['React'],'HorizontalAlignment','right', 'Fontsize', 25);
end
text(73.7,0.9,['Collision'],'HorizontalAlignment','right', 'Fontsize', 28, 'Color', [0.6 0.6 0.6]);
plot(ttt, 0.59*ones(1,length(ttt)), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 5.);
xlabel('time (s)')
ylabel('$x_{\textrm{rel}}/l_{\textrm{car}}$', 'Interpreter', 'latex')
xlim([71.5, 75.5])
ylim([0.,8])
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30);
legend('Human', 'RS-IRL, T=2', 'RN-IRL, T=2', 'Decision time', 'Disturbance time', 'Location', 'northwest');

%% Sequence 2
l_car = 4.2; x_rel_actual = x_expert(6,:) - x_expert(1,:);
figure()
ttt = 300:405; ttt = 0.1*ttt;
plot(ttt, (x_rel_actual(300:405))/l_car, 'g', 'LineWidth', 5.);
hold on;
for k=22:26
    t_k = (k-1)*K+1:k*K+1;
    t_k = 0.1*t_k;
    plot(t_k, (x_rel(1,:,k))/l_car, 'r', 'LineWidth', 5, 'MarkerSize', 3);
    hold on;
    plot(t_k, (x_rel_n(1,:,k))/l_car, 'b', 'LineWidth', 5, 'MarkerSize', 3);
    plot(0.1*(K*(k-1)+1)*ones(1,41), 3.:0.1:7, 'k', 'LineWidth', 4)
    plot(0.1*(K*(k-1)+9)*ones(1,41), 3.:0.1:7, 'k:', 'LineWidth', 4)
    plot(0.1*(K*k+1)*ones(1,41), 3.:0.1:7, 'k', 'LineWidth', 4)
end
text(35.7,0.9,['Collision'],'HorizontalAlignment','right', 'Fontsize', 28, 'Color', [0.6 0.6 0.6]);
plot(ttt, 0.59*ones(1,length(ttt)), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 5.);
xlabel('time (s)')
ylabel('$x_{\textrm{rel}}/l_{\textrm{car}}$', 'Interpreter', 'latex')
xlim([30,40.5])
ylim([0,8])
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30);
legend('Human', 'RS-IRL, T=2', 'RN-IRL, T=2', 'Decision time', 'Disturbance time', 'Location', 'southeast');

%% Sequence 3
l_car = 4.2; x_rel_actual = x_expert(6,:) - x_expert(1,:);
figure()
ttt = 190:250; ttt = 0.1*ttt;
plot(ttt, (x_rel_actual(190:250))/l_car, 'g', 'LineWidth', 5.);
hold on;
for k=14:16
    t_k = (k-1)*K+1:k*K+1;
    t_k = 0.1*t_k;
    plot(t_k, (x_rel(1,:,k) + 0.2)/l_car, 'r', 'LineWidth', 5, 'MarkerSize', 3);
    hold on;
    plot(t_k, (x_rel_n(1,:,k))/l_car, 'b', 'LineWidth', 5, 'MarkerSize', 3);
    plot(0.1*(K*(k-1)+1)*ones(1,51), 4.:0.1:9, 'k', 'LineWidth', 4)
    plot(0.1*(K*(k-1)+9)*ones(1,51), 4.:0.1:9, 'k:', 'LineWidth', 4)
    plot(0.1*(K*k+1)*ones(1,51), 4.:0.1:9, 'k', 'LineWidth', 4)
end
text(22,0.9,['Collision'],'HorizontalAlignment','right', 'Fontsize', 28, 'Color', [0.6 0.6 0.6]);
plot(ttt, 0.59*ones(1,length(ttt)), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 5.);
xlabel('time (s)')
ylabel('$x_{\textrm{rel}}/l_{\textrm{car}}$', 'Interpreter', 'latex')
xlim([19,25])
ylim([0,9.5])
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30);
legend('Human', 'RS-IRL, T=2', 'RN-IRL, T=2', 'Decision time', 'Disturbance time', 'Location', 'southeast');

%% Sequence 4
l_car = 4.2; x_rel_actual = x_expert(6,:) - x_expert(1,:);
figure()
ttt = 730:780; ttt = 0.1*ttt;
plot(ttt, (x_rel_actual(730:780))/l_car, 'g', 'LineWidth', 5.);
hold on;
for k=50:51
    t_k = (k-1)*K+1:k*K+1;
    t_k = 0.1*t_k;
    plot(t_k, (x_rel(1,:,k)+0.2)/l_car, 'r', 'LineWidth', 5, 'MarkerSize', 3);
    hold on;
    plot(t_k, (x_rel_n(1,:,k))/l_car, 'b', 'LineWidth', 5, 'MarkerSize', 3);
    plot(0.1*(K*(k-1)+1)*ones(1,51), 3.:0.1:8, 'k', 'LineWidth', 4)
    plot(0.1*(K*(k-1)+9)*ones(1,51), 3.:0.1:8, 'k:', 'LineWidth', 4)
    plot(0.1*(K*k+1)*ones(1,51), 3.:0.1:8, 'k', 'LineWidth', 4)
end
text(75.5,0.9,['Collision'],'HorizontalAlignment','right', 'Fontsize', 28, 'Color', [0.6 0.6 0.6]);
plot(ttt, 0.59*ones(1,length(ttt)), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 5.);
xlabel('time (s)')
ylabel('$x_{\textrm{rel}}/l_{\textrm{car}}$', 'Interpreter', 'latex')
%xlim([19,25])
ylim([0,10])
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30);
legend('Human', 'RS-IRL, T=2', 'RN-IRL, T=2', 'Decision time', 'Disturbance time', 'Location', 'southeast');


%% Polytope 3 dimensions
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
proj_idx = [2;3];
figure()
P_simp.projection(proj_idx).plot('alpha',0.); hold on;
%Q_w = P_w.projection(proj_idx);
%Q_w.plot('color','red', 'alpha', 0.7); hold on;
delta_y = 0.4797:0.0001:0.8536;
plot(0.0964*ones(1,length(delta_y)), delta_y, 'r', 'LineWidth', 5)
%xlabel('p(1)   (nothing)')
xlabel('p(2)   (acceleration)')
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
legend('RN-IRL, T=2', 'RS-IRL, T=2', 'Location', 'northeast');
xlim([0,52]);
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
%saveas(fig,['Plots/Bar_expected/abs_bar_expected_',name_str],'epsc');

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(1:numExps,err_vx_n/l_car,0.7,'facecolor','b','edgecolor','b','facealpha',0.6); hold on
bar(1:numExps,err_vx/l_car,0.5,'facecolor','r','edgecolor','r','facealpha',0.45);
xlabel('Demonstration Index in Test Trajectory'); 
legend('RN-IRL, T=2', 'RS-IRL, T=2', 'Location', 'northeast');
xlim([0,52]);
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
%saveas(fig,['Plots/Bar_expected/abs_bar_expected_vx_',name_str],'epsc');

%% Full trajectory plot with boxes
l_car = 4.2;
tt = 1:780; K = 15;

x_rel_actual = x_expert(6,:) - x_expert(1,:);
figure()
plot(0.1*tt, (x_rel_actual(tt))/l_car, 'g', 'LineWidth', 5);
hold on;
text(45,0.9,['Collision'],'HorizontalAlignment','right', 'Fontsize', 28, 'Color', [0.6 0.6 0.6]);
plot(0.1*tt, 0.59*ones(1,length(tt)), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 5.);
rectangle('LineStyle',':', 'LineWidth', 4, 'Position',[72.1 3. 3.0 2.]);
text(77.5,2.5,['Sequence 3'],'HorizontalAlignment','right', 'Fontsize', 26);
rectangle('LineStyle',':', 'LineWidth', 4, 'Position',[31.6 3.8 7.6 2.4]);
text(41,3.3,['Sequence 2'], 'LineWidth', 5, 'HorizontalAlignment','right', 'Fontsize', 26);
rectangle('LineStyle',':', 'LineWidth', 4, 'Position',[19.6 4.0 4.0 5.0]);
text(27,9.4,['Sequence 1'], 'LineWidth', 5, 'HorizontalAlignment','right', 'Fontsize', 26);
rectangle('LineStyle',':', 'LineWidth', 4, 'Position',[72.9 5.1 4.5 4.]);
text(77.9,9.6,['Sequence 4'], 'LineWidth', 5, 'HorizontalAlignment','right', 'Fontsize', 26);
xlim([0,78])
ylim([0,12])
xlabel('time (s)')
ylabel('$x_{\textrm{rel}}/l_{\textrm{car}}$', 'Interpreter', 'latex')
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30);

%% Percentage plot
fig = figure('units','normalized','outerposition',[0 0 1 1]);
ind_RS = find(err_x_n(1:end) > err_x(1:end));
ind_RN = find(err_x_n(1:end) < err_x(1:end));
improvement_RS = (err_x_n(ind_RS)'-err_x(ind_RS)')./err_x_n(ind_RS)';
bar(ind_RS,100*improvement_RS, 'r', 'BarWidth', 0.6);
hold on;
improvement_RN = (err_x_n(ind_RN)'-err_x(ind_RN)')./err_x_n(ind_RN)';
bar(ind_RN,100*improvement_RN, 'b', 'BarWidth', 0.3);
xlabel('Demonstration Index in Test Trajectory');
lim_sup = 100*max(max(improvement_RN), max(improvement_RS));
xlim([0,52])
ylim([-lim_sup-5,lim_sup+5])
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
%saveas(fig,['Plots/Bar_expected/bar_expected_percentage_',name_str],'epsc');


%% Percentage plot vx
fig = figure('units','normalized','outerposition',[0 0 1 1]);
ind_RS = find(err_vx_n(1:end) > err_vx(1:end));
ind_RN = find(err_vx_n(1:end) < err_vx(1:end));
improvement_RS = (err_vx_n(ind_RS)'-err_vx(ind_RS)')./err_vx_n(ind_RS)';
bar(ind_RS,100*improvement_RS, 'r', 'BarWidth', 0.6);
hold on;
improvement_RN = (err_vx_n(ind_RN)'-err_vx(ind_RN)')./err_vx_n(ind_RN)';
bar(ind_RN,100*improvement_RN, 'b', 'BarWidth', 0.3);
xlabel('Demonstration Index in Test Trajectory');
lim_sup = 100*max(max(improvement_RN), max(improvement_RS));
xlim([0,52])
ylim([-lim_sup-5,lim_sup+5])
set(findall(gcf,'type','text'),'FontSize',30);set(gca,'FontSize',30)
%saveas(fig,['Plots/Bar_expected/bar_expected_vx_percentage_',name_str],'epsc');











