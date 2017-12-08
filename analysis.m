%% Difference Bar plots

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(1:numExps-1,err_x_n(1:end-1)'-err_x(1:end-1)');
xlabel('Demonstration Index in Test Trajectory');
set(findall(gcf,'type','text'),'FontSize',28);set(gca,'FontSize',28)
saveas(fig,['Data/Inference/Performance/Plots/Bar_expected/bar_expected_',name_str],'epsc');

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(1:numExps,err_soft_x_n'-err_soft_x');
xlabel('Demonstration Index in Test Trajectory');
set(findall(gcf,'type','text'),'FontSize',28);set(gca,'FontSize',28)
saveas(fig,['Data/Inference/Performance/Plots/Bar_most_likely/bar_most_likely_',name_str],'epsc');

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(1:numExps,err_y_n'-err_y');
xlabel('Demonstration Index in Test Trajectory');
set(findall(gcf,'type','text'),'FontSize',28);set(gca,'FontSize',28)
saveas(fig,['Data/Inference/Performance/Plots/Bar_expected/bar_expected_',name_str,'_y'],'epsc');

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(1:numExps,err_soft_y_n'-err_soft_y');
xlabel('Demonstration Index in Test Trajectory');
set(findall(gcf,'type','text'),'FontSize',28);set(gca,'FontSize',28)
saveas(fig,['Data/Inference/Performance/Plots/Bar_most_likely/bar_most_likely_',name_str,'_y'],'epsc');

%% Absolute bar plots

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(1:numExps,err_x_n,0.7,'facecolor','b','edgecolor','b','facealpha',0.6); hold on
bar(1:numExps,err_x,0.5,'facecolor','r','edgecolor','r','facealpha',0.45);
xlabel('Demonstration Index in Test Trajectory'); 
legend('RN-IRL, T=2', 'RS-IRL, T=2', 'Location', 'northwest');
set(findall(gcf,'type','text'),'FontSize',28);set(gca,'FontSize',28)
saveas(fig,['Data/Inference/Performance/Plots/Bar_expected/abs_bar_expected_',name_str],'epsc');

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(1:numExps,err_soft_x_n,0.7,'facecolor','b','edgecolor','b','facealpha',0.6); hold on
bar(1:numExps,err_soft_x,0.5,'facecolor','r','edgecolor','r','facealpha',0.45);
xlabel('Demonstration Index in Test Trajectory'); 
legend('RN-IRL, T=2', 'RS-IRL, T=2', 'Location', 'northwest');
set(findall(gcf,'type','text'),'FontSize',28);set(gca,'FontSize',28)
saveas(fig,['Data/Inference/Performance/Plots/Bar_most_likely/abs_bar_most_likely_',name_str],'epsc');

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(1:numExps,err_y_n,0.7,'facecolor','b','edgecolor','b','facealpha',0.6); hold on
bar(1:numExps,err_y,0.5,'facecolor','r','edgecolor','r','facealpha',0.45);
xlabel('Demonstration Index in Test Trajectory'); 
set(findall(gcf,'type','text'),'FontSize',28);set(gca,'FontSize',28)
saveas(fig,['Data/Inference/Performance/Plots/Bar_expected/abs_bar_expected_',name_str,'_y'],'epsc');

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(1:numExps,err_soft_y_n,0.7,'facecolor','b','edgecolor','b','facealpha',0.6); hold on
bar(1:numExps,err_soft_y,0.5,'facecolor','r','edgecolor','r','facealpha',0.45);
xlabel('Demonstration Index in Test Trajectory'); 
set(findall(gcf,'type','text'),'FontSize',28);set(gca,'FontSize',28)
saveas(fig,['Data/Inference/Performance/Plots/Bar_most_likely/abs_bar_most_likely_',name_str,'_y'],'epsc');

%% Full trajectory plot

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
