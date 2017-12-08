load('means_karen.mat');

for k=1:5
    plot(0.1*(1:15), means_2(k,16:30), 'LineWidth', 4);
    hold all;
end
xlim([0.1, 1.5]);
xlabel('time (s)')
ylabel('$u_s$', 'Interpreter', 'latex')
set(gca,'FontSize',40);