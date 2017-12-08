P_simp = Polyhedron('V',eye(4));
proj_idx = [1;2;3];
figure()
P_simp.projection(proj_idx).plot('alpha',0.); hold on;
P_w.projection(proj_idx).plot('color','red', 'alpha', 0.7); hold on;
xlabel('p(1)')
ylabel('p(2)')
zlabel('p(3)')
set(gca,'FontSize', 30);
axis equal;
grid off;