name_str = 'p2';
load(['perf_N_2_',name_str,'.mat']);

delta_x = mean( (err_x_n-err_x)./err_x_n );
delta_y = mean( (err_y_n-err_y)./err_y_n );
delta_vx = mean( (err_vx_n-err_vx)./err_vx_n );
delta_vy = mean( (err_vy_n-err_vy)./err_vy_n );

disp(['delta x: ',num2str(delta_x)]);
disp(['delta y: ',num2str(delta_y)]);
disp(['delta vx: ',num2str(delta_vx)]);
disp(['delta vy: ',num2str(delta_vy)]);