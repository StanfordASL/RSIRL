dt = 0.1;
dynamics.f = @(x) [x(1)+x(4)*cos(x(3))*dt+(x(4)^2*sin(x(3))*tan(x(5))/l)*0.5*(dt^2);
    x(2)+x(4)*sin(x(3))*dt+(-x(4)^2*cos(x(3))*tan(x(5))/l)*0.5*(dt^2);
    x(3)+(1/l)*(-x(4)*tan(x(5)))*dt;
    x(4);
    x(5)];
dynamics.g = @(x) [cos(x(3))*0.5*(dt^2),0;
    sin(x(3))*0.5*(dt^2),0;
    -(tan(x(5))/l)*0.5*(dt^2),-(x(4)/l)*(sec(x(5))^2)*0.5*(dt^2);
    dt,0;
    0,dt];

dynamics.A = [1. dt 0 0 0; 0 1 0 0 0; 0 0 1 dt dt^2/2; 0 0 0 1 dt; 0 0 0 0 1];
dynamics.B = [dt^2/2 0; dt 0; 0 dt^3/6; 0 dt^2/2; 0 dt];