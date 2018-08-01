function dx_dt  = pvtol_dyn(t,x,u)

mass = 0.486;
J = 0.00383;
g = 9.81;
len = 0.25;

phi = wrapToPi(x(3));
phi_dot = x(6);
vx = x(4); vy = x(5); 

x_dot = vx*cos(phi) - vy*sin(phi);
y_dot = vx*sin(phi) + vy*cos(phi);

f = [x_dot;
     y_dot;
     phi_dot;
     vy*phi_dot - g*sin(phi);
    -vx*phi_dot - g*cos(phi);
     0];
 
B = [zeros(4,2);
     (1/mass), (1/mass);
     (len/J), -(len/J)];
 
dx_dt = f + B*u;

end