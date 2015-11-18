T = linspace(0, 4, 1001);
a2 = 17.0;
a3 = 0.8;
d3 = 5.5;
d4 = 17.1;

J1 = @(x) 8.4375 * x.^2 -5.625 * x.^3;
J2 = @(x) 4.6875 * x.^2 -3.125 * x.^3;
J3 = @(x) 1.875 * x.^2 +1.25 * x.^3;
J4 = @(x) 5.625 * x.^2 -3.75 * x.^3;
J5 = @(x) 16.875 * x.^2 -11.25 * x.^3;
J6 = @(x) 2.8125 * x.^2 -1.875 * x.^3;

j1 = J1(T);
j2 = J2(T);
j3 = J3(T);
j4 = J4(T);
j5 = J5(T);
j6 = J6(T);

jj1 =  2 * 8.4375 * T - 3 * 1.40625 * T .^ 2;
jj2 =  2 * 4.6875 * T - 3 * 0.78125 * T .^ 2;
jj3 = -2 * 1.875  * T + 3 * 0.3125  * T .^ 2;
jj4 =  2 * 5.625  * T - 3 * 0.9375  * T .^ 2;
jj5 =  2 * 16.875 * T - 3 * 2.8125  * T .^ 2;
jj6 =  2 * 2.8125 * T - 3 * 0.46875 * T .^ 2;

jjj1 =  2 * 8.4375 - 6 * 1.40625 * T;
jjj2 =  2 * 4.6875 - 6 * 0.78125 * T;
jjj3 = -2 * 1.875  + 6 * 0.3125  * T;
jjj4 =  2 * 5.625  - 6 * 0.9375  * T;
jjj5 =  2 * 16.875 - 6 * 2.8125  * T;
jjj6 =  2 * 2.8125 - 6 * 0.46875 * T;

px = cosd(j1).*(a2.*cosd(j2) + a3.*cosd(j2+j3) - d4.*sind(j2+j3)) - d3.*sind(j1);
py = sind(j1).*(a2.*cosd(j2) + a3.*cosd(j2+j3) - d4.*sind(j2+j3)) + d3.*cosd(j1);
pz = -a3.*sind(j2+j3) - a2.*sind(j2) - d4.*cosd(j2+j3);

subplot(3,1,1)
plot(T, px)
title('Joint Position in X direction')
xlabel('Time (s)')
ylabel('Position (in)')

subplot(3,1,2)
plot(T, py)
title('Joint Position in Y direction')
xlabel('Time (s)')
ylabel('Position (in)')

subplot(3,1,3)
plot(T, pz)
title('Joint Position in Z direction')
xlabel('Time (s)')
ylabel('Position (in)')

figure

subplot(3,1,1)
hold on
plot(T, j1)
plot(T, j2)
plot(T, j3)
plot(T, j4)
plot(T, j5)
plot(T, j6)
hold off
title('Joint Positions vs Time')
xlabel('Time (s)')
ylabel('$\Theta$','interpreter','latex')

subplot(3,1,2)
hold on
plot(T, jj1)
plot(T, jj2)
plot(T, jj3)
plot(T, jj4)
plot(T, jj5)
plot(T, jj6)
hold off
title('Joint Velocities vs Time')
xlabel('Time (s)')
ylabel('$\dot{\Theta}$','interpreter','latex')
legend('\theta_1','\theta_2','\theta_3','\theta_4','\theta_5','\theta_6')

subplot(3,1,3)
hold on
plot(T, jjj1)
plot(T, jjj2)
plot(T, jjj3)
plot(T, jjj4)
plot(T, jjj5)
plot(T, jjj6)
hold off
title('Joint Accelerations vs Time')
xlabel('Time (s)')
ylabel('$$\ddot{\Theta}$$','Interpreter','latex')

