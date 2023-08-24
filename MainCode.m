close all;
clear;
%% Initial Conditions
u0 = 21.04; v0 = 0; w0 = -0.9114; % m/s
x0 = 0; y0 = 0; z0 = 0; % m
p0 = 0; q0 = 0.001; r0 = 0; % rad/s
phi0 = 0; theta0 = 0; psi0 = 0; % rad

states0 = [u0, v0, w0, x0, y0, z0, p0, q0, r0, phi0, theta0, psi0];

[t,y] = rk4(@six_dof_model, 110, 0.01, states0);

%% calculate control inputs
c = zeros(numel(t), 4);
for i=1:numel(t)
    cc = cont(t(i), y(i));
    c(i,:) = cc;
end

figure();
subplot(5,2,1);
plot(t,y(:,1));
ylabel('u');

subplot(5,2,3);
plot(t,y(:,2));
ylabel('v');

subplot(5,2,5);
plot(t,y(:,3));
ylabel('w');

subplot(5,2,2);
plot(t,y(:,8));
ylabel('q');

subplot(5,2,4);
plot(t,y(:,9));
ylabel('r');

subplot(5,2,7);
plot(t,y(:,11));
ylabel('\theta');

subplot(5,2,9);
plot(t,y(:,12));
ylabel('\psi');
ylim([-0.08 0.08])

subplot(5,2,6);
plot(t,c(:,1));
ylabel('\delta e');

subplot(5,2,8);
plot(t,c(:,2));
ylabel('\delta a');

subplot(5,2,10);
plot(t,c(:,3));
ylabel('\delta r');

% V = sqrt((y(:,1).^2)+(y(:,2).^2)+(y(:,3).^2));
% figure()
% plot(t,V,'r-.')
