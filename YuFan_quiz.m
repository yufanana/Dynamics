clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% part 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial conditions
initial_x = -0.02;
initial_z = 0.01;
initial_xdot = 0;
initial_zdot = 0;
initial_conditions = [initial_x initial_z initial_xdot initial_zdot];

kom = 3;
t = 0:0.1:300;

[t,z] = ode45(@(t,z)ode_system(t,z,kom), t, initial_conditions, kom);

% plot x and z vs time
figure;
subplot(2,1,1); plot(t, z(:,1));
xlabel('t'); ylabel('x');
subplot(2,1,2); plot(t, z(:,2));
xlabel('t'); ylabel('z');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% part 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = kom*[2 1; 1 4];
[V,D] = eig(A);
w = sqrt(diag(D));

figure;
plot([1 2], V(:,1), 'o--', [1 2], V(:,2), 'o--');
xlim([0.5 2.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% part 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A0 = 0.8;
gamma = 0.3;

n = 100;
wd_vec = linspace(0, 6, n);
xMax = zeros(n,1);
zMax = zeros(n,1);

for k = 1:n
    wd = wd_vec(k);
    param = [kom A0 gamma wd];
    
    [t,z] = ode45(@(t,z)driven_system(t,z,param), t, initial_conditions, param);
    
    ll = max(size(t));
    z_eval = z(round(0.8*ll):ll, :);
    xMax(k) = max(z_eval(:,1));
    zMax(k) = max(z_eval(:,2));    
end

% plot xMax and zMax vs wd
figure;
plot(wd_vec, xMax, wd_vec, zMax);
xlabel('wd'); ylabel('Max');
legend('x','z');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function dzdt = ode_system(t,z,kom)
    dzdt_1 = z(3);
    dzdt_2 = z(4);
    dzdt_3 = -kom*((1+z(2))^2 + (1+z(1))^4 - 2);
    dzdt_4 = -kom*(2-exp(-2*z(2)) - exp(-z(1)));
    
    dzdt = [dzdt_1; dzdt_2; dzdt_3; dzdt_4];
end

function dzdt = driven_system(t,z,param)
    kom = param(1);
    A0 = param(2);
    gamma = param(3);
    wd = param(4);
    
    dzdt_1 = z(3);
    dzdt_2 = z(4);
    dzdt_3 = -gamma*z(3)/2 -kom*(4*z(1)+z(2));
    dzdt_4 = -gamma*z(4)/2 -kom*(2*z(2)+z(1)) + A0*sin(3*wd*t/2);
    
    dzdt = [dzdt_1; dzdt_2; dzdt_3; dzdt_4];
end