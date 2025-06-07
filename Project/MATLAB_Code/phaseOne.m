clc;
clear;

%% 1) system variables (state and input) and parameters
syms h1 h2 h3 h4 V1 V2
syms k1 k2 a1 a2 a3 a4 g A1 A2 A3 A4 landa1 landa2 kc

%% 2) nonlinear system equations (differential equations)
dh1_dt = (landa1*k1 * V1) / A1 + a3 * sqrt(2 * g * h3) / A1 - a1 * sqrt(2 * g * h1) / A1;
dh2_dt = (landa2*k2 * V2) / A2 + a4 * sqrt(2 * g * h4) / A2 - a2 * sqrt(2 * g * h2) / A2;
dh3_dt = (1 - landa2) * k2 * V2 / A3 - a3 * sqrt(2 * g * h3) / A3;
dh4_dt = (1 - landa1) * k1 * V1 / A4 - a4 * sqrt(2 * g * h4) / A4;

nonlinear_equations = {dh1_dt, dh2_dt, dh3_dt, dh4_dt};
disp('Nonlinear Equations:');
for i = 1:length(nonlinear_equations)
    disp(nonlinear_equations{i});
end

%% 3) output equations (here we consider h1 and h2 only but h3 and h4 can be observed too as the paper has analays them)
y1 = kc * h1;
y2 = kc * h2;

%% 4) state and input vectors 
x = [h1; h2; h3; h4];
u = [V1; V2];        

%% 5) symbolic Jacobian matrices
A_sym = jacobian([dh1_dt; dh2_dt; dh3_dt; dh4_dt], x);
B_sym = jacobian([dh1_dt; dh2_dt; dh3_dt; dh4_dt], u);
C_sym = jacobian([y1; y2], x);
D_sym = jacobian([y1; y2], u);

disp('Symbolic A Matrix:');
disp(A_sym);
disp('Symbolic B Matrix:');
disp(B_sym);
disp('Symbolic C Matrix:');
disp(C_sym);
disp('Symbolic D Matrix:');
disp(D_sym);

%% 6) numerical values for the parameters
k1_val    = 3.33e-5;
k2_val    = 3.35e-5;
a1_val    = 0.071e-4;
a2_val    = 0.057e-4;
a3_val    = 0.071e-4;
a4_val    = 0.057e-4;
A1_val    = 28e-4;
A2_val    = 32e-4;
A3_val    = 28e-4;
A4_val    = 32e-4;
g_val     = 9.81;
landa1_val= 0.7;
landa2_val= 0.6;
kc_val    = 1;

%% 7) Substitute numerical parameter values into the differential equations
dh1_dt = subs(dh1_dt, {k1, k2, a1, a2, a3, a4, A1, A2, A3, A4, g, landa1, landa2, kc}, ...
                     {k1_val, k2_val, a1_val, a2_val, a3_val, a4_val, A1_val, A2_val, A3_val, A4_val, g_val, landa1_val, landa2_val, kc_val});
dh2_dt = subs(dh2_dt, {k1, k2, a1, a2, a3, a4, A1, A2, A3, A4, g, landa1, landa2, kc}, ...
                     {k1_val, k2_val, a1_val, a2_val, a3_val, a4_val, A1_val, A2_val, A3_val, A4_val, g_val, landa1_val, landa2_val, kc_val});
dh3_dt = subs(dh3_dt, {k1, k2, a1, a2, a3, a4, A1, A2, A3, A4, g, landa1, landa2, kc}, ...
                     {k1_val, k2_val, a1_val, a2_val, a3_val, a4_val, A1_val, A2_val, A3_val, A4_val, g_val, landa1_val, landa2_val, kc_val});
dh4_dt = subs(dh4_dt, {k1, k2, a1, a2, a3, a4, A1, A2, A3, A4, g, landa1, landa2, kc}, ...
                     {k1_val, k2_val, a1_val, a2_val, a3_val, a4_val, A1_val, A2_val, A3_val, A4_val, g_val, landa1_val, landa2_val, kc_val});

%% 8) Equilibrium Calculation
% i first set the input v1 and v2 to zero but i faced division by zero
% because all the state and input variables became zero
% according to my research the only vay was to fiund another input other
% than zero 

V1_input = 1;
V2_input = 1;

eqns = [ subs(dh1_dt, {V1,V2}, {V1_input,V2_input}) == 0, ...
         subs(dh2_dt, {V1,V2}, {V1_input,V2_input}) == 0, ...
         subs(dh3_dt, {V1,V2}, {V1_input,V2_input}) == 0, ...
         subs(dh4_dt, {V1,V2}, {V1_input,V2_input}) == 0];

% (For square-root expressions we have these conditions h1,h2,h3,h4 > 0)
assume(h1 > 0);
assume(h2 > 0);
assume(h3 > 0);
assume(h4 > 0);

equilibrium_points = solve(eqns, [h1, h2, h3, h4], 'ReturnConditions', true); 
disp('Validity Conditions for the equilibrium solution:');
disp(equilibrium_points.conditions);

% Convert the computed equilibrium points into numerical values
eq_h1 = double(vpa(equilibrium_points.h1, 10));
eq_h2 = double(vpa(equilibrium_points.h2, 10));
eq_h3 = double(vpa(equilibrium_points.h3, 10));
eq_h4 = double(vpa(equilibrium_points.h4, 10));

disp('Equilibrium Points (numerical):');
disp(['h1 = ' num2str(eq_h1)]);
disp(['h2 = ' num2str(eq_h2)]);
disp(['h3 = ' num2str(eq_h3)]);
disp(['h4 = ' num2str(eq_h4)]);

x_eq = [eq_h1; eq_h2; eq_h3; eq_h4];
u_eq = [V1_input; V2_input];

%% 9) Substitute numeric equilibrium and parameter values into A, B, C, D
subsList = { h1,  h2,  h3,  h4,  k1,    V1,       V2,      k2,   a1,    a2,    a3,    a4,    g,     A1,    A2,    A3,    A4,   landa1,   landa2, kc };
valuesList = { eq_h1, eq_h2, eq_h3, eq_h4, k1_val, V1_input, V2_input, k2_val, a1_val, a2_val, a3_val, a4_val, g_val, A1_val, A2_val, A3_val, A4_val, landa1_val, landa2_val, kc_val };

A_numeric = double(subs(A_sym, subsList, valuesList));
B_numeric = double(subs(B_sym, subsList, valuesList));
C_numeric = double(subs(C_sym, subsList, valuesList));
D_numeric = double(subs(D_sym, subsList, valuesList));

disp('A Matrix at the equilibrium point:');
disp(A_numeric);
disp('B Matrix at the equilibrium point:');
disp(B_numeric);
disp('C Matrix at the equilibrium point:');
disp(C_numeric);
disp('D Matrix at the equilibrium point:');
disp(D_numeric);

%% 10) Controllability, Observability, and Stability Analysis
Controllability = ctrb(A_numeric, B_numeric);
Observability   = obsv(A_numeric, C_numeric);

disp('Controllability Matrix:');
disp(Controllability);
disp('Rank of Controllability Matrix:');
disp(rank(Controllability));
rank_C = rank(Controllability);
if rank_C == length(x) 
    disp('The system is controllable.'); 
else 
    disp('The system is not controllable.'); 
end 

disp('Observability Matrix:');
disp(Observability);
disp('Rank of Observability Matrix:');
disp(rank(Observability));
rank_O = rank(Observability);
if rank_O == length(x) 
    disp('The system is observable.'); 
else 
    disp('The system is not observable.'); 
end

eig_values = eig(A_numeric);
disp('Eigenvalues of A Matrix:');
disp(eig_values);

if all(real(eig_values) < 0)
    disp('The system is asymptotically stable.');
else
    disp('The system is unstable or marginally stable.');
end

%% 11) Linear Simulation using lsim
sys = ss(A_numeric, B_numeric, C_numeric, D_numeric);
t_sim = 0:1:10000;  
u_sim = ones(length(t_sim), 2); 
[y_linear, t_linear] = lsim(sys, u_sim, t_sim);

figure('Color', 'w','Position', [100, 100, 800, 600]);

subplot(2,1,1);

yyaxis left
plot(t_linear, y_linear(:,1), 'r-', 'LineWidth', 2);
hold on;
plot(t_linear, y_linear(:,2), 'b-', 'LineWidth', 2);
ylabel('Tank Levels (h1 & h2)', 'FontSize', 14);
yyaxis right
pump_input = @(t) 1 + 0.4*(t>=1);
plot(t_linear, pump_input(t_linear), 'k--', 'LineWidth', 2);
ylabel('Pump Input', 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 14);
title('Linear System Step Response with Pump Input', 'FontSize', 16);
legend('h1 (Tank 1)', 'h2 (Tank 2)', 'Pump Input', 'Location', 'best');
grid on;
set(gca, 'FontSize', 12);

%% 12) Nonlinear Simulation using ode45
pump_input = @(t) 1 + 0.4*(t>=1);
nonlinear_ode = @(t, x) [ (landa1_val*k1_val * pump_input(t)) / A1_val + a3_val * sqrt(2*g_val*x(3)) / A1_val - a1_val * sqrt(2*g_val*x(1)) / A1_val;
                         (landa2_val*k2_val * pump_input(t)) / A2_val + a4_val * sqrt(2*g_val*x(4)) / A2_val - a2_val * sqrt(2*g_val*x(2)) / A2_val;
                         (1 - landa2_val)* k2_val * pump_input(t) / A3_val - a3_val * sqrt(2*g_val*x(3)) / A3_val;
                         (1 - landa1_val)* k1_val * pump_input(t) / A4_val - a4_val * sqrt(2*g_val*x(4)) / A4_val ];
                     
initial_conditions = x_eq;
[t_nonlinear, y_nonlinear] = ode45(nonlinear_ode, t_sim, initial_conditions);

subplot(2,1,2);
yyaxis left
plot(t_nonlinear, y_nonlinear(:,1), 'r-', 'LineWidth', 2);
hold on;
plot(t_nonlinear, y_nonlinear(:,2), 'b-', 'LineWidth', 2);
ylabel('Tank Levels (h1 & h2)', 'FontSize', 14);
yyaxis right
plot(t_nonlinear, pump_input(t_nonlinear), 'k--', 'LineWidth', 2);
ylabel('Pump Input', 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 14);
title('Nonlinear System Step Response with Pump Input', 'FontSize', 16);
legend('h1 (Tank 1)', 'h2 (Tank 2)', 'Pump Input', 'Location', 'best');
grid on;
set(gca, 'FontSize', 12);