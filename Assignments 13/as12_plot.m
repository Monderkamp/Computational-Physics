clear all;
clc;
close all;

B = importdata('init_Nsweeps=10000_delta=0.100000.txt');

scatter(B(:,1),B(:,2),'.');
legend('initial position','location','northeastoutside')
axis([0 14 0 14])
daspect([1 1 1])

figure
B = importdata('finalpos_Nsweeps=10000_delta=0.100000.txt');

scatter(B(:,1),B(:,2),'.');
legend('final position','location','northeastoutside')
axis([0 14 0 14])
daspect([1 1 1])


T(1) = 0.001;
mu(1) = 0.000286533;

T(2) = 1;
mu(2) = -0.0787361;

T(3) = 1.5;
mu(3) = -0.192195;

T(4) = 2;
mu(4) = -0.326766;

T(5) = 2.5;
mu(5) = -0.501348;

T(6) = 3;
mu(6) = -0.67198;

T(7) = 3.5;
mu(7) = -0.88321;

T(8) = 4;
mu(8) = -1.01919;

T(9) = 6;
mu(9) = -1.92522;

T(10) = 8;
mu(10) = -2.73407;

T(11) = 9;
mu(11) =-3.52466;

T(12) = 10;
mu(12) = -4.12948;

figure
plot(T,mu)
legend('\mu(T)','location','southwest')

mu_ex(1,1) = 1;
mu_ex(1,2) = 1.12184;


mu_ex(2,1) = 1.5;
mu_ex(2,2) = 1.60867;

mu_ex(3,1) = 2;
mu_ex(3,2) = 2.07439;

mu_ex(4,1) = 2.5;
mu_ex(4,2) = 2.5001;

mu_ex(5,1) = 3;
mu_ex(5,2) = 2.92975;

mu_ex(5,1) = 3.5;
mu_ex(5,2) = 3.31881;

mu_ex_fit = fit(mu_ex(:,1),mu_ex(:,2),'poly1');
mu_ex_fit_coeffs = coeffvalues(mu_ex_fit);
f = @(x) mu_ex_fit_coeffs(2) + mu_ex_fit_coeffs(1) * x;
a = min(mu_ex(:,1)):0.5:max(mu_ex(:,1));

figure
plot(mu_ex(:,1),mu_ex(:,2),a,f(a))
legend('\mu_{ex}(T)','location','southeast')

disp('a * x + b')

a = mu_ex_fit_coeffs(1)
b = mu_ex_fit_coeffs(2)

f_mu = @(T) -T*log(21.9089*21.9089/(144+1)) + f(T)