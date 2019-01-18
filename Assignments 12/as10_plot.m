
clear all;
clc;
close all;

B = importdata('init_conf.txt');

figure
scatter(B(:,2),B(:,3),'.');
legend('initial position','location','northeastoutside')
axis([0 14 0 14])
daspect([1 1 1])

msd = importdata('msd_Nsweeps=10000_delta=0.100000.txt');
msd(:,1) = msd(:,1)/144;
linearfit = fit(msd(750:1029,1),msd(750:1029,2),'poly1')

linearfit_coeffs = coeffvalues(linearfit)
linearfit_function = @(x) linearfit_coeffs(1)*x + linearfit_coeffs(2)


figure 
plot(msd(:,1),msd(:,2),msd(:,1),linearfit_function(msd(:,1)))
legend('msd','location','southeast')

D = linearfit_coeffs(1)/4



%{
msd_ot = msd(:,2)./(msd(:,1))
figure 
plot(msd(:,1),msd_ot)
legend('msd_{ot}','location','southeast')
%}



% finalpos = importdata('finalpos_Nsweeps=100_delta=0.200000.txt');
% figure
% scatter(finalpos(:,1),finalpos(:,2),'.');
% legend('final position','location','northeastoutside')
% axis([0 14 0 14])
% daspect([1 1 1])

% acc_rate = zeros(4,2);
% 
% acc_rate(1,1) = 0.01; %delta
% acc_rate(1,2) = 0.995556; %rate
% 
% acc_rate(2,1) = 0.1; %delta
% acc_rate(2,2) =  0.693889; %rate
% 
% acc_rate(3,1) = 0.2; %delta
% acc_rate(3,2) = 0.446458; %rate
% 
% acc_rate(4,1) = 0.25; %delta
% acc_rate(4,2) = 0.372639; %rate
% 
% acc_rate(5,1) = 0.5; %delta
% acc_rate(5,2) = 0.140278; %rate
% 
% acc_rate(6,1) = 1; %delta
% acc_rate(6,2) = 0.0417361; %rate
% 
% acc_rate(7,1) = 2; %delta
% acc_rate(7,2) = 0.00722222; %rate
% 
% acc_rate(8,1) = 3; %delta
% acc_rate(8,2) = 0.00326389; %rate

% accrfit = fit(acc_rate(:,1),acc_rate(:,2),'exp1')
% coeffs = coeffvalues(accrfit)

% fit = @(x) coeffs(1)*exp(coeffs(2)*x)

% figure 
% x = 1:1000;
% x = x*3/1000;
% % plot(acc_rate(:,1),acc_rate(:,2),x,fit(x))
% legend('acceptence rate vs delta','exponential fit of acceptence rate vs delta','location','northeast')
% 
% % G = importdata('g_Nsweeps=10000_delta=0.100000.txt');
% figure
% plot(G(:,2),G(:,3))


