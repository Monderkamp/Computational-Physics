clear all;
clc;
close all;

density = importdata('dens_mu=1.000000.txt');
plot(density(:,1),density(:,2))
legend('density profile for \mu = 1','location','northeast');

figure 
N_evol = importdata('N_evolution_mu=1.000000.txt');
plot(N_evol(:,1),N_evol(:,2))
legend('N(t) for \mu = 5','location','northeast');


figure
density = importdata('dens_mu=5.000000.txt');
plot(density(:,1),density(:,2))
legend('density profile for \mu = 5','location','northeast');

figure 
N_evol = importdata('N_evolution_mu=5.000000.txt');
plot(N_evol(:,1),N_evol(:,2))
legend('N(t) for \mu = 5','location','northeast');


figure
density = importdata('dens_mu=20.000000.txt');
plot(density(:,1),density(:,2))
legend('density profile for \mu = 20','location','northeast');

figure 
N_evol = importdata('N_evolution_mu=20.000000.txt');
plot(N_evol(:,1),N_evol(:,2))
legend('N(t) for \mu = 20','location','northeast');

