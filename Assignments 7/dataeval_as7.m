clc;
clear all;
format long;

A = importdata('termodyn_Nsteps=20000_dt=0.000500.txt');

B = importdata("termodyn_fluc_Nsteps=40000_dt=0.000500.txt");

Nb = length(B(:,1))

E = B(:,6)
K = A(:,5)

plot(1:Nb,E)

deltaE2 = std(E)*std(E)

deltaK2 = std(K)^2

CVmicro = (2/2)*(1-(2/(2*144))*deltaK2/(144))^(-1) %C_v/ (N * (f/2)*k_B)

%http://www.nyu.edu/classes/tuckerman/stat.mech/lectures/lecture_6/node2.html