clear all;
clc;


A = importdata('problem_2_2.txt');

plot(A(:,2),A(:,3))
legend('sphase space coordinate','location','northeast')
legend('boxoff')
