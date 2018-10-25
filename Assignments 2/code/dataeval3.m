clear all;
clc;


A1 = importdata('root1.txt');
A2 = importdata('root2.txt');
A3 = importdata('root3.txt');
A4 = importdata('root4.txt');

plot(A1(:,1),A1(:,2),A2(:,1),A2(:,2),A3(:,1),A3(:,2),A4(:,1),A4(:,2))
legend('distance to root -.86','distance to root -.34','distance to root .34','distance to root .86','location','northeast')
legend('boxoff')
