close all;
clc;
clear all;
format long;
A = importdata("UvsN_markov.txt");

plot(A(:,1),A(:,3));
A(1,1);
N=A(:,1);
U =A(:,3);

N=N(1:100:1e4);
U=U(1:100:1e4);

B = zeros(4,2);
B(1,1) = 1e1;
B(2,1) = 1e2;
B(3,1) = 1e3;
B(4,1) = 1e4;

B(1,2) = A(B(1,1),2);
B(2,2) = A(B(2,1),2);
B(3,2) = A(B(3,1),2);
B(4,2) = A(B(4,1),2);

pi_vec=zeros(1,4);
pi_vec=pi_vec+pi;
plot(B(:,1),B(:,2),B(:,1),pi_vec);
legend("\pi_{approx}","\pi",'location','southeast')
figure
plot(N,U)
legend('|\pi_{approx}-\pi|')

C = importdata("delta_opt.txt");
figure
plot(C(:,1),C(:,3))
legend('|\pi_{approx}-\pi| vs delta with N = 1e4')

[M,I] = min(C(:,3))
C(I,1)

D = C(:,3)
D = smooth(smooth(smooth(D)))
figure
plot(C(:,1),D)
legend('|\pi_{approx}-\pi| vs delta with N = 1e4')

[M,I] = min(D)
C(I,1)

