clear all; close all; clc;

%% Measure the number of balanced and imbalanced triads in a network

A = randi([-1 1],6);

A = triu(A,1) + triu(A,1)';  %symmetric matrix

figure();
imagesc(A);

P = 1*(A>0);  % just the positive edges
N = -1*(A<0);  % just the negative edges


%% what do the entries in these matrices represent?
P*P
P*P*P
N*P
N*P;
N*N*P
N*N*N
%%
clear all; close all; clc;
%% Random matrix measure of imbalance

sz = 5
stat_i=[];
stat_f=[];
% for i = 1:200
%     imb = [];
    for j = 1:100
        X0 = 0.3*randn(sz);
        X0 = triu(X0,1) + triu(X0,1)';  % initial condition connectivity matrix (symmetric)
        X0_vec = reshape(X0,[],1);
        P0 = 1*(X0>0); 
        N0 = 1*(X0<0); 
        PPPi = P0*P0*P0;
        NNNi = N0*N0*N0;
        NNPi = N0*N0*P0;
        NPPi = N0*P0*P0;
        X0_imb = sum(diag(NPPi) + diag(NNNi))/sum(diag(NPPi) + diag(NNNi)+ diag(PPPi) + diag(NNPi));
        X0_bal(j) = 1-X0_imb;
        
        tspan = 0:0.01:2;  %timespan for ode integration
        [t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz), tspan, X0_vec);
        
        X_final = reshape(X_vec(end,:),sz,sz); %look at the final connectivity matrix
        P = 1*(X_final>0);  % just the positive edges
        N = 1*(X_final<0);
        PPP = P*P*P;
        NNN = N*N*N;
        NNP = N*N*P;
        NPP = N*P*P;
        Xf_imb = sum(diag(NPP) + diag(NNN))/sum(diag(NPP) + diag(NNN)+ diag(PPP) + diag(NNP));
        Xf_bal(j) = 1-X0_imb;
    end
%     %              mean imbalance    var imbalance    mean Pii          mean P^3ii          mean N^3ii         var N^3ii
%     stat_i(i,:) = [mean(imb(:,1)), var(imb(:,1)), mean(probs_i(:,1)), mean(probs_i(:,2)), mean(probs_i(:,3)), var(probs_i(:,3)), mean(pr obs_i(:,4)), mean(probs_i(:,5))];
%     stat_f(i,:) = [mean(imb(:,2)), var(imb(:,2)), mean(probs_f(:,1)), mean(probs_f(:,2)), mean(probs_f(:,3)), var(probs_f(:,3)), mean(probs_f(:,4)), mean(probs_f(:,5))];
%     bal(i,:) = [mean(X0_bal), mean(Xf_bal)];
% end
%%
subplot(1,2,1)
h = histfit(X0_bal);
title("balance measure of 50 random matrices of sz 5 (initial)");
subplot(1,2,2)
h = histfit(Xf_bal);
title("balance measure of 50 random matrices of sz 5");
%%
[values1, edges1] = histcounts(stat_i(:,1), 'Normalization', 'probability');
centers1 = (edges1(1:end-1)+edges1(2:end))/2;
[values2, edges2] = histcounts(stat_f(:,1), 'Normalization', 'probability');
centers2 = (edges2(1:end-1)+edges2(2:end))/2;
figure
subplot(1,2,1)
%hist(stat_i(:,1)); hold on;
plot(centers1, values1, 'k-')
xlim([-100 100]);
title('histogram of mean(X0)')
subplot(1,2,2)
%hist(stat_f(:,1)); hold on;
plot(centers2, values2, 'k-')
title('histogram of mean(X final)')

%%
for i=1:8
    figure()
    hist(stat_i(:,i));
end
%%
figure()
hist(bal(:,1));
figure()
hist(bal(:,2));
%% functions

function dXdt = ode_struc_bal(t,X,sz)

X = reshape(X,sz,sz);  %must reshape
dXdt = X^2;

dXdt = reshape(dXdt,[],1);
end