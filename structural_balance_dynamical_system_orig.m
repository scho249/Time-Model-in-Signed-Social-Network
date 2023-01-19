close all; clear all; clc;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%%

sz = 6;  % size of the connectivity matrix 4 nodes, so matrix is 4x4

X0 = 0.3*randn(sz);
 % initial condition connectivity matrix (symmetric)
X0 = triu(X0,1) + triu(X0,1)';
X0_vec = reshape(X0,[],1);
tspan = 0:0.01:2;  %timespan for ode integration
[t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz), tspan, X0_vec);

X_final = reshape(X_vec(end,:),sz,sz); %look at the final connectivity matrix

[V,D] = eig(X_final);   %eigenvalues and eigenvectors of final connectivity matrix
[m, i] = max(diag(D));
%% make figures

%timeseries of matrix weights
fig = figure('position', [0, 0, 300, 200]); hold on;
plot(t,X_vec);
ylim([-10,10]);
xlabel('time');
ylabel('$X_{ij}$')

% final connectivity matrix values
figure();
imagesc(X_final);
pbaspect([1 1 1]);
colorbar();
title('final matrix');

% eigenvalues of final connectivity matrix and sign of outer product of
% leading eigenvector
fig = figure('position', [0, 0, 600, 200]); hold on;
subplot(1,2,1);
plot(diag(D),'bo');
title('eigenvalues');


subplot(1,2,2);
plot(diag(D),'bo');
title('$u_1$ outer product');
imagesc(sign(V(:,i)*V(:,i)'));
pbaspect([1 1 1]);
colorbar();




%% functions

function dXdt = ode_struc_bal(t,X,sz)

    X = reshape(X,sz,sz);  %must reshape
    dXdt = X^2;

    dXdt = reshape(dXdt,[],1);
end