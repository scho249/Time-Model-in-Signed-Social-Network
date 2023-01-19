close all; clear all; clc;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%%

sz = [4,10,20,50,100];  % size of the connectivity matrix 4 nodes, so matrix is 4x4

for j = 1:length(sz)
    X0 = 0.3*randn(sz(j));
    X0 = triu(X0,1) + triu(X0,1)';  % initial condition connectivity matrix (symmetric)
    X0_vec = reshape(X0,[],1);
    
    tspan = 0:0.01:2;  %timespan for ode integration
    [t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz(j)), tspan, X0_vec);
    
    X_final = reshape(X_vec(end,:),sz(j),sz(j)); %look at the final connectivity matrix
    X_norm = X_final/norm(X_final);
    [V,D] = eig(X_final);
    [V_n,D_n] = eig(X_norm);
    
    % make figures
    
    %timeseries of matrix weights
%     fig = figure('position', [0, 0, 300, 200]);
%     subplot(1,2,1);
%     plot(t,X_vec);
%     ylim([-10,10]);
%     xlabel('time');ylabel('$X_{ij}$')
%     title(['size ', num2str(sz(j))]);
%     
%     subplot(1,2,2);
%     plot(t,X_vec/norm(X_vec));
%     ylim([-10,10]);
%     xlabel('time');ylabel('$X_{ij}$')
%     title(['size ', num2str(sz(j)),'(normalized)']);
%     
%     % final connectivity matrix values
%     figure();
%     imagesc(X_final);
%     pbaspect([1 1 1]);
%     colorbar();
%     title(['final matrix of size ', num2str(sz(j))]);
%     
%     % eigenvalues of final connectivity matrix and sign of outer product of
%     % leading eigenvector
%     fig = figure('position', [0, 0, 600, 200]); hold on;
%     subplot(1,2,1);
%     plot(diag(D),'bo');
%     title(['final matrix of size ', num2str(sz(j)), ' eigenvalues']);
%     
%     subplot(1,2,2);
%     plot(diag(D_n),'bo');
%     title(['final matrix (normalized) of size ', num2str(sz(j)), ' eigenvalues']);
%     
%     
%     subplot(1,2,2);
%     plot(diag(D),'bo');
%     title(['$u_1$ of size ', num2str(sz(j)), ' outer product']);
%     imagesc(sign(V(:,end)*V(:,end)'));
%     pbaspect([1 1 1]);
%     colorbar();
    
    
end

%% disconnected graph
sz = 6;  % size of the connectivity matrix 4 nodes, so matrix is 4x4

% X0 = 0.3*randn(sz);
%  % initial condition connectivity matrix (symmetric)

X0 = [0,-1,-1,0,0,0;
    -1,0,-1,0,0,0;
    -1,-1,0,1,0,0;
    0,0,1,0,1,1;
    0,0,0,1,0,1;
    0,0,0,1,1,0];
%% disconnected with 1 positive tie between
X0 = [0,1,1,0,0,0;
    1,0,1,0,0,0;
    1,1,0,-1,0,0;
    0,0,-1,0,1,1;
    0,0,0,1,0,1;
    0,0,0,1,1,0];
%% 3 disconneced clusters
sz = 9; 
X0 = [0 1 1 0 0 0 0 0 0;
      1 0 1 0 0 0 0 0 0;
      1 1 0 0 0 0 0 0 0;
      0 0 0 0 1 1 0 0 0;
      0 0 0 1 0 1 0 0 0;
      0 0 0 1 1 0 0 0 0;
      0 0 0 0 0 0 0 1 1;
      0 0 0 0 0 0 1 0 1;
      0 0 0 0 0 0 1 1 0];
%%
X0 = triu(X0,1) + triu(X0,1)';
X0_vec = reshape(X0,[],1);
tspan = 0:0.01:2;  %timespan for ode integration
[t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz), tspan, X0_vec);

X_final = reshape(X_vec(end,:),sz,sz); %look at the final connectivity matrix

[V,D] = eig(X_final);   %eigenvalues and eigenvectors of final connectivity matrix

% make figures

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
imagesc(sign(V(:,1)*V(:,1)'));
pbaspect([1 1 1]);
colorbar();

%% positive Mu
sz = 4;
X0 = 0.3*randn(sz);
X0 = X0 + max(X0(:));
X0 = triu(X0,1) + triu(X0,1)';
X0_vec = reshape(X0,[],1);
tspan = 0:0.01:2;  %timespan for ode integration
[t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz), tspan, X0_vec);

X_final = reshape(X_vec(end,:),sz,sz); %look at the final connectivity matrix

[V,D] = eig(X_final);   %eigenvalues and eigenvectors of final connectivity matrix

% make figures

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
imagesc(sign(V(:,end)*V(:,end)'));
pbaspect([1 1 1]);
colorbar();

%% negative mu
sz = 4;
X0 = 0.3*randn(sz);
X0 = X0 - max(X0(:));
X0 = triu(X0,1) + triu(X0,1)';
X0_vec = reshape(X0,[],1);
tspan = 0:0.01:2;  %timespan for ode integration
[t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz), tspan, X0_vec);

X_final = reshape(X_vec(end,:),sz,sz); %look at the final connectivity matrix

[V,D] = eig(X_final);   %eigenvalues and eigenvectors of final connectivity matrix

% make figures

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
imagesc(sign(V(:,end)*V(:,end)'));
pbaspect([1 1 1]);
colorbar();

%%  part 4
clear; close all; clc;
sz = 6;  % size of the connectivity matrix 4 nodes, so matrix is 4x4

% X0 = [0,-1,-1,1,1,1;
%     -1,0,-1,1,1,1;
%     -1,-1,0,1,1,1;
%     1,1,1,0,1,1;
%     1,1,1,1,0,1;
%     1,1,1,1,1,0];
% for i = 1:sz
%     X0(:,i) = 0.3*rand*X0(:,i);
% end
% X0 = [-1;-1;1;1;-1;1]*diag(0.3*rand)*[-1;-1;1;1;-1;1].'
X0 = triu(X0,1) + triu(X0,1)';
X0_vec = reshape(X0,[],1);
tspan = 0:0.01:2;  %timespan for ode integration
[t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz), tspan, X0_vec);

X_final = reshape(X_vec(end,:),sz,sz); %look at the final connectivity matrix

[V,D] = eig(X_final);   %eigenvalues and eigenvectors of final connectivity matrix

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
imagesc(sign(V(:,end)*V(:,end)'));
pbaspect([1 1 1]);
colorbar();

%% eigenvector


%% functions

function dXdt = ode_struc_bal(t,X,sz)

X = reshape(X,sz,sz);  %must reshape
dXdt = X^2;

dXdt = reshape(dXdt,[],1);
end