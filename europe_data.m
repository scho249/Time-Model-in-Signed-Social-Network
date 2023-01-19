clear; close all; clc;
%%

years = 1873:1:1922; %these are the years for the graph timeseries

load('WW1_cube_timeseries.mat'); %load the graph data stored in mat file

% GP: code names for the 5 great power in Europe - United Kingdom, France,
% Germany, Austria-Hungary, Russia

% mat_cube_ally: alliance values over time range from 0 to 4

%mat_cube_riv: rivalry values over time

% mat_cube_mid: militarized interstate disputes over time, negative values
% (-5 to 5) more negative = more hostile; positiv means cooperating to
% fight other
% are conflicts and positive values are cooperations during these conflicts

% NMC_mat_gp: National Material Capabilities of great powers over time
%% ally + riv + mid
X_ir = [];
for j = 1:50
    X_ir(:,:,j) = 2.5.*mat_cube_ally(:,:,j) + 5.*mat_cube_riv(:,:,j)+ mat_cube_mid(:,:,j);
end

%% ally + riv
X_ar = [];
for j = 1:50
    X_ar(:,:,j) = mat_cube_ally(:,:,j) + 2.*mat_cube_riv(:,:,j);
end
%% mid
X_mid = mat_cube_mid;

%% measure imbalance (europe data)
imb = [];
for j = 1:50
    X0_ir = X_ir(:,:,j);
    X0_ar = X_ar(:,:,j);
    X0_mid = X_mid(:,:,j);
    X0ir_vec = reshape(X0_ir,[],1);
    X0ar_vec = reshape(X0_ar,[],1);
    X0mid_vec = reshape(X0_mid,[],1);
    P0_ir = 1*(X0_ir>0); P0_ar = 1*(X0_ar>0); P0_mid = 1*(X0_mid>0);
    N0_ir = 1*(X0_ir<0); N0_ar = 1*(X0_ar<0); N0_mid = 1*(X0_mid<0);
    PPP_ir = P0_ir*P0_ir*P0_ir; PPP_ar = P0_ar*P0_ar*P0_ar; PPP_mid = P0_mid*P0_mid*P0_mid;
    NNN_ir = N0_ir*N0_ir*N0_ir; NNN_ar = N0_ar*N0_ar*N0_ar; NNN_mid = N0_mid*N0_mid*N0_mid;
    NNP_ir = N0_ir*N0_ir*P0_ir; NNP_ar = N0_ar*N0_ar*P0_ar; NNP_mid = N0_mid*N0_mid*P0_mid;
    NPP_ir = N0_ir*P0_ir*P0_ir; NPP_ar = N0_ar*P0_ar*P0_ar; NPP_mid = N0_mid*P0_mid*P0_mid;
    ir_imb = sum(diag(NPP_ir) + diag(NNN_ir))/sum(diag(NPP_ir) + diag(NNN_ir)+ diag(PPP_ir) + diag(NNP_ir));
    ar_imb = sum(diag(NPP_ar) + diag(NNN_ar))/sum(diag(NPP_ar) + diag(NNN_ar)+ diag(PPP_ar) + diag(NNP_ar));
    mid_imb = sum(diag(NPP_mid) + diag(NNN_mid))/sum(diag(NPP_mid) + diag(NNN_mid)+ diag(PPP_mid) + diag(NNP_mid));
    imb(j,1) = ir_imb;
end

%%
imb = [];
sz = 5;
for j = 1:50
    X0_ir = X_ir(:,:,j);
    X0ir_vec = reshape(X0_ir,[],1);
    tspan = 0:0.01:1;  %timespan for ode integration
    [t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz), tspan, X0ir_vec);
    
    Xf_ir = reshape(X_vec(end,:),sz,sz); %look at the final connectivity matrix
    P0_ir = 1*(Xf_ir>0);
    N0_ir = 1*(Xf_ir<0);
    %         P0_ir = 1*(X0_ir>0);
    %         N0_ir = 1*(X0_ir<0);
    PPP_ir = P0_ir*P0_ir*P0_ir;
    NNN_ir = N0_ir*N0_ir*N0_ir;
    NNP_ir = N0_ir*N0_ir*P0_ir;
    NPP_ir = N0_ir*P0_ir*P0_ir;
    ir_imb = sum(diag(NPP_ir) + diag(NNN_ir))/sum(diag(NPP_ir) + diag(NNN_ir)+ diag(PPP_ir) + diag(NNP_ir));
    eu_bal2(j,1) = 1-ir_imb;
end

%%
subplot(1,2,1)
h = histfit(X0_bal);
title("balance measure of 50 random matrices of sz 5 (initial)");
subplot(1,2,2)
h = histfit(Xf_bal);
title("balance measure of 50 random matrices of sz 5");
%%
subplot(3,1,1)
plot(imb(:,1))
subplot(3,1,2)
plot(imb(:,2))
subplot(3,1,3)
plot(imb(:,3))


%% measure imbalance random matrix

for j = 1:1000
    X0 = randn(5);
    X0 = triu(X0,1) + triu(X0,1)';
    X0_vec = reshape(X0,[],1);
    tspan = 0:0.01:2;  %timespan for ode integration
    [t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz), tspan, X0_vec);
    X = reshape(X_vec(end,:),sz,sz);
    
    P0 = 1*(X0>0);
    N0 = 1*(X0<0);
    PPPi = P0*P0*P0;
    NNNi = N0*N0*N0;
    NNPi = N0*N0*P0;
    NPPi = N0*P0*P0;
    X0_imb = sum(diag(NPPi) + diag(NNNi))/sum(diag(NPPi) + diag(NNNi)+ diag(PPPi) + diag(NNPi));
    bal1(j,1) = 1-X0_imb;
    
    P = 1*(X>0);
    N = 1*(X<0);
    PPP = P*P*P;
    NNN = N*N*N;
    NNP = N*N*P;
    NPP = N*P*P;
    X_imb = sum(diag(NPP) + diag(NNN))/sum(diag(NPP) + diag(NNN)+ diag(PPP) + diag(NNP));
    bal2(j,1) = 1-X_imb;
    
end

%%
subplot(1,2,1)
h = histogram(bal1)
h.FaceColor = [0.6 0.8 1];
hold on 
line([mean(bal1) mean(bal1)], ylim, 'Linewidth', 2, 'Color', 'r');
title("balance measure of 1000 random matrices of sz 5 (initial)");
subplot(1,2,2)
h2 = histogram(bal2)
h2.FaceColor = [0.6 0.8 1];
hold on 
line([mean(bal2) mean(bal2)], ylim, 'Linewidth', 2, 'Color', 'r');
title("balance measure of 1000 random matrices of sz 5 (final)");

%%
p2 = prctile(1-imb, [0:25:100], 'all');
% h2 = histfit(1-imb)
% hold on
% line([mean(eu_bal(1:end-1)) mean(eu_bal(1:end-1))], ylim, 'LineWidth', 2, 'Color', 'r');
% hold on
% line([eu_bal(1) eu_bal(1)], ylim, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
% hold on
% line([eu_bal(end-1) eu_bal(end-1)], ylim, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
% h2(1).FaceColor = [0.6 0.8 1];
% h2(2).Color = [0.1 0.5 0.7];
% xticks([0 0.25 0.3125 0.3333 0.5 0.6875 0.75 1]); xtickangle(90);
% xticklabels(["0", "10", "20", "30", "50", "70", "90", "100"]);
% xlabel('percentile')

%% histogram (percentile)
p = prctile(1-rand_imb, [0:5:100], 'all');
label = compose("%.2f", 0:0.05:1);
h = histfit(1-rand_imb)
hold on
line([mean(1-rand_imb) mean(1-rand_imb)], ylim, 'LineWidth', 2, 'Color', [0.1 0.5 0.7])
% hold on
% line([eu_bal eu_bal], ylim);
h(1).FaceColor = [0.6 0.8 1];
h(2).Color = [0.1 0.5 0.7];
xticks(p);
xticklabels(label); xtickangle(90);
xlabel('percentile')
title('mean balance distribution of random sz 5 matrix');
%% make figures (ally + riv)
sz = 5;
X0_vec = reshape(X_ar(:,:,1),[],1);
tspan = 0:0.01:1;  %timespan for ode integration
[t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz), tspan, X0_vec);

X_final = reshape(X_vec(end,:),sz,sz); %look at the final connectivity matrix

[V,D] = eig(X_final);   %eigenvalues and eigenvectors of final connectivity matrix
[m, i] = max(diag(D));
%timeseries of matrix weights
fig = figure('position', [0, 0, 300, 200]); hold on;
plot(t,X_vec);
%ylim([-10,10]);
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

%% make figures (ally + riv + mid)
sz = 5;
X0_vec = reshape(X_ir(:,:,1),[],1);
tspan = 0:0.005:1;  %timespan for ode integration
[t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz), tspan, X0_vec);

X_final = reshape(X_vec(end,:),sz,sz); %look at the final connectivity matrix

[V,D] = eig(X_final);   %eigenvalues and eigenvectors of final connectivity matrix
[m, i] = max(diag(D));
%timeseries of matrix weights
fig = figure('position', [0, 0, 300, 200]); hold on;
plot(t,X_vec);
%ylim([-10,10]);
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

%% make figures (ally + riv)
sz = 5;
X0_vec = reshape(X_ar(:,:,1),[],1);
tspan = 1:1:50;  %timespan for ode integration
[t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz), tspan, X0_vec);

X_final = reshape(X_vec(end,:),sz,sz); %look at the final connectivity matrix

[V,D] = eig(X_final);   %eigenvalues and eigenvectors of final connectivity matrix
[m, i] = max(diag(D));
%timeseries of matrix weights
fig = figure('position', [0, 0, 300, 200]); hold on;
plot(t,X_vec);
%ylim([-10,10]);
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

%% make figures (ally + riv)
sz = 5;
X0_vec = reshape(X_ar(:,:,1),[],1);
tspan = 1:1:50;  %timespan for ode integration
[t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz), tspan, X0_vec);

X_final = reshape(X_vec(end,:),sz,sz); %look at the final connectivity matrix

[V,D] = eig(X_final);   %eigenvalues and eigenvectors of final connectivity matrix
[m, i] = max(diag(D));
%timeseries of matrix weights
fig = figure('position', [0, 0, 300, 200]); hold on;
plot(t,X_vec);
%ylim([-10,10]);
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

%% make figures (mid)
sz = 5;
X0_vec = reshape(X_mid(:,:,15),[],1);
tspan = 0:0.01:1;  %timespan for ode integration
[t,X_vec] = ode45(@(t,X) ode_struc_bal(t,X,sz), tspan, X0_vec);

X_final = reshape(X_vec(end,:),sz,sz); %look at the final connectivity matrix

[V,D] = eig(X_final);   %eigenvalues and eigenvectors of final connectivity matrix
[m, i] = max(diag(D));
%timeseries of matrix weights
fig = figure('position', [0, 0, 300, 200]); hold on;
plot(t,X_vec);
%ylim([-10,10]);
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