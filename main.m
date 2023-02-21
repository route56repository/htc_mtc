%% --------------------------------------------------------------------- %%
%                 Energy Driven Transmission Schemes for                  %
%                     Coexistence between HTC and MTC                     %
%% --------------------------------------------------------------------- %%
   
clear; close all;
delete(gcp('nocreate'));

% cvx_setup;
% parpool;

rng(1);

%% Simulation Options
densities   = input('Select sensor densities, case i (1) or case ii (2): ');
duplex      = input('Select UTs communication, half duplex (1) or full duplex (2): ');

%% System Parameters
% 1. General
T_f         = 1e0;                                                         % Frame duration
T_d         = T_f/2;                                                       % DL
T_u         = T_f - T_d;                                                   % UL
N_s         = 10;                                                          % Number of slots
T_s         = T_u/N_s;                                                     % Slot duration          

% 2. BS
M           = 2e1;                                                         % Number of antennas
f_o         = 1.9e9;                                                       % Carrier frequency
c_o         = 3e8;                                                         % Light speed
lambda_o    = c_o/f_o;                                                     % Wavelength
K_o         = 2*pi/lambda_o;                                               % Wavenumber
D           = lambda_o/2;                                                  % Distance between array elements
phi         = 2*pi/M;                                                      % Separation between elements
radius      = D/phi;                                                       % Radius of UCA
alpha       = 2;                                                           % Decay factor
Ptx         = 10;                                                          % BS power
G_BS        = 2.15;                                                        % BS antenna gain [dB]
sigma_c     = 1;                                                           % Power of BS-to-sensor channel (fading)

% 3. Sensors
K           = 1e1;                                                         % Number of clusters
R           = 4e1;                                                         % Radius of clusters
C           = zeros(K,2);                                                  % Center of clusters
S           = 75*R/K;                                                      % Region of clusters
d           = 2*R;                                                         % Minimum distance between clusters
d_min       = R;                                                           % Minimum distance to BS and UTs
area        = pi*R^2;                                                      % Cluster area
Lambda      = 0.05;                                                        % Mean density
lambda      = zeros(1,K);
m_dist      = 1;                                                           % Minimum distance between sensors 
p_act       = 0.1;                                                         % Probability factor
P_s         = 0.01;                                                        % Transmit power of sensors
sigma_g     = 1;                                                           % Power of sensor-to-sensor channel (fading)

% 4. UTs
UT_NF       = 7;                                                           % UT noise figure [dB]
Thrm_Noise  = -174;                                                        % Thermal noise level [dBm/Hz]
BW          = 94*180e3;                                                    % Bandwidth
thm_noise   = UT_NF + Thrm_Noise + 10*log10(BW);
thm_noise   = 10^(thm_noise*0.1)*1e-3;                                                                
sigma_w     = sqrt(thm_noise);                                             % Noise power at UT (normalized)
sigma_x     = 1;                                                           % Power of UT-to-sensor channel (fading)
G_UT        = 2.15;                                                        % UT antenna gain [dB]
AG_DL       = 10^(G_BS*0.1)*10^(G_UT*0.1);                                 % Total antenna gain (DL)
AG_UL       = 10^(G_UT*0.1)*10^(G_UT*0.1);                                 % Total gain UL (same at RX and TX: UT)

%% Magnitudes
% 1. Sensors 
% Generate disjoint clusters
C(1,:)              = (S - 2*R)*rand(1,2) + R*ones(1,2);
while norm(C(1,:) - [S S]/2) < d_min
    C(1,:)          = (S - 2*R)*rand(1,2) + R*ones(1,2);
end
k                   = 1;
counter             = 0;
while k < K
    c               = (S - 2*R)*rand(1,2) + R*ones(1,2);
    n               = zeros(k,1);
    for i = 1:k
       n(i)         = norm(c - C(i,:));
    end
    if sum(n >= d) == k && norm(c - [S S]/2) >= d_min
        k           = k + 1;
        C(k,:)      = c;
        counter     = 0;
    else
        counter     = counter + 1;
        if counter == 1000  
            C       = zeros(K,2);
            C(1,:)  = (S - 2*R)*rand(1,2) + R*ones(1,2);
            k       = 1; 
            counter = 0;
        end
    end
end

% Generate AoAs and distances between BS-sensors
aoa_s           = zeros(K,1);                                              % AoA w.r.t. x-axis
dist_s          = zeros(K,1);                                              % Distances to BS
for k = 1:K
    aoa_s(k)    = atan2(C(k,2) - S/2,C(k,1) - S/2); 
    dist_s(k)   = norm(C(k,:) - [S S]/2);
end

% Generate angles and distances between clusters
a_cl                = zeros(K,K);                                          % Argument
d_cl                = zeros(K,K);                                          % Modulus
theta               = zeros(K,K);                                          % Front angle entangled between the two tangent lines
for j = 1:K
    vector          = 1:K;
    vector          = vector(vector ~= j);
    for k = vector
        d_k         = norm(C(k,:) - C(j,:));
        a_k         = atan2(C(k,2) - C(j,2),C(k,1) - C(j,1));              % X-axis  
        
        d_cl(j,k)   = d_k;
        a_cl(j,k)   = a_k;
        theta(j,k)  = 2*asin(R/d_k);
    end
end

% 2. UTs
% Generate set of cellular UTs
N                   = 10;                                                  % Number of cellular UTs such that M > N - 1 
if M <= N - 1
   error('Too many UTs. BD Processing KO!'); 
end

U                   = zeros(N,2);
i                   = 1;
while i <= N
    U(i,:)          = S*rand(1,2);
    d_k             = zeros(1,K);
    d_BS            = norm(U(i,:) - [S S]/2);
    for k = 1:K
        d_k(k)      = norm(C(k,:) - U(i,:));
    end
    if sum(d_k < d_min) == 0 && d_BS <= S/2
        i           = i + 1;
    end
end

% AoA and distances to BS
aoa_u           = zeros(N,1);                                              % AoA w.r.t. x-axis
dist_u          = zeros(N,1);                                              % Distances to BS
for i = 1:N
    aoa_u(i)    = atan2(U(i,2) - S/2,U(i,1) - S/2); 
    dist_u(i)   = norm(U(i,:) - [S S]/2);
end

% AoA and distances to clusters
a_us                = zeros(N,K);                                          % Argument
d_us                = zeros(N,K);                                          % Modulus
for k = 1:K
    for i = 1:N
        d_k         = norm(C(k,:) - U(i,:));
        a_k         = atan2(C(k,2) - U(i,2),C(k,1) - U(i,1));              % X-axis          
        d_us(i,k)   = d_k;
        a_us(i,k)   = a_k;
    end
end

% Channels (BS to UTs)
P_z         = 10^2.4*1e-3;                                                 % Transmit power of UTs
sigma_h     = 1;                                                           % Power of BS-to-UT channel

G           = zeros(M,N);
for i = 1:N
    s       = exp(1i*K_o*radius*cos(aoa_u(i) - (0:M-1)*phi));    
    G(:,i)  = sqrt(AG_DL) * dist_u(i)^(-alpha/2) * sigma_h * s' / (lambda_o/(4*pi))^(-alpha/2);
end

% BD Processing
V_o             = zeros(M,M - N + 1,N);
for i = 1:N
    H           = G(:,[(1:i - 1) (i + 1:N)])';
    [~,~,V]     = svd(H);
    V_o(:,:,i)  = V(:,N:end);
end

%% Energy Statistics
% 1. Mean of harvested energy
% Sensors own cluster (eq. 14)
aux         = ComputeAverageDistance(R,m_dist,alpha);
aux         = aux(2);
E_main      = p_act*(1 - p_act)*sigma_g^2*P_s*aux;                         % Harvest power (UL)
E_main      = E_main*AG_UL;                                                % Add antenna gains
E_main      = E_main/(lambda_o/(4*pi))^(-alpha);                           % Distance normalization factor

E_main      = T_u*E_main;                                                  % Power to energy

% Add Other Clusters (eq. 14)
E_other                     = zeros(K,K);
for j = 1:K
    for k = j:K
        if j ~= k
            aux             = ComputeRealIntegralSG(R,a_cl(j,k),d_cl(j,k),theta(j,k),alpha);
            E_k             = p_act*(1 - p_act)*sigma_g^2*P_s*aux;                    
            E_k             = E_k*AG_UL;
            E_k             = E_k/(lambda_o/(4*pi))^(-alpha);
            E_other(j,k)    = E_k;
            E_other(k,j)    = E_other(j,k);
        end
    end
end
E_other                     = T_u*E_other;

% UTs (eq. 15 and 17)
E_UTs_o                 = zeros(N,K);
for k = 1:K
    for i = 1:N
        angle           = 2*asin(R/d_us(i,k));
        aux             = ComputeRealIntegralSG(R,a_us(i,k),d_us(i,k),angle,alpha);
        E_UTs_o(i,k)    = sigma_x^2*P_z*(1 - p_act)/area*aux;
    end
end
E_UTs_o                 = E_UTs_o*AG_UL;
E_UTs_o                 = E_UTs_o/(lambda_o/(4*pi))^(-alpha);

E_UTs_FD                = (T_f - p_act*T_u)*E_UTs_o/(1 - p_act);           % Power to energy
E_UTs_HD                = T_u*E_UTs_o;

% 2. Variance of harvested energy
% Sensors (eq. 22)
aux                         = ComputeAverageDistance(R,m_dist,2*alpha);
aux                         = aux(2);
V_main                      = p_act*(1 - p_act)*sigma_g^4*P_s^2*aux;       % Harvest power (UL)
b_1                         = 9 + 3*(N_s - 1)*(1 - p_act)*p_act;           % Include cross-products
V_main                      = V_main*AG_UL^2;                              % Add antenna gains
V_main                      = V_main/(lambda_o/(4*pi))^(-2*alpha);         % Distance normalization factor

V_cross_m                   = 0;                                                                
b_2                         = p_act + (N_s - 1)*(1 - p_act)*p_act;         % Include cross-products
V_main                      = b_1*V_main + b_2*V_cross_m;
V_main                      = T_u*T_s*V_main;

V_other                     = zeros(K,K);
for j = 1:K
    for k = j:K
        if j ~= k
            aux             = ComputeRealIntegralSG(R,a_cl(j,k),d_cl(j,k),theta(j,k),2*alpha);
            V_k             = p_act*(1 - p_act)*sigma_g^4*P_s^2*aux;       
            V_k             = V_k*AG_UL^2;
            V_k             = V_k/(lambda_o/(4*pi))^(-2*alpha);
            V_other(j,k)    = V_k;
            V_other(k,j)    = V_other(j,k);
        end
    end
end
V_cross_o                   = 0;                                                            
V_other                     = b_1*V_other + b_2*V_cross_o;
V_other                     = T_u*T_s*V_other;

% UTs (eq. 23)
c_1_HD                          = 9*T_u*T_s*(1 - p_act)*(1 + (N_s - 1)*(1 - p_act));
c_2_HD                          = T_u*T_s*(1 - p_act) + T_u*T_s*(N_s - 1)*(1 - p_act)^2;
c_1_FD                          = c_1_HD + 9*T_d^2 + 18*T_d*T_u*(1 - p_act);
c_2_FD                          = c_2_HD + T_d^2 + 2*T_d*T_u*(1 - p_act);

V_UTs_HD                        = zeros(1,K);
V_UTs_FD                        = zeros(1,K);
for k = 1:K
   Var_main                     = zeros(1,N);
   for i = 1:N
       angle                    = 2*asin(R/d_us(i,k));
       aux                      = ComputeRealIntegralSG(R,a_us(i,k),d_us(i,k),angle,2*alpha);
       Var_main(i)              = 1/area*aux;
   end
   Var_main                     = sum(Var_main);

   Var_cross                    = zeros(N,N);
   for j = 1:N
       angle                    = 2*asin(R/d_us(j,k));
       aux2                     = ComputeRealIntegralSG(R,a_us(j,k),d_us(j,k),angle,alpha);
       for i = 1:N
           if j ~= i 
               angle            = 2*asin(R/d_us(i,k));
               aux1             = ComputeRealIntegralSG(R,a_us(i,k),d_us(i,k),angle,alpha);
               Var_cross(i,j)   = 1/area*aux1 * 1/area*aux2;
           end
       end 
   end
   V_UTs_HD(k)                  = c_1_HD*Var_main + c_2_HD*sum(sum(Var_cross));
   V_UTs_FD(k)                  = c_1_FD*Var_main + c_2_FD*sum(sum(Var_cross));
   
end
V_UTs_HD                        = P_z^2*sigma_x^4*V_UTs_HD;
V_UTs_HD                        = V_UTs_HD*AG_UL^2/(lambda_o/(4*pi))^(-2*alpha);
V_UTs_FD                        = P_z^2*sigma_x^4*V_UTs_FD;
V_UTs_FD                        = V_UTs_FD*AG_UL^2/(lambda_o/(4*pi))^(-2*alpha);

V_UTs_HD                        = V_UTs_HD - sum(E_UTs_HD).^2;                                
V_UTs_FD                        = V_UTs_FD - sum(E_UTs_FD).^2;

if duplex == 1
    E_UTs   = E_UTs_HD;
    V_UTs   = V_UTs_HD;
elseif duplex == 2
    E_UTs   = E_UTs_FD;
    V_UTs   = V_UTs_FD;
end

% Covariance matrix of BS-to-sensors channels (eq. 9)
CM_t                = zeros(M,M,K);
for k = 1:K
    CM              = zeros(M,M);
    angle           = 2*asin(R/dist_s(k));
    for j = 1:M
        for i = j:M
            s       = @(psi) exp(1i*K_o*radius*cos(psi - (i - 1)*phi))*exp(-1i*K_o*radius*cos(psi - (j - 1)*phi));
            aux     = ComputeRealIntegralSG_MISO(R,aoa_s(k),dist_s(k),angle,alpha,s);
            CM(i,j) = sigma_c^2/area*aux;
            CM(j,i) = conj(CM(i,j));
        end
    end
    CM              = CM*AG_DL;
    CM              = CM/(lambda_o/(4*pi))^(-alpha);
    CM_t(:,:,k)     = CM;
end

% Variance of BS-to-sensors channels (eq. 19 and 20)
VM          = zeros(K,1);                      
for k = 1:K
    aux     = ComputeRealIntegralSG(R,aoa_s(k),dist_s(k),angle,2*alpha);
    VM(k)   = 1/area*M^2*aux;
end
VM          = VM*AG_DL^2/(lambda_o/(4*pi))^(-2*alpha);

%% Optimization 
% 1. Energy constraints
Y       = 8; 
delta_o = linspace(1e-7,8e-7,Y);
T       = 10;
mu_o    = linspace(1e-9,1e-7,T);

% 2. Densities
if densities == 1
    lambda  = linspace(0.5*Lambda,2*Lambda,K);
elseif densities == 2
    lambda  = linspace(0.1*Lambda,Lambda,K);
end
lambda      = lambda(randperm(K));
     
% 3. Final mean and variance
E_main      = (area*lambda - 1)*E_main;
E_other     = diag(lambda)*E_other;

V_main      = (area*lambda - 1)*V_main;
V_other     = diag(lambda)*V_other;
V_sensors   = V_main + sum(V_other);
    
% 4. Mean and variance of energy threshold
aux_E       = E_main + sum(E_other) + sum(E_UTs);
aux_V       = V_sensors + V_UTs;

%% Sweeps (over delta and mu)
W           = T;
Q_total     = cell(Y,W);
W_total     = cell(Y,W);
R_sum       = zeros(Y,W);
E_total     = zeros(Y,W,K);
V_total     = zeros(Y,W,K);
V_real      = zeros(Y,W,4);
Percentage  = zeros(Y,W,Y);
delta_total = zeros(Y,W);
mu_total    = zeros(Y,W);

for y = 1:Y    
    disp(['Delta ' num2str(delta_o(y)) ' (' num2str(y) ' of '  num2str(Y) ')']);
    for w = W:-1:1
        disp(['Mu ' num2str(mu_o(w)) ' (' num2str(w) ' of '  num2str(W) ')']); newline;
        
        % Thresholds w.r.t. sensors' statistics
        delta   = delta_o(y) - aux_E;
        delta   = max(delta,0);        
        mu      = mu_o(w)/3 - aux_V;
        if sum(mu < 0) > 0
            break;
        end        
        disp(['Thresholds: Delta = ' num2str(max(delta)) ' (' num2str(delta_o(y)) ') and Mu = ' num2str(min(mu)) ' (' num2str(mu_o(y)) ')']);

        % DPC        
        if duplex == 1
            [Q_opt,R_opt,E_opt,V_opt]   = SolveDCP_HD(K,M,N,T_d,sigma_c,sigma_w,Ptx,VM,delta,mu,G,CM_t,V_o);
        else
            [QW,R_opt,E_opt,V_opt]      = SolveDCP_FD(K,M,N,T_d,T_u,N_s,p_act,sigma_c,sigma_w,Ptx,VM,{delta delta_o(y)},mu,G,CM_t,V_o);
            Q_opt                       = QW{1};
            W_opt                       = QW{2};
        end
        
        % Convert to bps
        R_opt = R_opt/log(2);        
        if isnan(R_opt)
            break;
        end        
        E_opt = E_opt + aux_E;
        V_opt = 3*(V_opt + aux_V);
               
        % Save results
        if duplex == 1
           W_total{y,w}     = W_opt; 
        end
        Q_total{y,w}        = Q_opt;
        R_sum(y,w)          = R_opt;                   
        E_total(y,w,:)      = E_opt;
        V_total(y,w,:)      = V_opt;        
        delta_total(y,w)    = delta_o(y);
        mu_total(y,w)       = mu_o(w);        

    end
end