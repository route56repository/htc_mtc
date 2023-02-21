function [ QW, R_sum_opt, E_opt, V_opt ] = SolveDCP_FD( K, M, N, T_d, T_u, N_s, p_act, sigma_chi, sigma_w, P_tx, E_d, delta, mu, G, C, V_o )

rng(1);

%% Parameters
% Thresholds
delta_o         = delta{2};
delta           = delta{1};

% Temporal/auxiliar
T_f             = T_u + T_d;
T_s             = T_u/N_s;
T_a             = T_f - p_act*T_u;
a_1             = 3*(T_d^2 + T_u*(1 - p_act)*(T_s*(1 + (N_s - 1)*(1 - p_act)) + 2*T_d));
a_2             = 3*T_d^2;
a_3             = 3*(T_d^2 + T_d*T_u*(1 - p_act));

%% Initialization
% Feasible points
Q_opt           = zeros(M,M,N);
Q_o             = zeros(M - N + 1,M - N + 1,N);
for i = 1:N 
   Q_o(:,:,i)   = eye(M - N + 1); 
end

V_sum           = zeros(M,M);
for i           = 1:N
    V_sum       = V_sum + V_o(:,:,i)*V_o(:,:,i)';
end
alpha           = 0.5 * P_tx/trace(V_sum);
Q_o             = alpha * Q_o;

W_o             = eye(M);           
alpha           = 0.5 * T_f * P_tx/(T_d * M);
W_o             = alpha * W_o;

% Sum rate, mean energy, and energy variance
R               = zeros(N,1);
for i = 1:N
    R(i)        = log(1 + 1/sigma_w^2 * real(G(:,i)'*V_o(:,:,i)*Q_o(:,:,i)*V_o(:,:,i)'*G(:,i)));
end
R_sum_new       = sum(R);

Q_sum           = zeros(M,M);
for i = 1:N
    Q_sum       = Q_sum + V_o(:,:,i)*Q_o(:,:,i)*V_o(:,:,i)';
end
P_new           = real(T_f*trace(Q_sum) + T_d*trace(W_o));

%% Iterative Optimization
% Constraint set (scaled)
order_P     = round(log10(P_tx));
order_E     = round(log10(max(delta(delta > 0))));
if isempty(order_E)
    order_E = 0;
end
order_V     = round(log10(min(mu)));
order_R     = round(log10(R_sum_new));

scale_P     = 10^order_R/10^order_P;
scale_E     = 10^order_R/10^order_E;
scale_V     = 10^order_R/10^order_V;

% Variables
R_sum_old       = 0;
P_old           = Inf;
counter         = 0;
while R_sum_new - R_sum_old > 1e-3 && counter < 1e3
% Q (TX)
m                       = 0;
W_opt                   = W_o;
while R_sum_new - R_sum_old > 1e-3 && m < 1e3
    if counter == 1 && m == 0
        R_sum_old       = 0;
    else
        R_sum_old       = R_sum_new;
    end
    
    tic
    cvx_begin
    cvx_quiet(true)
        cvx_precision high
        cvx_solver sedumi
        variable Q(M - N + 1,M - N + 1,N) hermitian semidefinite
        R_sum           = 0;
        for i = 1:N
            R_sum       = R_sum + log(1 + 1/sigma_w^2 * real(G(:,i)'*V_o(:,:,i)*Q(:,:,i)*V_o(:,:,i)'*G(:,i)));
        end
        maximize R_sum
        subject to
            % Constraint 1
            power       = 0;
            for i = 1:N
                power   = power + T_f*real(trace(Q(:,:,i)));
            end
            power       = power + T_d*real(trace(W_opt));
            power       = scale_P*power;
            power       <= scale_P*T_f*P_tx;
            % Constraint 2 & 3
            Q_sum       = zeros(M,M);
            Q_sum_o     = zeros(M,M);
            fro_nrm = 0;
            for i = 1:N
                Q_sum   = Q_sum + V_o(:,:,i)*Q(:,:,i)*V_o(:,:,i)';
                Q_sum_o = Q_sum_o + V_o(:,:,i)*Q_o(:,:,i)*V_o(:,:,i)';
            end
            for n = 1:M
                for l = 1:M
                    fro_nrm = fro_nrm + square_abs(Q_sum(l,n));
                end 
            end
            for k = 1:K
                % Mean
                E_tx    = T_a*real(vec(C(:,:,k))'*vec(Q_sum));
                E_wpt   = T_d*real(trace(C(:,:,k)*W_opt));
                E(k)    = E_tx + E_wpt;
                E(k)    = scale_E*E(k);
                E(k)    >= scale_E*delta(k);
                % Variance
                V_tot   = a_1*fro_nrm + a_2*real(norm(W_opt,'fro')^2) + 2*a_3*real(vec(W_opt)'*vec(Q_sum));
                F       = real(T_a^2*(vec(C(:,:,k))'*vec(Q_sum_o))^2 + 2*T_a^2*vec(C(:,:,k))'*vec(Q_sum_o)*vec(C(:,:,k))'*vec(Q_sum - Q_sum_o));
                F       = F + E_wpt^2 + 2*E_tx*E_wpt;
                V(k)    = sigma_chi^4*E_d(k)*V_tot;
                V(k)    = V(k) - F;
                V(k)    = scale_V*V(k);
                V(k)    <= scale_V*mu(k);            
            end
    cvx_end    
    time                = toc;
    status              = cvx_status;
      
    % Re-scale mean and variance
    E = E/scale_E;
    V = V/scale_V;
    
    % Update feasible point and sum rate    
    Q_o             = Q;
    R               = zeros(N,1);
    for i = 1:N
        R(i)        = log(1 + 1/sigma_w^2 * real(G(:,i)'*V_o(:,:,i)*Q(:,:,i)*V_o(:,:,i)'*G(:,i)));
    end
    R_sum_new       = sum(R);
    m               = m + 1;
    
    disp([status '. Old Sum Rate: ' num2str(R_sum_old) ' - New Sum Rate: ' num2str(R_sum_new) ' (' num2str(time) ' s)']);

    if (R_sum_old - R_sum_new)/R_sum_old > 1e-3 || sum(E - delta < -1e-3) > 0 || sum(mu - V < -1e-3) > 0
       R_sum_new    = NaN;
       Q            = NaN*ones(M - N + 1,M - N + 1,N);
       W            = NaN*ones(M,M);
       break;
    else
       clear E V
    end

end

if isnan(R_sum_new)
    R_sum_new       = NaN;
    Q               = NaN*ones(M - N + 1,M - N + 1,N);
    W               = NaN*ones(M,M);
    break;
else
    Q_sum_opt       = zeros(M,M);
    for i = 1:N
        Q_sum_opt   = Q_sum_opt + V_o(:,:,i)*Q(:,:,i)*V_o(:,:,i)';
    end
end

% W (WPT)
m               = 0;
while P_old - P_new > 1e-3 && m < 1e3
    if counter == 1 && m == 0
        P_old   = Inf;
    else
        P_old   = P_new;
    end

    tic    
    cvx_begin
    cvx_quiet(true)
        cvx_precision high
        cvx_solver sedumi
        variable P_total
        variable W(M,M) hermitian semidefinite
        minimize P_total
        subject to
            % Constraint 1            
            power   = T_f*real(trace(Q_sum_opt)) + T_d*real(trace(W));
            power   = scale_P*power;
            power   <= scale_P*T_f*P_total;
            % Constraints 2 & 3
            fro_nrm         = 0;
            for n = 1:M
                for l = 1:M
                    fro_nrm = fro_nrm + square_abs(W(l,n));
                end 
            end
            for k = 1:K
                % Mean
                E_tx    = T_a*real(vec(C(:,:,k))'*vec(Q_sum_opt));
                E_wpt   = T_d*real(vec(C(:,:,k))'*vec(W));
                E(k)    = E_tx + E_wpt;
                E(k)    = scale_E*E(k);
                E(k)    >= scale_E*delta(k);
                % Variance
                V_tot   = a_1*real(norm(Q_sum_opt,'fro')^2) + a_2*fro_nrm + 2*a_3*real(vec(Q_sum_opt)'*vec(W));                
                F       = real(T_d^2*(vec(C(:,:,k))'*vec(W_o))^2 + 2*T_d^2*vec(C(:,:,k))'*vec(W_o)*vec(C(:,:,k))'*vec(W - W_o));
                F       = F + E_tx^2 + 2*E_tx*E_wpt;
                V(k)    = sigma_chi^4*E_d(k)*V_tot;
                V(k)    = V(k) - F; 
                V(k)    = scale_V*V(k);
                V(k)    <= scale_V*mu(k);                
            end
    cvx_end    
    time                = toc;
    status              = cvx_status;
        
    % Re-scale mean and variance
    E = E/scale_E;
    V = V/scale_V;
    
    % Update initial feasible point and TX power    
    W_o             = W;
    P_new           = P_total;
    
    m               = m + 1;
    disp([status '. Old Power: ' num2str(P_old) ' - New Power: ' num2str(P_new) ' (' num2str(time) ' s)']);

    if sum(E - delta < -1e-3) > 0 || sum(mu - V < -1e-3) > 0
       R_sum_new    = NaN;
       Q            = NaN*ones(M - N + 1,M - N + 1,N);
       W            = NaN*ones(M,M);
       break;
    else
       clear E V
    end

end

counter             = counter + 1;

end

%% Save Optimal Values
for i = 1:N
    Q_opt(:,:,i)    = V_o(:,:,i)*Q(:,:,i)*V_o(:,:,i)';
end
R_sum_opt           = R_sum_new;

Q_sum               = zeros(M,M);
for i = 1:N
    Q_sum           = Q_sum + V_o(:,:,i)*Q(:,:,i)*V_o(:,:,i)';
end

W_opt               = W;

a_1                 = 3*T_d^2 + T_u*(1 - p_act)*(T_s*(3 + (N_s - 1)*(1 - p_act)) + 2*T_d);
a_2                 = 3*T_d^2;
a_3                 = 3*T_d^2 + T_d*T_u*(1 - p_act);

E_opt               = zeros(1,K);
V_opt               = zeros(1,K);
for k = 1:K
    E_tx            = (T_f - p_act*T_u)*real(trace(C(:,:,k)*Q_sum));
    E_wpt           = T_d*real(trace(C(:,:,k)*W));
    E_opt(k)        = E_wpt + E_tx;
    V_tx            = a_1*sigma_chi^4*E_d(k)*norm(Q_sum,'fro')^2;
    V_wpt           = a_2*sigma_chi^4*E_d(k)*norm(W,'fro')^2;
    V_cross         = 2*a_3*sigma_chi^4*E_d(k)*real(trace(W*Q_sum));
    V_opt(k)        = V_wpt + V_tx + V_cross;
    V_opt(k)        = V_opt(k) - E_opt(k)^2;
end

if sum(E_opt - delta < -1e-3) > 0 || sum(mu - V_opt < -1e-3) > 0
    error('Constraints on E or V not satisfied!');
end

QW{1}               = Q_opt;
QW{2}               = W_opt;

end