function [ Q_opt, R_sum_opt, E_opt, V_opt ] = SolveDCP_HD( K, M, N, T_d, sigma_chi, sigma_w, P_tx, E_d, delta, mu, G, C, V_o )

rng(1);

%% Initialization
% Feasible point
Q_o             = zeros(M - N + 1,M - N + 1,N);
for i = 1:N 
   Q_o(:,:,i)   = eye(M - N + 1); 
end

V_sum           = zeros(M,M);
for i           = 1:N
    V_sum       = V_sum + V_o(:,:,i)*V_o(:,:,i)';
end
alpha           = P_tx/trace(V_sum);
Q_o             = alpha * Q_o;

% Sum rate, mean energy, and energy variance
R           = zeros(N,1);
for i = 1:N
    R(i)    = log(1 + 1/sigma_w^2 * real(G(:,i)'*V_o(:,:,i)*Q_o(:,:,i)*V_o(:,:,i)'*G(:,i)));
end
R_sum_new   = sum(R);

Q_sum       = zeros(M,M);
for i = 1:N
    Q_sum   = Q_sum + V_o(:,:,i)*Q_o(:,:,i)*V_o(:,:,i)';
end

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
R_sum_old           = 0;
m                   = 0;
Q_opt               = zeros(M,M,N);
while R_sum_new - R_sum_old > 1e-3 && m < 1e3
    disp(['Step ' num2str(m + 1)]);
    R_sum_old       = R_sum_new;

    tic    
    cvx_begin
    cvx_quiet(true)
        cvx_precision high
        cvx_solver sedumi
        variable Q(M - N + 1,M - N + 1,N) hermitian semidefinite
        R_sum       = 0;
        for i = 1:N
            R_sum   = R_sum + log(1 + 1/sigma_w^2 * real(G(:,i)'*V_o(:,:,i)*Q(:,:,i)*V_o(:,:,i)'*G(:,i)));
        end
        maximize R_sum
        subject to
            % Constraint 1
            power       = 0;
            for i = 1:N
                power   = power + real(trace(Q(:,:,i)));
            end
            power       = scale_P*power;
            power       <= scale_P*P_tx;
            
            % Constraints 2 & 3
            Q_sum           = zeros(M,M);
            Q_sum_o         = zeros(M,M);
            fro_nrm         = 0;
            for i = 1:N
                Q_sum       = Q_sum + V_o(:,:,i)*Q(:,:,i)*V_o(:,:,i)';
                Q_sum_o     = Q_sum_o + V_o(:,:,i)*Q_o(:,:,i)*V_o(:,:,i)';
            end
            for n = 1:M
                for l = 1:M
                    fro_nrm = fro_nrm + square_abs(Q_sum(l,n));
                end 
            end    
            for k = 1:K
                % Mean
                E(k)    = T_d*real(vec(C(:,:,k))'*vec(Q_sum));
                E(k)    = scale_E*E(k);
                E(k)    >= scale_E*delta(k);
                % Variance
                F       = real(T_d^2*(vec(C(:,:,k))'*vec(Q_sum_o))^2 + 2*T_d^2*vec(C(:,:,k))'*vec(Q_sum_o)*vec(C(:,:,k)).'*vec(Q_sum - Q_sum_o));               
                V(k)    = 3*T_d^2*sigma_chi^4*E_d(k)*fro_nrm - F; 
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
       break;
    else
       clear E V
    end

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

E_opt               = zeros(1,K);
V_opt               = zeros(1,K);
for k = 1:K
    E_opt(k)        = T_d*real(trace(C(:,:,k)*Q_sum));
    V_opt(k)        = 3*T_d^2*sigma_chi^4*E_d(k)*norm(Q_sum,'fro')^2 - E_opt(k)^2;
end

if sum(E_opt - delta < -1e-3) > 0 || sum(mu - V_opt < -1e-3) > 0
    error('Constraints on E or V not satisfied!');
end

end