%%%%
function P = KBLTD_ARD(Y,Su,Sv, parameter)


[I,J,K] = size(Y);

Su = double(Su);
Sv = double(Sv);

R   = parameter.maxRank;


%% 初始化
Y = tensor(Y);

a_lambda0 = 1;  %lata
b_lambda0 = 1;  %lata

lambdas = ones(R,1);  %

lambdag_indices = repmat(logical(eye(R)), [1, 1, I]);  
lambdah_indices = repmat(logical(eye(R)), [1, 1, J]);

%
randn('state', 1);
U.mu = randn(I,R);
U.sigma = repmat(eye(I), [1, 1, R]);     %

randn('state', 1);
V.mu = randn(J,R);
V.sigma = repmat(eye(J), [1, 1, R]);     %


%
% G.mu = Su*U.mu;
randn('state', 1);
G.mu = randn(I,R);
G.sigma = repmat(eye(R), [1 1 I]);     %

% H.mu = Sv*V.mu;
randn('state', 1);
H.mu = randn(J,R);
H.sigma = repmat(eye(R), [1 1 J]);     %


randn('state', 1);
W.mu = randn(K,R);
W.sigma = repmat(eye(R), [1 1 K]);     %

% --------- E(aa') = cov(a,a) + E(a)E(a')----------------
%
Bg = (reshape(G.sigma, [R*R, I]))';
Bh = (reshape(H.sigma, [R*R, J]))';
Bw = (reshape(W.sigma, [R*R, K]))';


SuSu = Su * Su';
SvSv = Sv * Sv';


%%%
kx = ones(I,J,K);  %
c = 5;
Aijk = (c*Y-1+Y)/2;
lam_kx = jisuan_lamb(kx);
Bijk = (c*Y+1-Y).*lam_kx;


%% Model learning
for it = 1:50
    disp(it)
    
    %%%%%%%%%%%%%%%%%
    sigma_g_a = (I*R)/2;
    sigma_g_b = trace((G.mu-Su*U.mu)*(G.mu-Su*U.mu)')+...
        sum(G.sigma(lambdag_indices));
    for r = 1:R
        sigma_g_b = sigma_g_b+trace(Su*U.sigma(:,:,r)*Su');
    end
    sigma_g_b = sigma_g_b/2;    
    sigma_g = sigma_g_a/sigma_g_b;
    
    
    sigma_h_a = (J*R)/2;
    sigma_h_b = trace((H.mu-Sv*V.mu)*(H.mu-Sv*V.mu)')+...
        sum(H.sigma(lambdah_indices));
    for r = 1:R
        sigma_h_b = sigma_h_b+trace(Sv*V.sigma(:,:,r)*Sv');
    end
    sigma_h_b = sigma_h_b/2;
    sigma_h = sigma_h_a/sigma_h_b;
      
%     [sigma_g,sigma_h]
    
    %% Update factor matrices(G,H,W)
    Aw = diag(lambdas);
    %
    ENZZT = reshape(khatrirao_fast({Bg,Bh},'r')' * double(tenmat(Bijk,3)'), [R,R,K]);
    FslashY = khatrirao_fast({G.mu,H.mu},'r')' * tenmat(Aijk, 3)';
    for k = 1:K
        W.sigma(:,:,k) = (2*ENZZT(:,:,k) + Aw )^(-1);  %
        %             ZSigma{n}(:,:,i) = (beta * ENZZT(:,:,i) + Aw )\eye(R);
        W.mu(k,:) = (W.sigma(:,:,k) * FslashY(:,k))';  %
    end
    Bw = (reshape(W.sigma, [R*R, K]) + khatrirao_fast(W.mu',W.mu'))';
   
    
    %G
    ENZZT = reshape(khatrirao_fast({Bh,Bw},'r')' * double(tenmat(Bijk,1)'), [R,R,I]);
    FslashY = khatrirao_fast({H.mu,W.mu},'r')' * tenmat(Aijk, 1)';
    for i=1:I
        G.sigma(:,:,i) = (2*ENZZT(:,:,i) + sigma_g*eye(R))^(-1);  %
        %             ZSigma{n}(:,:,i) = (beta * ENZZT(:,:,i) + Aw )\eye(R);
        G.mu(i,:) = (G.sigma(:,:,i) * (FslashY(:,i)+ sigma_g*U.mu'*Su(i,:)'))';  %
    end
    Bg = (reshape(G.sigma, [R*R, I]) + khatrirao_fast(G.mu',G.mu'))';
 
    %H
    ENZZT = reshape(khatrirao_fast({Bg,Bw},'r')' * double(tenmat(Bijk,2)'), [R,R,J]);
    FslashY = khatrirao_fast({G.mu,W.mu},'r')' * tenmat(Aijk, 2)';
    for j = 1:J
        H.sigma(:,:,j) = (2 * ENZZT(:,:,j) + sigma_h*eye(R))^(-1);  %表示后验协方差矩阵
        %             ZSigma{n}(:,:,i) = (beta * ENZZT(:,:,i) + Aw )\eye(R);
        H.mu(j,:) = (H.sigma(:,:,j) * (FslashY(:,j) + sigma_h*V.mu'*Sv(j,:)'))';  %这个表示后验期望
    end
    Bh = (reshape(H.sigma, [R*R, J]) + khatrirao_fast(H.mu',H.mu'))';
    
    
    
    %U
    for r = 1:R
        U.sigma(:, :, r) = (sigma_g*SuSu + lambdas(r) *eye(I))^(-1);
        U.mu(:, r) = sigma_g*U.sigma(:, :, r) * Su' * G.mu(:,r);
    end
    
    %V
    for r = 1:R
        V.sigma(:, :, r) = (sigma_h*SvSv + lambdas(r) *eye(J))^(-1);
        V.mu(:, r) = sigma_h*V.sigma(:, :, r) * Sv' * H.mu(:, r);
    end
    
    %% λ  ******************************
    a_gammaN = (0.5*(I+J+K) + a_lambda0)*ones(R,1);
    b_gammaN = zeros(R,1);    
    for r=1:R        
        b_gammaN(r) = V.mu(:,r)'*V.mu(:,r) + trace(V.sigma(:,:,r))+...
            U.mu(:,r)'* U.mu(:,r) + trace(U.sigma(:,:,r))+...
            W.mu(:,r)'* W.mu(:,r);
        for k=1:K
            b_gammaN(r) = b_gammaN(r) + W.sigma(r,r,k);  
        end
    end
    b_gammaN = b_lambda0 + 0.5.*b_gammaN;
    lambdas = a_gammaN ./ b_gammaN;   %gama
    %     lambdas(lambdas < 10^-6) = 10^-6;
   
%     lambdas'

    %%%% updateξ(kx)
    kx2 = reshape(khatrirao_fast({Bg,Bh,Bw},'r') * ones(R*R,1),I,J,K);
    if length(find(kx2<0))>0
        disp(it)
        error('错误')
    end
    kx = sqrt(kx2);
    
    [lam_kx,sig_kx] = jisuan_lamb(kx);
    Bijk = (c*Y+1-Y).*lam_kx;      
   
end


%% Prepare the results
P0 = ktensor({G.mu,H.mu,W.mu});
P0 = double(arrange(P0));
P = exp(P0)./(1+exp(P0));


end


