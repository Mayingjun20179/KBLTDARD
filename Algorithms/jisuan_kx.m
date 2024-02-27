%%%%计算变分参数ξ(kx)
function kx = jisuan_kx(U,V)

% %%%%测试
% M = 500;
% N = 600;
% R = 5;
% U.mu = randn(M,R);
% U.sigma = zeros(R,R,M);
% for m = 1:M
%     B = randn(100,R);
%     U.sigma(:,:,m) = B'*B+0.01*eye(R);
% end
% U.mu1 = khatrirao_fast(U.mu',U.mu'); %
% U.sigma1 = reshape(U.sigma,R*R,M);   %表示将U.sigma展开
% 
% V.mu = randn(N,R);
% V.sigma = zeros(R,R,N);
% for n = 1:N
%     B = randn(100,R);
%     V.sigma(:,:,n) = B'*B+0.01*eye(R);
% end
% V.mu1 = khatrirao_fast(V.mu',V.mu'); %
% V.sigma1 = reshape(V.sigma,R*R,N);   %表示将U.sigma展开
% 
% 
% tic
% U.mu1 = khatrirao_fast(U.mu',U.mu'); %
% U.sigma1 = reshape(U.sigma,R*R,M);   %表示将U.sigma展开
% V.mu1 = khatrirao_fast(V.mu',V.mu'); %
% V.sigma1 = reshape(V.sigma,R*R,N);   %表示将U.sigma展开
% kx2 = U.mu1'* V.mu1 + U.sigma1'* V.mu1+...
%     U.mu1'* V.sigma1 +  U.sigma1'* V.sigma1;
% toc
% length(find(kx2<0))
% 
% tic
% M = size(U.mu,1);
% N = size(V.mu,1);
% kx = zeros(M,N);
% for m = 1:M
%     for n = 1:N
%       kx(m,n) =  (U.mu(m,:) * V.mu(n,:)')^2+...
%           sum(sum(U.sigma(:,:,m).*(V.mu(n,:)'*V.mu(n,:))))+...
%           sum(sum(V.sigma(:,:,n).*(U.mu(m,:)'*U.mu(m,:))))+...
%           sum(sum(U.sigma(:,:,m).*V.sigma(:,:,n)));
%     end
% end
% toc
% 
% length(find(kx<0))
% max(max(abs(kx-kx2)))

%%%%判断协方差矩阵是否为半正定

M = size(U.mu,1);
N = size(V.mu,1);
kx = zeros(M,N);
for m = 1:M
    [~,value1] = eig(U.sigma(:,:,m)) ;
    if length(find(diag(value1)<0))
        disp('m')
        error('出错')
    end
end



% for m = 1:M
%     for n = 1:N
%       kx =  (U.mu(m,:) * V.mu(n,:)')^2;
%       kx1 =  sum(sum((U.mu(m,:)'*U.mu(m,:)).*(V.mu(n,:)'*V.mu(n,:))));
%       if abs(kx-kx1)>10^-10
%           disp([m,n])
%           error('错误')
%       end
%     end
% end






end