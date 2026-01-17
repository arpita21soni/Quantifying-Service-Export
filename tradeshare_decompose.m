function [pi_BM,pi_T,pi_u,pi_tao] = tradeshare_decompose(I,id,theta,T,u,d)

%--------------------------------------------------------------------------
% Decompose the trade share into T, u, tao
%--------------------------------------------------------------------------

% benchmark case
pi_BM = repmat(T(id,:),I,1).*( repmat(u(id,:),I,1).*squeeze(d(:,id,:)) ).^(-theta);
sum = 0;
for i=1:I
    sum = sum+repmat(T(i,:),I,1).*( repmat(u(i,:),I,1).*squeeze(d(:,i,:)) ).^(-theta);
end    
pi_BM = pi_BM./sum;

% constant T for treated country
T_temp = T;
T_temp(id,:) = T(id,1);

pi_T = repmat(T_temp(id,:),I,1).*( repmat(u(id,:),I,1).*squeeze(d(:,id,:)) ).^(-theta);
sum = 0;
for i=1:I
    sum = sum+repmat(T_temp(i,:),I,1).*( repmat(u(i,:),I,1).*squeeze(d(:,i,:)) ).^(-theta);
end    
pi_T = pi_T./sum;

% constant u for treated country
u_temp = u;
u_temp(id,:) = u(id,1);

pi_u = repmat(T(id,:),I,1).*( repmat(u_temp(id,:),I,1).*squeeze(d(:,id,:)) ).^(-theta);
sum = 0;
for i=1:I
    sum = sum+repmat(T(i,:),I,1).*( repmat(u_temp(i,:),I,1).*squeeze(d(:,i,:)) ).^(-theta);
end    
pi_u = pi_u./sum;

% constant tao for treated country
d_temp = d;
for i=1:I
    d_temp(i,id,:) = d(i,id,1);
end

pi_tao = repmat(T(id,:),I,1).*( repmat(u(id,:),I,1).*squeeze(d_temp(:,id,:)) ).^(-theta);
sum = 0;
for i=1:I
    sum = sum+repmat(T(i,:),I,1).*( repmat(u(i,:),I,1).*squeeze(d_temp(:,i,:)) ).^(-theta);
end    
pi_tao = pi_tao./sum;





