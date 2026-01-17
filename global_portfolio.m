function [fi,diff] = global_portfolio(I,GDP,NX,L)

%--------------------------------------------------------------------------
% Computes share of income sent to global portfolio
%--------------------------------------------------------------------------

% However, directly compute the share will generate infeasible solution
LHS = zeros(I,I);
TL = sum(L);
small = 1e-5;

for i=1:I
    
   %LHS(i,:) = -L(i)/TL*GDP';
   LHS(i,:) = -GDP(i)/sum(GDP)*GDP';
   %LHS(i,:) = -1/I*GDP';
   LHS(i,i) = LHS(i,i)+GDP(i); 
end    
% fi = LHS\NX;

%fun = @(x)sum( ((LHS*x-NX)./NX).^2 );
%fun = @(x)sum( (LHS*x-NX).^2 );
fun = @(x)sum( ((LHS*x-NX)./GDP).^2 );
%fun = @(x)max( ((LHS*x-NX)./GDP).^2 );
x0 = 0.1*ones(I,1);
A = [];
b = [];
Aeq = [];
beq = [];
lb = ones(1,I)*small;
ub = ones(1,I)*(1-small);
%lb(5) = 0.05;
%lb(2) = 0.8; 
%ub(4) = 0.8;
options = optimoptions('fmincon', 'Display', 'off');
fi = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, [], options);





diff = LHS*fi-NX;



