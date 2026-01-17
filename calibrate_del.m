function [del,Ax_t,Ax_t_raw,Eq,A_x_guess] = calibrate_del(x0,I,beta,lambda,num_year,PC_t,PX_t,P_x_t,k_t,...
    alpha_g_t,alpha_l_t,alpha_h_t,PY_g_t,PY_l_t,PY_h_t,nu_g_t,nu_l_t,nu_h_t)

%--------------------------------------------------------------------------
% Computes counry specific del to match Ax
%--------------------------------------------------------------------------

%f = fsolve(@find_solution,x0);
%f = fminsearch(@find_solution,x0);
%f = lsqnonlin(@find_solution,x0);

% there is some issue for computing this, unless we use time-varying del

A = [];
b = [];
Aeq = [];
beq = [];
lb = ones(1,I)*0.01;
ub = ones(1,I)*0.99;
% lb = [0.01,0.15,0.05,0.03,0.08];
% ub = [0.03,0.20,0.07,0.05,0.10];
options = optimoptions('fmincon', 'Display', 'off');
f = fmincon(@find_solution, x0, A, b, Aeq, beq, lb, ub, [], options);

    function f = find_solution(x)
        
        del = x;
        
        X_t = PX_t./P_x_t;
        Ax_t_raw = ( k_t(:,2:17)-(1-del).*k_t(:,1:17-1) )./(X_t(:,1:17-1).^lambda.*k_t(:,1:17-1).^(1-lambda));
        A_x_guess = Ax_t_raw(:,1);        
        
        K_t = repmat(k_t(:,1),1,num_year+1); % input K2000    
        K_t(:,2) = (1-del).*K_t(:,1)+A_x_guess.*( X_t(:,1).^lambda.*K_t(:,1).^(1-lambda) );
        r_t(:,1) = ( alpha_g_t(:,1).*PY_g_t(:,1).*nu_g_t(:,1)+alpha_l_t(:,1).*PY_l_t(:,1).*nu_l_t(:,1)+...
            alpha_h_t(:,1).*PY_h_t(:,1).*nu_h_t(:,1) )./K_t(:,1);
        for i=1:num_year-1
            
            r_t(:,i+1) = ( alpha_g_t(:,i+1).*PY_g_t(:,i+1).*nu_g_t(:,i+1)+alpha_l_t(:,i+1).*PY_l_t(:,i+1).*nu_l_t(:,i+1)+...
                alpha_h_t(:,i+1).*PY_h_t(:,i+1).*nu_h_t(:,i+1))./K_t(:,i+1);
            
            LHS = PC_t(:,i+1)./PC_t(:,i).*PX_t(:,i)./( K_t(:,i+1)-(1-del).*K_t(:,i) );
            temp = LHS/beta./lambda-r_t(:,i+1)-(1-lambda)./lambda.*PX_t(:,i+1)./K_t(:,i+1);
            K_t(:,i+2) = (1-del)./lambda.*PX_t(:,i+1)./temp+(1-del).*K_t(:,i+1);
        end
        Ax_t = ( K_t(:,2:num_year+1)-(1-del).*K_t(:,1:num_year) )./(X_t.^lambda.*K_t(:,1:num_year).^(1-lambda));
                  
        %Eq = ( Ax_t(:,1:7)-Ax_t_raw )./Ax_t_raw;
        Eq = Ax_t(:,1:num_year -1)-Ax_t_raw;
        %Eq = Ax_t(:,1:7)./Ax_t(:,1)-Ax_t_raw./Ax_t_raw(:,1);
        f =  sum(Eq(:).^2);
    end

end



