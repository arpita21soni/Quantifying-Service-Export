function w = compeqbm_planner(w0,A_x,T_g,T_l,T_h,d_g,d_l,d_h,L,K_1,A_1,I,TT,...
    theta,eta,del,lambda,epsilon_c,epsilon_x,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,...
    rho,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,...
    fi,omega_bar,omega_cg,omega_cl,omega_ch,omega_xg,omega_xl,omega_xh)
%------------------------------------------------------------------------
% Solves the world goods market clearing condition by iterating on wages.
% The output is the equilibrium matrix of wages and the world interest rate
% from which everything else can be backed out.
%------------------------------------------------------------------------

% Initial guess
w=w0;
dif=1e5;
iter=0;
scl = 1e-1*10;

% Compute allocations as a funtion of wages

% -------------------------------------------------------------------
% Iterate on wage until deficit equals net borrowing.
while dif>1e-3 %& iter<300
    iter=iter+1;
    
    
    [w,r,P_g,P_l,P_h,P_c,P_x,C,X,K,A,F, ...
        K_g,K_l,K_h,L_g,L_l,L_h,Y_g,Y_l,Y_h,pi_g,pi_l,pi_h] = ...
        alloc_planner(w,A_x,T_g,T_l,T_h,d_g,d_l,d_h,L,K_1,A_1,I,TT,...
        theta,eta,del,lambda,epsilon_c,epsilon_x,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,...
        rho,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,...
        fi,omega_bar,omega_cg,omega_cl,omega_ch,omega_xg,omega_xl,omega_xh);
    
    Borrow = -A;
    Deficit = F;
    
    Z = Borrow-Deficit; % already normalized variables in WPC, no need to further normalize
    Z = Z./(w.*L + r.*K(:,1:TT)); % check the difference of Z in terms of share of GDP
    Tw = w.*(1 + scl*Z); 
    
    % Compute distance between 'new' and 'old' wage vectors.
    dif = max( abs( Z(:)) );
    
    % Update wage using smoothing parameter (convex combo b/w new and old).
    smooth = .2*rand + .8;
    w = smooth*Tw + (1-smooth)*w;
    
    
    if mod(iter,5)==0
        telapsed=toc; sec=mod(telapsed,60); min=floor(telapsed/60);
        fprintf('\t\t Iterations completed: %6.0f\n',iter);
        fprintf('\t\t\t time elapsed: %6.0f min, %2.0f sec\n',min,sec);
        fprintf('\t\t\t dif wage: %3.10f\n',dif);
    end
    
end




