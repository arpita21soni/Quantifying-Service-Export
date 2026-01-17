function w = SS_planner(w0,A_x,T_g,T_l,T_h,d_g,d_l,d_h,L,I,...
    theta,eta,del,lambda,epsilon_c,epsilon_x,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,...
    beta,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,...
    fi2,omega_bar,omega_cg,omega_cl,omega_ch,omega_xg,omega_xl,omega_xh)

%------------------------------------------------------------------------
% Solves the world goods market clearing condition by iterating on wages.
% The output is the equilibrium matrix of wages, from which everything else
% can be backed out.
%------------------------------------------------------------------------

% Initial guess
w=w0;
dif=10; iter=0;
scl = 1e-1*10;

% display('---------------------------------------------------------------');
while dif>1e-6
    iter=iter+1;
    
    %w=w./sum(w.*L./(1-alpha),1);
    
    % Compute allocations as a funtion of wages
    [w,r,q,P_g,P_l,P_h,P_c,P_x,C,X,K,A,F,...
        K_g,K_l,K_h,L_g,L_l,L_h,Y_g,Y_l,Y_h,pi_g,pi_l,pi_h] = ...
        SS_alloc_planner(w,A_x,T_g,T_l,T_h,d_g,d_l,d_h,L,I,...
        theta,eta,del,lambda,epsilon_c,epsilon_x,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,...
        beta,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,...
        fi2,omega_bar,omega_cg,omega_cl,omega_ch,omega_xg,omega_xl,omega_xh);
    
    Borrow = -A; % net foreign asset, or net transfer
    Deficit = F; % trade deficit
        
    Z = Borrow-Deficit;
    Z = Z./(w.*L + r.*K); % check the difference of Z in terms of share of GDP
    Tw = w.*(1 + scl*Z);
    
    % Compute distance between 'new' and 'old' wage vectors.
    dif = max( abs( Z(:)) );
    
    % Update wage using smoothing parameter (convex combo b/w new and old).
    smooth = .2*rand + .8;
    w = smooth*Tw + (1-smooth)*w;
    
    if mod(iter,500)==0
        telapsed=toc; sec=mod(telapsed,60); min=floor(telapsed/60);
        fprintf('\t\t SS iterations on w completed: %6.0f\n',iter);
        fprintf('\t\t\t time elapsed: %6.0f min, %2.0f sec\n',min,sec);
        fprintf('\t\t\t dif: %3.10f\n',dif);
    end
end



