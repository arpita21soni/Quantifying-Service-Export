function [rho,w0,Zr] ...
    = comp_inv_labor_planner(rho0,w0,A_x,T_g,T_l,T_h,d_g,d_l,d_h,L,K_1,A_1,I,TT,...
    theta,eta,del,lambda,epsilon_c,epsilon_x,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,...
    beta,sigma,r2,P_c2,P_x2,C2,X2,K2,A_x2,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,...
    fi,omega_bar,omega_cg,omega_cl,omega_ch,omega_xg,omega_xl,omega_xh)

% Iterate over labor supply, within each loop, solve for equilibrium taking
% the labor supply as given.
dif=10;
iter=0;
rho=rho0; Trho = rho;
scl = 3e-1;
tol = 1e-2;
while dif>tol*0+1e-3
    iter=iter+1;
    
    smooth = 0.2.*rand(I,TT) + 0.8;
    rho = smooth.*Trho + (1-smooth).*rho;
    rho(rho<0.0001) = 0.0001;
    rho(rho>0.9999) = 0.9999;
    
    %w0=w0./repmat(sum(w0.*L./(1-alpha),1),I,1);
    
    % Compute equilibrium wages and interest rate along transition.
    w = compeqbm_planner(w0,A_x,T_g,T_l,T_h,d_g,d_l,d_h,L,K_1,A_1,I,TT,...
        theta,eta,del,lambda,epsilon_c,epsilon_x,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,...
        rho,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,...
        fi,omega_bar,omega_cg,omega_cl,omega_ch,omega_xg,omega_xl,omega_xh);
    w0=w;
    
    % Solve for remaining allocations as a function of wages and interest rate
    % along transition.
    [w,r,P_g,P_l,P_h,P_c,P_x,C,X,K,A,F, ...
        K_g,K_l,K_h,L_g,L_l,L_h,Y_g,Y_l,Y_h,pi_g,pi_l,pi_h] = ...
        alloc_planner(w,A_x,T_g,T_l,T_h,d_g,d_l,d_h,L,K_1,A_1,I,TT,...
        theta,eta,del,lambda,epsilon_c,epsilon_x,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,...
        rho,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,...
        fi,omega_bar,omega_cg,omega_cl,omega_ch,omega_xg,omega_xl,omega_xh);
        
        Zr = lambda.*beta.*( [r(:,2:TT),r2]+(1-lambda)./lambda.*[P_x(:,2:TT),P_x2].*[X(:,2:TT),X2]./[K(:,2:TT),K2]+...
        (1-del)./lambda./[A_x(:,2:TT),A_x2].*[P_x(:,2:TT),P_x2].*([X(:,2:TT),X2]./[K(:,2:TT),K2]).^(1-lambda) )./...
        ( P_x(:,1:TT)./A_x(:,1:TT).*( X(:,1:TT)./K(:,1:TT) ).^(1-lambda) )./...
        ( [P_c(:,2:TT),P_c2].*[C(:,2:TT),C2]./( P_c(:,1:TT).*C(:,1:TT) ) )-1;
    
    
    Zr(rho<=0.0001 & Zr<0)=-1e-6;
    Zr(rho>=0.9999 & Zr>0)=1e-6;
    
    dif = max(max( abs(Zr) ));
    Trho = rho.*(1 + scl.*Zr);
    
    if mod(iter,1)==0
        telapsed=toc; sec=mod(telapsed,60); mnt=floor(telapsed/60);
        fprintf('\t Iterations completed: %6.0f\n',iter);
        fprintf('\t\t time elapsed: %6.0f min, %2.0f sec\n',mnt,sec);
        fprintf('\t\t dif labor & inv rate: %3.10f\n',dif);
    end
end

% telapsed=toc; sec=mod(telapsed,60); mnt=floor(telapsed/60);
% display(sprintf('\t Iterations completed: %6.0f',iter));
% display(sprintf('\t\t time elapsed: %6.0f min, %2.0f sec',mnt,sec));
% display(sprintf('\t\t dif labor & inv rate: %3.10f',dif));




