function [Y_g,Y_l,Y_h,w,r,P_g,P_l,P_h,pi_g,pi_l,pi_h]=prodn_alloc_iter_transition(...
    K,rho,w,T_g,T_l,T_h,d_g,d_l,d_h,L,I,theta,eta,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,Pguess_g,Pguess_l,Pguess_h,...
    mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,...
    fi,omega_bar,omega_cg,omega_cl,omega_ch,omega_xg,omega_xl,omega_xh)
%--------------------------------------------------------------------------
% Given prices, trade shares, consumption, and investment this function
% computes production of intermediates at a point in time in a country.
%--------------------------------------------------------------------------

PY_g = ones(I,1);
PY_l = ones(I,1);
PY_h = ones(I,1);

r = (alpha_g.*nu_g.*PY_g+alpha_l.*nu_l.*PY_l+alpha_h.*nu_h.*PY_h)./K;
PX = rho.*(w.*L+r.*K);
% WPC = sum(w.*L + r.*K - PX); % sum of net transfer is zero
% PC = WPC*omega_bar.*fi; % net transfer is in terms of WPC
%transfer = sum( fi.*(w.*L + r.*K) ).*L./sum(L);
transfer = sum( fi.*(w.*L + r.*K) ).*(w.*L + r.*K)./sum(w.*L + r.*K);
%transfer = sum( fi.*(w.*L + r.*K) )/I;

PC = (1-fi).*(w.*L + r.*K)+transfer-PX;


% compute sector price
[P_g,P_l,P_h]=priceindex_planner(w,r,...
    T_g,T_l,T_h,d_g,d_l,d_h,I,theta,eta,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,Pguess_g,Pguess_l,Pguess_h,...
    mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh);

% compute trade share
[pi_g,pi_l,pi_h] = tradeshare_planner(w,r,...
    P_g,P_l,P_h,T_g,T_l,T_h,d_g,d_l,d_h,I,theta,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,...
    mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh);


dif=10;
iter=0;
while dif>1e-6
    iter=iter+1;
    
    Ups_g =  kron( ones(I,1) , mu_gg'.*(1-nu_g') ) .* pi_g'; 
    V_g = pi_g' * ( mu_lg.*(1-nu_l).*PY_l + mu_hg.*(1-nu_h).*PY_h ); 
    F_g = pi_g' * (omega_cg.*PC+omega_xg.*PX); % final consumption demand
    PY_gnew = (eye(I)-Ups_g)\(V_g+F_g);
    
    Ups_l =  kron( ones(I,1) , mu_ll'.*(1-nu_l') ) .* pi_l';
    V_l = pi_l' * ( mu_gl.*(1-nu_g).*PY_g + mu_hl.*(1-nu_h).*PY_h );
    F_l = pi_l' * (omega_cl.*PC+omega_xl.*PX); % final consumption demand
    PY_lnew = (eye(I)-Ups_l)\(V_l+F_l);
    
    Ups_h =  kron( ones(I,1) , mu_hh'.*(1-nu_h') ) .* pi_h';
    V_h = pi_h' * ( mu_gh.*(1-nu_g).*PY_g + mu_lh.*(1-nu_l).*PY_l );
    F_h = pi_h' * (omega_ch.*PC+omega_xh.*PX); % final consumption demand
    PY_hnew = (eye(I)-Ups_h)\(V_h+F_h);
    
         
    dif = max(abs([PY_gnew; PY_lnew; PY_hnew]-[PY_g; PY_l; PY_h])./[PY_g; PY_l; PY_h]);
    smooth=.1*rand+.9;
    PY_g=smooth*PY_gnew+(1-smooth)*PY_g;
    PY_l=smooth*PY_lnew+(1-smooth)*PY_l;
    PY_h=smooth*PY_hnew+(1-smooth)*PY_h;
    
    % update prices
    r = (alpha_g.*nu_g.*PY_g+alpha_l.*nu_l.*PY_l+alpha_h.*nu_h.*PY_h)./K;
    PX = rho.*(w.*L+r.*K);
%    WPC = sum(w.*L + r.*K - PX); 
%     PC = WPC*omega_bar.*fi; 
    %transfer = sum( fi.*(w.*L + r.*K) ).*L./sum(L);
    transfer = sum( fi.*(w.*L + r.*K) ).*(w.*L + r.*K)./sum(w.*L + r.*K);
    %transfer = sum( fi.*(w.*L + r.*K) )/I;
    PC = (1-fi).*(w.*L + r.*K)+transfer-PX;
    
    % compute sector price
    [P_g,P_l,P_h]=priceindex_planner(w,r,...
        T_g,T_l,T_h,d_g,d_l,d_h,I,theta,eta,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,Pguess_g,Pguess_l,Pguess_h,...
        mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh);
    
    % compute trade share
    [pi_g,pi_l,pi_h] = tradeshare_planner(w,r,...
        P_g,P_l,P_h,T_g,T_l,T_h,d_g,d_l,d_h,I,theta,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,...
        mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh);
        
end

% sector output
Y_g = PY_g./P_g;
Y_l = PY_l./P_l;
Y_h = PY_h./P_h;

% denominate in WGDP
WGDP = sum(w.*L + r.*K); 
numeraire = WGDP;
w = w/numeraire;
r = r/numeraire;
P_g = P_g/numeraire;
P_l = P_l/numeraire;
P_h = P_h/numeraire;



