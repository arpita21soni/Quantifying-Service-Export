function [P_g,P_l,P_h]=priceindex_planner(w,r,...
    T_g,T_l,T_h,d_g,d_l,d_h,I,theta,eta,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,Pguess_g,Pguess_l,Pguess_h,...
    mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh)

%--------------------------------------------------------------------------
% Computes a relevant price index for each country.
%--------------------------------------------------------------------------

gam=gamma(1+1/theta.*(1-eta)).^(1/(1-eta));

P_g = Pguess_g;
P_l = Pguess_l;
P_h = Pguess_h;

dif=10;
iter=0;
while dif>1e-6
    iter=iter+1;
    
    u_g = (r./(alpha_g.*nu_g)).^(alpha_g.*nu_g) ...
        .* (w./((1-alpha_g).*nu_g)).^((1-alpha_g).*nu_g) ...
        .* (P_g./( mu_gg.*(1-nu_g) )).^( mu_gg.*(1-nu_g) )...
        .* (P_l./( mu_gl.*(1-nu_g) )).^( mu_gl.*(1-nu_g) )...
        .* (P_h./( mu_gh.*(1-nu_g) )).^( mu_gh.*(1-nu_g) ); % dimension 3*1
    part1_g = u_g.^(-theta).*T_g;
    together_g = d_g.^(-theta).*repmat(part1_g',I,1); % (d_a 3*3) .* (3*3)
    P_gnew = gam.*sum(together_g,2).^(-1/theta); % sum over column for each row, (1,1),(1,2), destination first, source second
    
    u_l = (r./(alpha_l.*nu_l)).^(alpha_l.*nu_l) ...
        .* (w./((1-alpha_l).*nu_l)).^((1-alpha_l).*nu_l) ...
        .* (P_g./( mu_lg.*(1-nu_l) )).^( mu_lg.*(1-nu_l) )...
        .* (P_l./( mu_ll.*(1-nu_l) )).^( mu_ll.*(1-nu_l) )...
        .* (P_h./( mu_lh.*(1-nu_l) )).^( mu_lh.*(1-nu_l) );
    part1_l = u_l.^(-theta).*T_l;
    together_l = d_l.^(-theta).*repmat(part1_l',I,1);
    P_lnew = gam.*sum(together_l,2).^(-1/theta);
    
    u_h = (r./(alpha_h.*nu_h)).^(alpha_h.*nu_h) ...
        .* (w./((1-alpha_h).*nu_h)).^((1-alpha_h).*nu_h) ...
        .* (P_g./( mu_hg.*(1-nu_h) )).^( mu_hg.*(1-nu_h) )...
        .* (P_l./( mu_hl.*(1-nu_h) )).^( mu_hl.*(1-nu_h) )...
        .* (P_h./( mu_hh.*(1-nu_h) )).^( mu_hh.*(1-nu_h) );
    part1_h = u_h.^(-theta).*T_h;
    together_h = d_h.^(-theta).*repmat(part1_h',I,1);
    P_hnew = gam.*sum(together_h,2).^(-1/theta);
    
    dif = max(abs([P_gnew; P_lnew; P_hnew]-[P_g; P_l; P_h])./[P_g; P_l; P_h]);
    smooth=.1*rand+.9;
    P_g=smooth*P_gnew+(1-smooth)*P_g;
    P_l=smooth*P_lnew+(1-smooth)*P_l;
    P_h=smooth*P_hnew+(1-smooth)*P_h;
end





