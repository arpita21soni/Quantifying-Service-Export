function [pi_g,pi_l,pi_h] = tradeshare_planner(w,r,P_g,P_l,P_h,T_g,T_l,T_h,d_g,d_l,d_h,I,theta,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh)

%--------------------------------------------------------------------------
% Computes trade shares pi_ij.
%--------------------------------------------------------------------------

u_g = (r./(alpha_g.*nu_g)).^(alpha_g.*nu_g) .* (w./((1-alpha_g).*nu_g)).^((1-alpha_g).*nu_g) ...
    .* (P_g./( mu_gg.*(1-nu_g) )).^( mu_gg.*(1-nu_g) )...
    .* (P_l./( mu_gl.*(1-nu_g) )).^( mu_gl.*(1-nu_g) )...
    .* (P_h./( mu_gh.*(1-nu_g) )).^( mu_gh.*(1-nu_g) ); % dimension 3*1
part1_g = u_g.^(-theta).*T_g; % 3*1
together_g = d_g.^(-theta).*repmat(part1_g',I,1);
denom_g = sum(together_g,2); % 3*1
pi_g = together_g.*kron(ones(1,I),1./denom_g); % kronecker product, sth like repmat, but with product

u_l = (r./(alpha_l.*nu_l)).^(alpha_l.*nu_l) .* (w./((1-alpha_l).*nu_l)).^((1-alpha_l).*nu_l) ...
    .* (P_g./( mu_lg.*(1-nu_l) )).^( mu_lg.*(1-nu_l) )...
    .* (P_l./( mu_ll.*(1-nu_l) )).^( mu_ll.*(1-nu_l) )...
    .* (P_h./( mu_lh.*(1-nu_l) )).^( mu_lh.*(1-nu_l) );
part1_l = u_l.^(-theta).*T_l;
together_l = d_l.^(-theta).*repmat(part1_l',I,1);
denom_l = sum(together_l,2); 
pi_l = together_l.*kron(ones(1,I),1./denom_l); 

u_h = (r./(alpha_h.*nu_h)).^(alpha_h.*nu_h) .* (w./((1-alpha_h).*nu_h)).^((1-alpha_h).*nu_h) ...
    .* (P_g./( mu_hg.*(1-nu_h) )).^( mu_hg.*(1-nu_h) )...
    .* (P_l./( mu_hl.*(1-nu_h) )).^( mu_hl.*(1-nu_h) )...
    .* (P_h./( mu_hh.*(1-nu_h) )).^( mu_hh.*(1-nu_h) );
part1_h = u_h.^(-theta).*T_h;
together_h = d_h.^(-theta).*repmat(part1_h',I,1);
denom_h = sum(together_h,2); 
pi_h = together_h.*kron(ones(1,I),1./denom_h); 



