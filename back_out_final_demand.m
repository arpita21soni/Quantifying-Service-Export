function [F_g,F_l,F_h,PQ_g,PQ_l,PQ_h]=back_out_final_demand(PY_g,PY_l,PY_h,pi_g,pi_l,pi_h,I,nu_g,nu_l,nu_h,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh)
%--------------------------------------------------------------------------
% Compute the model consistent sector final demand
%--------------------------------------------------------------------------

% sector final demand
Ups_g = kron( ones(I,1) , mu_gg'.*(1-nu_g') ) .* pi_g'; 
V_g = pi_g' * ( mu_lg.*(1-nu_l).*PY_l + mu_hg.*(1-nu_h).*PY_h );
F_g = pi_g'\((eye(I)-Ups_g)*PY_g-V_g);

Ups_l =  kron( ones(I,1) , mu_ll'.*(1-nu_l') ) .* pi_l';
V_l = pi_l' * ( mu_gl.*(1-nu_g).*PY_g + mu_hl.*(1-nu_h).*PY_h );
F_l = pi_l'\((eye(I)-Ups_l)*PY_l-V_l);

Ups_h = kron( ones(I,1) , mu_hh'.*(1-nu_h') ) .* pi_h';
V_h = pi_h' * ( mu_gh.*(1-nu_g).*PY_g + mu_lh.*(1-nu_l).*PY_l );
F_h = pi_h'\((eye(I)-Ups_h)*PY_h-V_h);

% sector absorption
PQ_g = F_g + mu_gg.*(1-nu_g).*PY_g + mu_lg.*(1-nu_l).*PY_l + mu_hg.*(1-nu_h).*PY_h ;
PQ_l = F_l + mu_gl.*(1-nu_g).*PY_g + mu_ll.*(1-nu_l).*PY_l + mu_hl.*(1-nu_h).*PY_h ;
PQ_h = F_h + mu_gh.*(1-nu_g).*PY_g + mu_lh.*(1-nu_l).*PY_l + mu_hh.*(1-nu_h).*PY_h ;


