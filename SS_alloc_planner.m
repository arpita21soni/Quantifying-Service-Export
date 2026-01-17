function [w,r,q,P_g,P_l,P_h,P_c,P_x,C,X,K,A,F,...
            K_g,K_l,K_h,L_g,L_l,L_h,Y_g,Y_l,Y_h,pi_g,pi_l,pi_h] = ...
        SS_alloc_planner(w,A_x,T_g,T_l,T_h,d_g,d_l,d_h,L,I,...
                theta,eta,del,lambda,epsilon_c,epsilon_x,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,...
                                beta,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,...
                                fi2,omega_bar,omega_cg,omega_cl,omega_ch,omega_xg,omega_xl,omega_xh)

%--------------------------------------------------------------------------
% This function computes prices and allocations as a function of wages 
% in steady state.
%--------------------------------------------------------------------------

% Given wages, update remaining prices and trade shares.
[P_g,P_l,P_h,r]=SS_priceindex_planner(w,A_x,T_g,T_l,T_h,d_g,d_l,d_h,I,theta,eta,lambda,del,epsilon_x,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,beta,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,omega_xg,omega_xl,omega_xh);
             
q = 1./beta-1;
[pi_g,pi_l,pi_h] = tradeshare_planner(w,r,P_g,P_l,P_h,T_g,T_l,T_h,d_g,d_l,d_h,I,theta,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh);

% Given prices, compute sector allocation
if epsilon_c == 1
    P_c = (P_g./omega_cg).^omega_cg.*(P_l./omega_cl).^omega_cl.*(P_h./omega_ch).^omega_ch;
else
    P_c = ( omega_cg.*P_g.^(1-epsilon_c)+omega_cl.*P_l.^(1-epsilon_c)+omega_ch.*P_h.^(1-epsilon_c) );
end

if epsilon_x == 1
    P_x = (P_g./omega_xg).^omega_xg.*(P_l./omega_xl).^omega_xl.*(P_h./omega_xh).^omega_xh;
else
    P_x = ( omega_xg.*P_g.^(1-epsilon_x)+omega_xl.*P_l.^(1-epsilon_x)+omega_xh.*P_h.^(1-epsilon_x) );
end

[Y_g,Y_l,Y_h,K,X,C]=prodn_alloc_iter_planner(A_x,pi_g,pi_l,pi_h,w,r,P_g,P_l,P_h,P_c,P_x,L,I,lambda,del,epsilon_c,epsilon_x,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,...
    fi2,omega_bar,omega_cg,omega_cl,omega_ch,omega_xg,omega_xl,omega_xh);
 
% denominate prices in WGDP
numeraire = sum(P_c.*C+P_x.*X);
w = w/numeraire;
r = r/numeraire;
P_g = P_g/numeraire;
P_l = P_l/numeraire;
P_h = P_h/numeraire;

P_c = P_c/numeraire;
P_x = P_x/numeraire;

% net foreign asset or net transfer
A = w.*L+r.*K-P_c.*C-P_x.*X;

% Labor used in each sector
L_g = (1-alpha_g).*nu_g .* P_g.*Y_g ./ w;
L_l = (1-alpha_l).*nu_l .* P_l.*Y_l ./ w;
L_h = (1-alpha_h).*nu_h .* P_h.*Y_h ./ w;


% Capital stock used in each sector
K_g = alpha_g.*nu_g .* P_g.*Y_g ./ r;
K_l = alpha_l.*nu_l .* P_l.*Y_l ./ r;
K_h = alpha_h.*nu_h .* P_h.*Y_h ./ r;



% Gross spending in each sector for intermediate and final use
Spnd_g = mu_gg.*(1-nu_g).*P_g.*Y_g+mu_lg.*(1-nu_l).*P_l.*Y_l+mu_hg.*(1-nu_h).*P_h.*Y_h+omega_cg.*(P_g./P_c).^(1-epsilon_c).*P_c.*C+omega_xg.*( P_g./P_x ).^(1-epsilon_x).*P_x.*X;
Spnd_l = mu_gl.*(1-nu_g).*P_g.*Y_g+mu_ll.*(1-nu_l).*P_l.*Y_l+mu_hl.*(1-nu_h).*P_h.*Y_h+omega_cl.*(P_l./P_c).^(1-epsilon_c).*P_c.*C+omega_xl.*( P_l./P_x ).^(1-epsilon_x).*P_x.*X;
Spnd_h = mu_gh.*(1-nu_g).*P_g.*Y_g+mu_lh.*(1-nu_l).*P_l.*Y_l+mu_hh.*(1-nu_h).*P_h.*Y_h+omega_ch.*(P_h./P_c).^(1-epsilon_c).*P_c.*C+omega_xh.*( P_h./P_x ).^(1-epsilon_x).*P_x.*X;


% Trade deficit
F = Spnd_g + Spnd_l + Spnd_h - P_g.*Y_g - P_l.*Y_l - P_h.*Y_h ;


% GDP = P_c.*C + P_x.*X;
% % Diagnostics
% ca_share = A ./ GDP;
% tb_share = F ./ GDP;
% gap_share = (A - F) ./ GDP;
% 
% % Assert feasibility
% if any(~isfinite([A;F;GDP])) || any(GDP <= 0)
%     error('NaN/Inf or nonpositive GDP encountered.');
% end
% 
% % Soft warnings (comment out once stable)
% if max(abs(ca_share)) > 0.5 || max(abs(tb_share)) > 0.5
%     warning('Large external shares: max |A/GDP| = %.3f, |F/GDP| = %.3f', ...
%             max(abs(ca_share)), max(abs(tb_share)));
% end
% 
% if max(abs(gap_share)) > 1e-4
%     % not fatal mid-iteration, but should go to ~1e-6 at convergence
%     % fprintf('Gap A-F (share of GDP) max = %.2e\n', max(abs(gap_share)));
% end




