function [P_g,P_l,P_h,r]=SS_priceindex_planner(w,A_x,T_g,T_l,...
                    T_h,d_g,d_l,d_h,I,theta,eta,lambda,del,epsilon_x,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,beta,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,...
                    omega_xg,omega_xl,omega_xh)
                
%--------------------------------------------------------------------------
% Computes a relevant price index for each country.
%--------------------------------------------------------------------------

gam=gamma(1+1/theta.*(1-eta)).^(1/(1-eta));

P_g = ones(I,1);
P_l = ones(I,1);
P_h = ones(I,1);

PHI1 = 1./lambda.*del.^( (1-lambda)./lambda ).*A_x.^(-1./lambda);
PHI2 = del.^( (1-lambda)./lambda ).*(del-1./lambda).*A_x.^(-1./lambda);

if epsilon_x == 1
    P_x = (P_g./omega_xg).^omega_xg.*(P_l./omega_xl).^omega_xl.*(P_h./omega_xh).^omega_xh;
else
    P_x = (omega_xg.*P_g.^(1-epsilon_x)+omega_xl.*P_l.^(1-epsilon_x)+omega_xh.*P_h.^(1-epsilon_x)).^(1/(1-epsilon_x));
end
r = ( PHI1/beta + PHI2 ) .* P_x; % rental rate in the ss
   
dif=10;
iter=0;
% damp = 0.2;
% prev_dif = dif;
while dif>1e-6
    iter=iter+1;
    
% % --- inside the while dif>1e-6 loop, replace the three u_* blocks with: ---
% 
% % Clamp shares to avoid exact 0/1 edge cases
% epsS = 1e-12;
% alpha_g = min(max(alpha_g,epsS),1-epsS);
% alpha_h = min(max(alpha_h,epsS),1-epsS);
% alpha_l = min(max(alpha_l,epsS),1-epsS);
% nu_g    = min(max(nu_g,epsS),1-epsS);
% nu_h    = min(max(nu_h,epsS),1-epsS);
% nu_l    = min(max(nu_l,epsS),1-epsS);
% 
% Sanity on prices & factor prices
% if any(r<=0|~isfinite(r)) || any(w<=0|~isfinite(w)) ...
%    || any(P_g<=0|~isfinite(P_g)) || any(P_h<=0|~isfinite(P_h)) || any(P_l<=0|~isfinite(P_l))
% 
%     dif
%     P_g
%     P_h
%     P_l
%     error('Nonpositive/NaN factor or sector price in unit-cost step.');
% end
% 
% % ---- g sector ----
% aK = alpha_g.*nu_g;            aL = (1-alpha_g).*nu_g;
% aG = mu_gg.*(1-nu_g);          aH = mu_gh.*(1-nu_g);   aLw = mu_gl.*(1-nu_g);
% 
% log_u_g = aK.*log(r./(alpha_g.*nu_g)) ...
%         + aL.*log(w./((1-alpha_g).*nu_g)) ...
%         + aG.*log(P_g) + aH.*log(P_h) + aLw.*log(P_l);
% u_g = exp(log_u_g);
% 
% part1_g   = u_g.^(-theta).*T_g;               % I×1
% P_gnew    = gam.*sum( d_g.^(-theta).*part1_g', 2 ).^(-1/theta);
% 
% % ---- l sector ----
% aK = alpha_l.*nu_l;            aL = (1-alpha_l).*nu_l;
% aG = mu_lg.*(1-nu_l);          aH = mu_lh.*(1-nu_l);   aLw = mu_ll.*(1-nu_l);
% 
% log_u_l = aK.*log(r./(alpha_l.*nu_l)) ...
%         + aL.*log(w./((1-alpha_l).*nu_l)) ...
%         + aG.*log(P_g) + aH.*log(P_h) + aLw.*log(P_l);
% u_l = exp(log_u_l);
% 
% part1_l   = u_l.^(-theta).*T_l;
% P_lnew    = gam.*sum( d_l.^(-theta).*part1_l', 2 ).^(-1/theta);
% 
% % ---- h sector ----
% aK = alpha_h.*nu_h;            aL = (1-alpha_h).*nu_h;
% aG = mu_hg.*(1-nu_h);          aH = mu_hh.*(1-nu_h);   aLw = mu_hl.*(1-nu_h);
% 
% log_u_h = aK.*log(r./(alpha_h.*nu_h)) ...
%         + aL.*log(w./((1-alpha_h).*nu_h)) ...
%         + aG.*log(P_g) + aH.*log(P_h) + aLw.*log(P_l);
% u_h = exp(log_u_h);
% 
% part1_h   = u_h.^(-theta).*T_h;
% P_hnew    = gam.*sum( d_h.^(-theta).*part1_h', 2 ).^(-1/theta);
% 
% % ---- investment price & rental rate update (keep your aggregator) ----
% if epsilon_x == 1
%     % (If omega are CD budget shares, the standard CD price is prod P^omega;
%     % your (P/omega)^omega form is acceptable only if you're using a normalized
%     % expenditure function with taste shifters. Keep it consistent with elsewhere.)
%     P_xnew = (P_gnew./omega_xg).^omega_xg .* (P_lnew./omega_xl).^omega_xl .* (P_hnew./omega_xh).^omega_xh;
% else
%     P_xnew = ( omega_xg.*P_gnew.^(1-epsilon_x) + omega_xl.*P_lnew.^(1-epsilon_x) + omega_xh.*P_hnew.^(1-epsilon_x) ).^(1/(1-epsilon_x));
% end
% rnew = ( PHI1/beta + PHI2 ) .* P_xnew;





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
    
    if epsilon_x == 1
        P_xnew = (P_gnew./omega_xg).^omega_xg.*(P_lnew./omega_xl).^omega_xl.*(P_hnew./omega_xh).^omega_xh;
    else
        P_xnew = ( omega_xg.*P_gnew.^(1-epsilon_x)+omega_xl.*P_lnew.^(1-epsilon_x)+omega_xh.*P_hnew.^(1-epsilon_x)).^(1/(1-epsilon_x));
    end
    rnew = ( PHI1/beta + PHI2 ) .* P_xnew;
    
    dif = max(abs([P_gnew; P_lnew; P_hnew; rnew]-[P_g; P_l; P_h; r])./[P_g; P_l; P_h; r]); 
    smooth=0.05;
    P_g=smooth*P_gnew+(1-smooth)*P_g; 
    P_l=smooth*P_lnew+(1-smooth)*P_l;
    P_h=smooth*P_hnew+(1-smooth)*P_h;
    r=smooth*rnew+(1-smooth)*r; 

% --- inside the loop, after computing P_anew,..., rnew ---
% logP_g = log(P_g);  
% logP_l = log(P_l);  
% logP_h = log(P_h);  
% logr = log(r);
% logP_gnew = log(P_gnew); 
% logP_lnew = log(P_lnew); 
% logP_hnew = log(P_hnew); 
% logrnew = log(rnew);
% 
% logP_g = (1-damp)*logP_g + damp*logP_gnew;
% logP_l = (1-damp)*logP_l + damp*logP_lnew;
% logP_h = (1-damp)*logP_h + damp*logP_hnew;
% logr   = (1-damp)*logr   + damp*logrnew;
% 
% P_g = exp(logP_g);
% P_l = exp(logP_l); 
% P_h = exp(logP_h); 
% r = exp(logr);
% 
% dif = max(abs([logP_gnew;logP_lnew;logP_hnew;logrnew] - ...
%               [logP_g;   logP_l;   logP_h;   logr]));
% 
% 
% if iter>1 && dif > 1.05*prev_dif
%     damp = max(0.05, 0.5*damp);
% end
% prev_dif = dif;


end

              