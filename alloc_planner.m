function [w,r,P_g,P_l,P_h,P_c,P_x,C,X,K,A,F, ...
    K_g,K_l,K_h,L_g,L_l,L_h,Y_g,Y_l,Y_h,pi_g,pi_l,pi_h] = ...
    alloc_planner(w,A_x,T_g,T_l,T_h,d_g,d_l,d_h,L,K_1,A_1,I,TT,...
    theta,eta,del,lambda,epsilon_c,epsilon_x,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,...
    rho,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,...
    fi,omega_bar,omega_cg,omega_cl,omega_ch,omega_xg,omega_xl,omega_xh)

%--------------------------------------------------------------------------
% This function computes prices and allocations as a function of wages at a
% point in time along the transition path. It calls the functions:
%--------------------------------------------------------------------------

% Given wages, update remaining prices, trade shares, investment and path
% for captal stocks.
K=zeros(I,TT+1);
X=zeros(I,TT);
r=zeros(I,TT);
P_g=zeros(I,TT);
P_l=zeros(I,TT);
P_h=zeros(I,TT);
P_x=zeros(I,TT);
pi_g=zeros(I,I,TT);
pi_l=zeros(I,I,TT);
pi_h=zeros(I,I,TT);
Y_g=zeros(I,TT);
Y_l=zeros(I,TT);
Y_h=zeros(I,TT);

K(:,1) = K_1;

for tt=1:TT
    if tt==1
        Pguess_g=ones(I,1);
        Pguess_l=ones(I,1);
        Pguess_h=ones(I,1);
    else
        Pguess_g=P_g(:,tt-1);
        Pguess_l=P_l(:,tt-1);
        Pguess_h=P_h(:,tt-1);
    end
    
    [Y_g(:,tt),Y_l(:,tt),Y_h(:,tt),w(:,tt),r(:,tt),P_g(:,tt),P_l(:,tt),P_h(:,tt),pi_g(:,:,tt),pi_l(:,:,tt),pi_h(:,:,tt)]=prodn_alloc_iter_transition(...
        K(:,tt),rho(:,tt),w(:,tt),T_g(:,tt),T_l(:,tt),T_h(:,tt),d_g(:,:,tt),d_l(:,:,tt),d_h(:,:,tt),L(:,tt),I,theta,eta,...
        alpha_g(:,tt),alpha_l(:,tt),alpha_h(:,tt),nu_g(:,tt),nu_l(:,tt),nu_h(:,tt),Pguess_g,Pguess_l,Pguess_h,...
        mu_gg(:,tt),mu_gl(:,tt),mu_gh(:,tt),mu_lg(:,tt),mu_ll(:,tt),mu_lh(:,tt),...
        mu_hg(:,tt),mu_hl(:,tt),mu_hh(:,tt),...
        fi(:,tt),omega_bar,omega_cg(:,tt),omega_cl(:,tt),omega_ch(:,tt),omega_xg(:,tt),omega_xl(:,tt),omega_xh(:,tt));
    
    if epsilon_x == 1
        P_x(:,tt) = (P_g(:,tt)./omega_xg(:,tt)).^omega_xg(:,tt).*(P_l(:,tt)./omega_xl(:,tt)).^omega_xl(:,tt).*(P_h(:,tt)./omega_xh(:,tt)).^omega_xh(:,tt);
    else
        P_x(:,tt) = ( omega_xg(:,tt).*P_g(:,tt).^(1-epsilon_x)+omega_xl(:,tt).*P_l(:,tt).^(1-epsilon_x)+omega_xh(:,tt).*P_h(:,tt).^(1-epsilon_x)).^(1/(1-epsilon_x));
    end
    X(:,tt) = rho(:,tt).* (w(:,tt).*L(:,tt)+r(:,tt).*K(:,tt))./ P_x(:,tt);
    K(:,tt+1) = (1-del).*K(:,tt) + A_x(:,tt).*X(:,tt).^lambda.*K(:,tt).^(1-lambda);
    
end


% Intertemporal houehold optimization
if epsilon_c == 1
    P_c = (P_g./omega_cg).^omega_cg.*(P_l./omega_cl).^omega_cl.*(P_h./omega_ch).^omega_ch;
else
    P_c = ( omega_cg.*P_g.^(1-epsilon_c)+omega_cl.*P_l.^(1-epsilon_c)+omega_ch.*P_h.^(1-epsilon_c) ).^(1/(1-epsilon_c));
end
%C = repmat(omega_bar.*fi,1,TT)./P_c;
%C = repmat(omega_bar,1,TT).*fi./P_c;

%transfer = repmat( sum( fi.*( w.*L + r.*K(:,1:TT) ) ),I,1 ).*L./repmat( sum(L),I,1 );
transfer = repmat( sum( fi.*( w.*L + r.*K(:,1:TT) ) ),I,1 ).*( w.*L + r.*K(:,1:TT) )./repmat( sum( w.*L + r.*K(:,1:TT) ),I,1 );
%transfer = repmat( sum( fi.*( w.*L + r.*K(:,1:TT) ) ),I,1 )/I;
C = ( (1-fi).*(w.*L + r.*K(:,1:TT))+transfer-P_x.*X )./P_c;
A = w.*L + r.*K(:,1:TT) - P_x.*X - P_c.*C;


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
F = Spnd_g + Spnd_l + Spnd_h - P_g.*Y_g - P_l.*Y_l - P_h.*Y_h;






