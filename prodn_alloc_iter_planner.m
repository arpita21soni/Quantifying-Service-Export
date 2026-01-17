function [Y_g,Y_l,Y_h,K,X,C]=prodn_alloc_iter_planner(A_x,pi_g,pi_l,pi_h,w,r,P_g,P_l,P_h,P_c,P_x,L,I,lambda,del,epsilon_c,epsilon_x,alpha_g,alpha_l,alpha_h,nu_g,nu_l,nu_h,mu_gg,mu_gl,mu_gh,mu_lg,mu_ll,mu_lh,mu_hg,mu_hl,mu_hh,...
    fi2,omega_bar,omega_cg,omega_cl,omega_ch,omega_xg,omega_xl,omega_xh)
%--------------------------------------------------------------------------
% Given prices, trade shares, consumption, and investment this function
% computes production of intermediates at a point in time in a country.
%--------------------------------------------------------------------------

Y_g = ones(I,1);
Y_l = ones(I,1);
Y_h = ones(I,1);

K = (alpha_g.*nu_g.*P_g.*Y_g+alpha_l.*nu_l.*P_l.*Y_l+alpha_h.*nu_h.*P_h.*Y_h)./r;
X = (del./A_x).^(1./lambda).*K;
% WPC = sum(w.*L + r.*K - P_x.*X); % sum of net transfer is zero
% C = WPC*omega_bar.*fi2./P_c; % net transfer is in terms of WPC
%transfer = sum( fi2.*(w.*L + r.*K) ).*L./sum(L);
transfer = sum( fi2.*(w.*L + r.*K) ).*(w.*L + r.*K)./sum(w.*L + r.*K);
%transfer = sum( fi2.*(w.*L + r.*K) )/I;
C = ( (1-fi2).*(w.*L + r.*K)+transfer-P_x.*X )./P_c;

dif=10;
iter=0;
while dif>1e-6
    iter=iter+1;
% 
%     bad = ~isfinite(P_g) | P_g<=0;
% if any(bad)
%     disp('Indices with bad P_g:'); disp(find(bad).');
%     disp('Offending P_g values:'); disp(P_g(bad).');
% end
% epsP = 1e-8;
% P_g(~isfinite(P_g) | P_g<=0) = epsP;
% disp(P_g);
      
    Ups_g =  kron( ones(1,I) , 1./P_g ) ...
        .* kron( ones(I,1) , mu_gg'.*(1-nu_g').*P_g' ) .* pi_g'; % 3*3
    V_g = ( kron( ones(1,I) , 1./P_g ) .* pi_g' ) * ( mu_lg.*(1-nu_l).*P_l.*Y_l + mu_hg.*(1-nu_h).*P_h.*Y_h ); % 3*1
    F_g = ( kron( ones(1,I) , 1./P_g ) .* pi_g') * (omega_cg.*(P_g./P_c).^(1-epsilon_c).*P_c.*C+omega_xg.*( P_g./P_x ).^(1-epsilon_x).*P_x.*X); % final consumption demand
    Y_gnew = (eye(I)-Ups_g)\(V_g+F_g);
    
    Ups_l =  kron( ones(1,I) , 1./P_l ) ...
        .* kron( ones(I,1) , mu_ll'.*(1-nu_l').*P_l' ) .* pi_l';
    V_l = ( kron( ones(1,I) , 1./P_l ) .* pi_l' ) * ( mu_gl.*(1-nu_g).*P_g.*Y_g + mu_hl.*(1-nu_h).*P_h.*Y_h);
    F_l = ( kron( ones(1,I) , 1./P_l ) .* pi_l') * (omega_cl.*(P_l./P_c).^(1-epsilon_c).*P_c.*C+omega_xl.*( P_l./P_x ).^(1-epsilon_x).*P_x.*X); % final consumption demand
    Y_lnew = (eye(I)-Ups_l)\(V_l+F_l);
    
    Ups_h =  kron( ones(1,I) , 1./P_h ) ...
        .* kron( ones(I,1) , mu_hh'.*(1-nu_h').*P_h' ) .* pi_h';
    V_h = ( kron( ones(1,I) , 1./P_h ) .* pi_h' ) * ( mu_gh.*(1-nu_g).*P_g.*Y_g + mu_lh.*(1-nu_l).*P_l.*Y_l );
    F_h = ( kron( ones(1,I) , 1./P_h ) .* pi_h') * (omega_ch.*(P_h./P_c).^(1-epsilon_c).*P_c.*C+omega_xh.*( P_h./P_x ).^(1-epsilon_x).*P_x.*X); % final consumption demand
    Y_hnew = (eye(I)-Ups_h)\(V_h+F_h);
         
    dif = max(abs([Y_gnew; Y_lnew; Y_hnew]-[Y_g; Y_l; Y_h])./[Y_g; Y_l; Y_h]);
    smooth=.1*rand+.9;
    Y_g=smooth*Y_gnew+(1-smooth)*Y_g;
    Y_l=smooth*Y_lnew+(1-smooth)*Y_l;
    Y_h=smooth*Y_hnew+(1-smooth)*Y_h;
    
    K = (alpha_g.*nu_g.*P_g.*Y_g+alpha_l.*nu_l.*P_l.*Y_l+alpha_h.*nu_h.*P_h.*Y_h)./r;
    X = (del./A_x).^(1./lambda).*K;
%     WPC = sum(w.*L + r.*K - P_x.*X);
%     C = WPC*omega_bar.*fi2./P_c;
    %transfer = sum( fi2.*(w.*L + r.*K) ).*L./sum(L);
    transfer = sum( fi2.*(w.*L + r.*K) ).*(w.*L + r.*K)./sum(w.*L + r.*K);
    %transfer = sum( fi2.*(w.*L + r.*K) )/I;
    C = ( (1-fi2).*(w.*L + r.*K)+transfer-P_x.*X )./P_c;
    
end



