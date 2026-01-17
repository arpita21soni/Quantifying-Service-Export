function TFP_decompose(I,num_year,theta,eta,r_t,w_t,P_g_t,P_l_t,P_h_t,bt_g_t,bt_l_t,bt_h_t,...
    alpha_g_t,alpha_l_t,alpha_h_t,nu_g_t,nu_l_t,nu_h_t,...
    mu_gg_t,mu_gl_t,mu_gh_t,mu_lg_t,mu_ll_t,mu_lh_t,...
    mu_hg_t,mu_hl_t,mu_hh_t)

%--------------------------------------------------------------------------
% Decompose the trade share into T, u, tao
%--------------------------------------------------------------------------

% turn on benchmark
u_g_t = ( r_t./(alpha_g_t.*nu_g_t) ).^(alpha_g_t.*nu_g_t).*...
    ( w_t./( (1-alpha_g_t).*nu_g_t) ).^( (1-alpha_g_t).*nu_g_t).*...
    ( P_g_t./(mu_gg_t.*(1-nu_g_t)) ).^(mu_gg_t.*(1-nu_g_t)).*...
    ( P_l_t./(mu_gl_t.*(1-nu_g_t)) ).^(mu_gl_t.*(1-nu_g_t)).*...
    ( P_h_t./(mu_gh_t.*(1-nu_g_t)) ).^(mu_gh_t.*(1-nu_g_t));

u_l_t = ( r_t./(alpha_l_t.*nu_l_t) ).^(alpha_l_t.*nu_l_t).*...
    ( w_t./( (1-alpha_l_t).*nu_l_t) ).^( (1-alpha_l_t).*nu_l_t).*...
    ( P_g_t./(mu_lg_t.*(1-nu_l_t)) ).^(mu_lg_t.*(1-nu_l_t)).*...
    ( P_l_t./(mu_ll_t.*(1-nu_l_t)) ).^(mu_ll_t.*(1-nu_l_t)).*...
    ( P_h_t./(mu_lh_t.*(1-nu_l_t)) ).^(mu_lh_t.*(1-nu_l_t));

u_h_t = ( r_t./(alpha_h_t.*nu_h_t) ).^(alpha_h_t.*nu_h_t).*...
    ( w_t./( (1-alpha_h_t).*nu_h_t) ).^( (1-alpha_h_t).*nu_h_t).*...
    ( P_g_t./(mu_hg_t.*(1-nu_h_t)) ).^(mu_hg_t.*(1-nu_h_t)).*...
    ( P_l_t./(mu_hl_t.*(1-nu_h_t)) ).^(mu_hl_t.*(1-nu_h_t)).*...
    ( P_h_t./(mu_hh_t.*(1-nu_h_t)) ).^(mu_hh_t.*(1-nu_h_t));

gam=gamma(1+1/theta.*(1-eta)).^(1/(1-eta));
T_g_t = zeros(I,num_year);
T_l_t = zeros(I,num_year);
T_h_t = zeros(I,num_year);

for i=1:num_year
    T_g_t(:,i) = (gam.^theta).*diag(bt_g_t(:,:,i)).*( u_g_t(:,i)./P_g_t(:,i) ).^theta;
    T_l_t(:,i) = (gam.^theta).*diag(bt_l_t(:,:,i)).*( u_l_t(:,i)./P_l_t(:,i) ).^theta;
    T_h_t(:,i) = (gam.^theta).*diag(bt_h_t(:,:,i)).*( u_h_t(:,i)./P_h_t(:,i) ).^theta;
end

% % decompose TFP
% figure('Position', get(0, 'Screensize'));
% 
% % Create year vector for 1995–2011
% years = 1995:2011;
% num_year = length(years);
% 
% subplot(2,2,1)
% plot(years, T_g_t(1,1:num_year)/T_g_t(1,1), '-b', ...
%      years, T_g_t(2,1:num_year)/T_g_t(2,1), '-r', ...
%      years, T_g_t(3,1:num_year)/T_g_t(3,1), '-g', ...
%      years, T_g_t(4,1:num_year)/T_g_t(4,1), '-c', ...
%      years, T_g_t(5,1:num_year)/T_g_t(5,1), '-m', ...
%      years, T_g_t(6,1:num_year)/T_g_t(6,1), '-k', ...
%      years, T_g_t(7,1:num_year)/T_g_t(7,1), '--b', ...
%      years, T_g_t(8,1:num_year)/T_g_t(8,1), '--r', ...
%      years, T_g_t(9,1:num_year)/T_g_t(9,1), '--g', ...
%      'LineWidth', 2.0);
% 
% title('Goods Sector')
% legend('AUS','BRA','CHN','DEU','IND','JPN','MEX','USA','ROW','Location','best')
% grid on
% 
% subplot(2,2,2)
% plot(years, T_g_t(1,1:num_year)/T_g_t(1,2), '-b', ...
%      years, T_g_t(2,1:num_year)/T_g_t(2,2), '-r', ...
%      years, T_g_t(3,1:num_year)/T_g_t(3,2), '-g', ...
%      years, T_g_t(4,1:num_year)/T_g_t(4,2), '-c', ...
%      years, T_g_t(5,1:num_year)/T_g_t(5,2), '-m', ...
%      years, T_g_t(6,1:num_year)/T_g_t(6,2), '-k', ...
%      years, T_g_t(7,1:num_year)/T_g_t(7,2), '--b', ...
%      years, T_g_t(8,1:num_year)/T_g_t(8,2), '--r', ...
%      years, T_g_t(9,1:num_year)/T_g_t(9,2), '--g', ...
%      'LineWidth', 2.0);
% 
% title('Low-skill service Sector')
% legend('AUS','BRA','CHN','DEU','IND','JPN','MEX','USA','ROW','Location','best')
% grid on
% 
% subplot(2,2,3)
% plot(years, T_g_t(1,1:num_year)/T_g_t(1,3), '-b', ...
%      years, T_g_t(2,1:num_year)/T_g_t(2,3), '-r', ...
%      years, T_g_t(3,1:num_year)/T_g_t(3,3), '-g', ...
%      years, T_g_t(4,1:num_year)/T_g_t(4,3), '-c', ...
%      years, T_g_t(5,1:num_year)/T_g_t(5,3), '-m', ...
%      years, T_g_t(6,1:num_year)/T_g_t(6,3), '-k', ...
%      years, T_g_t(7,1:num_year)/T_g_t(7,3), '--b', ...
%      years, T_g_t(8,1:num_year)/T_g_t(8,3), '--r', ...
%      years, T_g_t(9,1:num_year)/T_g_t(9,3), '--g', ...
%      'LineWidth', 2.0);
% 
% title('High-skill service Sector')
% legend('AUS','BRA','CHN','DEU','IND','JPN','MEX','USA','ROW','Location','best')
% grid on
% 
% sgtitle('Figure. Productivity by sector and country relative to 2000')
% saveas(gcf, 'plots_transition/final_2/Decompose_TFP.png');
% close(gcf)
% 
% % turn on trade_share
% T_g_t = zeros(I,num_year);
% T_l_t = zeros(I,num_year);
% T_h_t = zeros(I,num_year);
% 
% for i=1:num_year
%     T_g_t(:,i) = (gam.^theta).*diag(bt_g_t(:,:,i)).*( u_g_t(:,i)./P_g_t(:,i) ).^0;
%     T_l_t(:,i) = (gam.^theta).*diag(bt_l_t(:,:,i)).*( u_l_t(:,i)./P_l_t(:,i) ).^0;
%     T_h_t(:,i) = (gam.^theta).*diag(bt_h_t(:,:,i)).*( u_h_t(:,i)./P_h_t(:,i) ).^0;
% end
% 
% % decompose TFP
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+num_year,T_g_t(1,1:num_year)/T_g_t(1,1),'-b',2000:1999+num_year,T_g_t(2,1:num_year)/T_g_t(2,1),'-r',2000:1999+num_year,T_g_t(3,1:num_year)/T_g_t(3,1),'-g',2000:1999+num_year,T_g_t(4,1:num_year)/T_g_t(4,1),'-c',2000:1999+num_year,T_g_t(5,1:num_year)/T_g_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','IND','ROW','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+num_year,T_l_t(1,1:num_year)/T_l_t(1,1),'-b',2000:1999+num_year,T_l_t(2,1:num_year)/T_l_t(2,1),'-r',2000:1999+num_year,T_l_t(3,1:num_year)/T_l_t(3,1),'-g',2000:1999+num_year,T_l_t(4,1:num_year)/T_l_t(4,1),'-c',2000:1999+num_year,T_l_t(5,1:num_year)/T_l_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+num_year,T_h_t(1,1:num_year)/T_h_t(1,1),'-b',2000:1999+num_year,T_h_t(2,1:num_year)/T_h_t(2,1),'-r',2000:1999+num_year,T_h_t(3,1:num_year)/T_h_t(3,1),'-g',2000:1999+num_year,T_h_t(4,1:num_year)/T_h_t(4,1),'-c',2000:1999+num_year,T_h_t(5,1:num_year)/T_h_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% sgtitle('Figure. Productivity by sector and country relative to 2000')
% saveas(gcf, 'plots_transition/final_2/Decompose_TFP_turn_on_trade_share.png');
% close(gcf)


% % turn_on_unit_cost_over_price
% T_g_t = zeros(I,num_year);
% T_l_t = zeros(I,num_year);
% T_h_t = zeros(I,num_year);
% for i=1:num_year
%     T_g_t(:,i) = (gam.^theta).*diag(bt_g_t(:,:,i).^0).*( u_g_t(:,i)./P_g_t(:,i) ).^theta;
%     T_l_t(:,i) = (gam.^theta).*diag(bt_l_t(:,:,i).^0).*( u_l_t(:,i)./P_l_t(:,i) ).^theta;
%     T_h_t(:,i) = (gam.^theta).*diag(bt_h_t(:,:,i).^0).*( u_h_t(:,i)./P_h_t(:,i) ).^theta;
% end
% 
% % decompose TFP
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+num_year,T_g_t(1,1:num_year)/T_g_t(1,1),'-b',2000:1999+num_year,T_g_t(2,1:num_year)/T_g_t(2,1),'-r',2000:1999+num_year,T_g_t(3,1:num_year)/T_g_t(3,1),'-g',2000:1999+num_year,T_g_t(4,1:num_year)/T_g_t(4,1),'-c',2000:1999+num_year,T_g_t(5,1:num_year)/T_g_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','IND','ROW','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+num_year,T_l_t(1,1:num_year)/T_l_t(1,1),'-b',2000:1999+num_year,T_l_t(2,1:num_year)/T_l_t(2,1),'-r',2000:1999+num_year,T_l_t(3,1:num_year)/T_l_t(3,1),'-g',2000:1999+num_year,T_l_t(4,1:num_year)/T_l_t(4,1),'-c',2000:1999+num_year,T_l_t(5,1:num_year)/T_l_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+num_year,T_h_t(1,1:num_year)/T_h_t(1,1),'-b',2000:1999+num_year,T_h_t(2,1:num_year)/T_h_t(2,1),'-r',2000:1999+num_year,T_h_t(3,1:num_year)/T_h_t(3,1),'-g',2000:1999+num_year,T_h_t(4,1:num_year)/T_h_t(4,1),'-c',2000:1999+num_year,T_h_t(5,1:num_year)/T_h_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% 
% sgtitle('Figure. Productivity by sector and country relative to 2000')
% saveas(gcf, 'plots_transition/final_2/Decompose_TFP_turn_on_unit_cost_over_price.png');
% close(gcf)
% 
% 
% % turn_on_wage_rental
% u_g_t = ( r_t./(alpha_g_t.*nu_g_t) ).^(alpha_g_t.*nu_g_t).*...
%     ( w_t./( (1-alpha_g_t).*nu_g_t) ).^( (1-alpha_g_t).*nu_g_t).*...
%     ( P_g_t.^0./(mu_gg_t.*(1-nu_g_t)) ).^(mu_gg_t.*(1-nu_g_t)).*...
%     ( P_l_t.^0./(mu_gl_t.*(1-nu_g_t)) ).^(mu_gl_t.*(1-nu_g_t)).*...
%     ( P_h_t.^0./(mu_gh_t.*(1-nu_g_t)) ).^(mu_gh_t.*(1-nu_g_t));
% 
% u_l_t = ( r_t./(alpha_l_t.*nu_l_t) ).^(alpha_l_t.*nu_l_t).*...
%     ( w_t./( (1-alpha_l_t).*nu_l_t) ).^( (1-alpha_l_t).*nu_l_t).*...
%     ( P_g_t.^0./(mu_lg_t.*(1-nu_l_t)) ).^(mu_lg_t.*(1-nu_l_t)).*...
%     ( P_l_t.^0./(mu_ll_t.*(1-nu_l_t)) ).^(mu_ll_t.*(1-nu_l_t)).*...
%     ( P_h_t.^0./(mu_lh_t.*(1-nu_l_t)) ).^(mu_lh_t.*(1-nu_l_t));
% 
% u_h_t = ( r_t./(alpha_h_t.*nu_h_t) ).^(alpha_h_t.*nu_h_t).*...
%     ( w_t./( (1-alpha_h_t).*nu_h_t) ).^( (1-alpha_h_t).*nu_h_t).*...
%     ( P_g_t.^0./(mu_hg_t.*(1-nu_h_t)) ).^(mu_hg_t.*(1-nu_h_t)).*...
%     ( P_l_t.^0./(mu_hl_t.*(1-nu_h_t)) ).^(mu_hl_t.*(1-nu_h_t)).*...
%     ( P_h_t.^0./(mu_hh_t.*(1-nu_h_t)) ).^(mu_hh_t.*(1-nu_h_t));
% 
% 
% 
% T_g_t = zeros(I,num_year);
% T_l_t = zeros(I,num_year);
% T_h_t = zeros(I,num_year);
% 
% for i=1:num_year
%     T_g_t(:,i) = (gam.^theta).*diag(bt_g_t(:,:,i).^0).*( u_g_t(:,i)./P_g_t(:,i).^0 ).^theta;
%     T_l_t(:,i) = (gam.^theta).*diag(bt_l_t(:,:,i).^0).*( u_l_t(:,i)./P_l_t(:,i).^0  ).^theta;
%     T_h_t(:,i) = (gam.^theta).*diag(bt_h_t(:,:,i).^0).*( u_h_t(:,i)./P_h_t(:,i).^0  ).^theta;
% end
% 
% % decompose TFP
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+num_year,T_g_t(1,1:num_year)/T_g_t(1,1),'-b',2000:1999+num_year,T_g_t(2,1:num_year)/T_g_t(2,1),'-r',2000:1999+num_year,T_g_t(3,1:num_year)/T_g_t(3,1),'-g',2000:1999+num_year,T_g_t(4,1:num_year)/T_g_t(4,1),'-c',2000:1999+num_year,T_g_t(5,1:num_year)/T_g_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','IND','ROW','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+num_year,T_l_t(1,1:num_year)/T_l_t(1,1),'-b',2000:1999+num_year,T_l_t(2,1:num_year)/T_l_t(2,1),'-r',2000:1999+num_year,T_l_t(3,1:num_year)/T_l_t(3,1),'-g',2000:1999+num_year,T_l_t(4,1:num_year)/T_l_t(4,1),'-c',2000:1999+num_year,T_l_t(5,1:num_year)/T_l_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+num_year,T_h_t(1,1:num_year)/T_h_t(1,1),'-b',2000:1999+num_year,T_h_t(2,1:num_year)/T_h_t(2,1),'-r',2000:1999+num_year,T_h_t(3,1:num_year)/T_h_t(3,1),'-g',2000:1999+num_year,T_h_t(4,1:num_year)/T_h_t(4,1),'-c',2000:1999+num_year,T_h_t(5,1:num_year)/T_h_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% sgtitle('Figure. Productivity by sector and country relative to 2000')
% saveas(gcf, 'plots_transition/final_2/Decompose_TFP_turn_on_wage_rental.png');
% close(gcf)
% 
% 
% % turn_on_price
% u_g_t = ( r_t.^0./(alpha_g_t.*nu_g_t) ).^(alpha_g_t.*nu_g_t).*...
%     ( w_t.^0./( (1-alpha_g_t).*nu_g_t) ).^( (1-alpha_g_t).*nu_g_t).*...
%     ( P_g_t./(mu_gg_t.*(1-nu_g_t)) ).^(mu_gg_t.*(1-nu_g_t)).*...
%     ( P_l_t./(mu_gl_t.*(1-nu_g_t)) ).^(mu_gl_t.*(1-nu_g_t)).*...
%     ( P_h_t./(mu_gh_t.*(1-nu_g_t)) ).^(mu_gh_t.*(1-nu_g_t));
% 
% u_l_t = ( r_t.^0./(alpha_l_t.*nu_l_t) ).^(alpha_l_t.*nu_l_t).*...
%     ( w_t.^0./( (1-alpha_l_t).*nu_l_t) ).^( (1-alpha_l_t).*nu_l_t).*...
%     ( P_g_t./(mu_lg_t.*(1-nu_l_t)) ).^(mu_lg_t.*(1-nu_l_t)).*...
%     ( P_l_t./(mu_ll_t.*(1-nu_l_t)) ).^(mu_ll_t.*(1-nu_l_t)).*...
%     ( P_h_t./(mu_lh_t.*(1-nu_l_t)) ).^(mu_lh_t.*(1-nu_l_t));
% 
% u_h_t = ( r_t.^0./(alpha_h_t.*nu_h_t) ).^(alpha_h_t.*nu_h_t).*...
%     ( w_t.^0./( (1-alpha_h_t).*nu_h_t) ).^( (1-alpha_h_t).*nu_h_t).*...
%     ( P_g_t./(mu_hg_t.*(1-nu_h_t)) ).^(mu_hg_t.*(1-nu_h_t)).*...
%     ( P_l_t./(mu_hl_t.*(1-nu_h_t)) ).^(mu_hl_t.*(1-nu_h_t)).*...
%     ( P_h_t./(mu_hh_t.*(1-nu_h_t)) ).^(mu_hh_t.*(1-nu_h_t));
% 
% T_g_t = zeros(I,num_year);
% T_l_t = zeros(I,num_year);
% T_h_t = zeros(I,num_year);
% 
% for i=1:num_year
%     T_g_t(:,i) = (gam.^theta).*diag(bt_g_t(:,:,i).^0).*( u_g_t(:,i)./P_g_t(:,i) ).^theta;
%     T_l_t(:,i) = (gam.^theta).*diag(bt_l_t(:,:,i).^0).*( u_l_t(:,i)./P_l_t(:,i) ).^theta;
%     T_h_t(:,i) = (gam.^theta).*diag(bt_h_t(:,:,i).^0).*( u_h_t(:,i)./P_h_t(:,i) ).^theta;
% end
% 
% % decompose TFP
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+num_year,T_g_t(1,1:num_year)/T_g_t(1,1),'-b',2000:1999+num_year,T_g_t(2,1:num_year)/T_g_t(2,1),'-r',2000:1999+num_year,T_g_t(3,1:num_year)/T_g_t(3,1),'-g',2000:1999+num_year,T_g_t(4,1:num_year)/T_g_t(4,1),'-c',2000:1999+num_year,T_g_t(5,1:num_year)/T_g_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','IND','ROW','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+num_year,T_l_t(1,1:num_year)/T_l_t(1,1),'-b',2000:1999+num_year,T_l_t(2,1:num_year)/T_l_t(2,1),'-r',2000:1999+num_year,T_l_t(3,1:num_year)/T_l_t(3,1),'-g',2000:1999+num_year,T_l_t(4,1:num_year)/T_l_t(4,1),'-c',2000:1999+num_year,T_l_t(5,1:num_year)/T_l_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+num_year,T_h_t(1,1:num_year)/T_h_t(1,1),'-b',2000:1999+num_year,T_h_t(2,1:num_year)/T_h_t(2,1),'-r',2000:1999+num_year,T_h_t(3,1:num_year)/T_h_t(3,1),'-g',2000:1999+num_year,T_h_t(4,1:num_year)/T_h_t(4,1),'-c',2000:1999+num_year,T_h_t(5,1:num_year)/T_h_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% 
% sgtitle('Figure. Productivity by sector and country relative to 2000')
% saveas(gcf, 'plots_transition/final_2/Decompose_TFP_turn_on_price.png');
% close(gcf)
% 
% % turn_on_wage
% u_g_t = ( r_t.^0./(alpha_g_t.*nu_g_t) ).^(alpha_g_t.*nu_g_t).*...
%     ( w_t./( (1-alpha_g_t).*nu_g_t) ).^( (1-alpha_g_t).*nu_g_t).*...
%     ( P_g_t.^0./(mu_gg_t.*(1-nu_g_t)) ).^(mu_gg_t.*(1-nu_g_t)).*...
%     ( P_l_t.^0./(mu_gl_t.*(1-nu_g_t)) ).^(mu_gl_t.*(1-nu_g_t)).*...
%     ( P_h_t.^0./(mu_gh_t.*(1-nu_g_t)) ).^(mu_gh_t.*(1-nu_g_t));
% 
% u_l_t = ( r_t.^0./(alpha_l_t.*nu_l_t) ).^(alpha_l_t.*nu_l_t).*...
%     ( w_t./( (1-alpha_l_t).*nu_l_t) ).^( (1-alpha_l_t).*nu_l_t).*...
%     ( P_g_t.^0./(mu_lg_t.*(1-nu_l_t)) ).^(mu_lg_t.*(1-nu_l_t)).*...
%     ( P_l_t.^0./(mu_ll_t.*(1-nu_l_t)) ).^(mu_ll_t.*(1-nu_l_t)).*...
%     ( P_h_t.^0./(mu_lh_t.*(1-nu_l_t)) ).^(mu_lh_t.*(1-nu_l_t));
% 
% u_h_t = ( r_t.^0./(alpha_h_t.*nu_h_t) ).^(alpha_h_t.*nu_h_t).*...
%     ( w_t./( (1-alpha_h_t).*nu_h_t) ).^( (1-alpha_h_t).*nu_h_t).*...
%     ( P_g_t.^0./(mu_hg_t.*(1-nu_h_t)) ).^(mu_hg_t.*(1-nu_h_t)).*...
%     ( P_l_t.^0./(mu_hl_t.*(1-nu_h_t)) ).^(mu_hl_t.*(1-nu_h_t)).*...
%     ( P_h_t.^0./(mu_hh_t.*(1-nu_h_t)) ).^(mu_hh_t.*(1-nu_h_t));
% 
% T_g_t = zeros(I,num_year);
% T_l_t = zeros(I,num_year);
% T_h_t = zeros(I,num_year);
% 
% for i=1:num_year
%     T_g_t(:,i) = (gam.^theta).*diag(bt_g_t(:,:,i).^0).*( u_g_t(:,i)./P_g_t(:,i).^0 ).^theta;
%     T_l_t(:,i) = (gam.^theta).*diag(bt_l_t(:,:,i).^0).*( u_l_t(:,i)./P_l_t(:,i).^0  ).^theta;
%     T_h_t(:,i) = (gam.^theta).*diag(bt_h_t(:,:,i).^0).*( u_h_t(:,i)./P_h_t(:,i).^0  ).^theta;
% end
% 
% % decompose TFP
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+num_year,T_g_t(1,1:num_year)/T_g_t(1,1),'-b',2000:1999+num_year,T_g_t(2,1:num_year)/T_g_t(2,1),'-r',2000:1999+num_year,T_g_t(3,1:num_year)/T_g_t(3,1),'-g',2000:1999+num_year,T_g_t(4,1:num_year)/T_g_t(4,1),'-c',2000:1999+num_year,T_g_t(5,1:num_year)/T_g_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','IND','ROW','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+num_year,T_l_t(1,1:num_year)/T_l_t(1,1),'-b',2000:1999+num_year,T_l_t(2,1:num_year)/T_l_t(2,1),'-r',2000:1999+num_year,T_l_t(3,1:num_year)/T_l_t(3,1),'-g',2000:1999+num_year,T_l_t(4,1:num_year)/T_l_t(4,1),'-c',2000:1999+num_year,T_l_t(5,1:num_year)/T_l_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+num_year,T_h_t(1,1:num_year)/T_h_t(1,1),'-b',2000:1999+num_year,T_h_t(2,1:num_year)/T_h_t(2,1),'-r',2000:1999+num_year,T_h_t(3,1:num_year)/T_h_t(3,1),'-g',2000:1999+num_year,T_h_t(4,1:num_year)/T_h_t(4,1),'-c',2000:1999+num_year,T_h_t(5,1:num_year)/T_h_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% 
% sgtitle('Figure. Productivity by sector and country relative to 2000')
% saveas(gcf, 'plots_transition/final_2/Decompose_TFP_turn_on_wage.png');
% close(gcf)
% 
% % turn_on_rental
% u_g_t = ( r_t./(alpha_g_t.*nu_g_t) ).^(alpha_g_t.*nu_g_t).*...
%     ( w_t.^0./( (1-alpha_g_t).*nu_g_t) ).^( (1-alpha_g_t).*nu_g_t).*...
%     ( P_g_t.^0./(mu_gg_t.*(1-nu_g_t)) ).^(mu_gg_t.*(1-nu_g_t)).*...
%     ( P_l_t.^0./(mu_gl_t.*(1-nu_g_t)) ).^(mu_gl_t.*(1-nu_g_t)).*...
%     ( P_h_t.^0./(mu_gh_t.*(1-nu_g_t)) ).^(mu_gh_t.*(1-nu_g_t));
% 
% u_l_t = ( r_t./(alpha_l_t.*nu_l_t) ).^(alpha_l_t.*nu_l_t).*...
%     ( w_t.^0./( (1-alpha_l_t).*nu_l_t) ).^( (1-alpha_l_t).*nu_l_t).*...
%     ( P_g_t.^0./(mu_lg_t.*(1-nu_l_t)) ).^(mu_lg_t.*(1-nu_l_t)).*...
%     ( P_l_t.^0./(mu_ll_t.*(1-nu_l_t)) ).^(mu_ll_t.*(1-nu_l_t)).*...
%     ( P_h_t.^0./(mu_lh_t.*(1-nu_l_t)) ).^(mu_lh_t.*(1-nu_l_t));
% 
% u_h_t = ( r_t./(alpha_h_t.*nu_h_t) ).^(alpha_h_t.*nu_h_t).*...
%     ( w_t.^0./( (1-alpha_h_t).*nu_h_t) ).^( (1-alpha_h_t).*nu_h_t).*...
%     ( P_g_t.^0./(mu_hg_t.*(1-nu_h_t)) ).^(mu_hg_t.*(1-nu_h_t)).*...
%     ( P_l_t.^0./(mu_hl_t.*(1-nu_h_t)) ).^(mu_hl_t.*(1-nu_h_t)).*...
%     ( P_h_t.^0./(mu_hh_t.*(1-nu_h_t)) ).^(mu_hh_t.*(1-nu_h_t));
% 
% T_g_t = zeros(I,num_year);
% T_l_t = zeros(I,num_year);
% T_h_t = zeros(I,num_year);
% 
% for i=1:num_year
%     T_g_t(:,i) = (gam.^theta).*diag(bt_g_t(:,:,i).^0).*( u_g_t(:,i)./P_g_t(:,i).^0 ).^theta;
%     T_l_t(:,i) = (gam.^theta).*diag(bt_l_t(:,:,i).^0).*( u_l_t(:,i)./P_l_t(:,i).^0  ).^theta;
%     T_h_t(:,i) = (gam.^theta).*diag(bt_h_t(:,:,i).^0).*( u_h_t(:,i)./P_h_t(:,i).^0  ).^theta;
% end

% decompose TFP
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+num_year,T_g_t(1,1:num_year)/T_g_t(1,1),'-b',2000:1999+num_year,T_g_t(2,1:num_year)/T_g_t(2,1),'-r',2000:1999+num_year,T_g_t(3,1:num_year)/T_g_t(3,1),'-g',2000:1999+num_year,T_g_t(4,1:num_year)/T_g_t(4,1),'-c',2000:1999+num_year,T_g_t(5,1:num_year)/T_g_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','IND','ROW','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+num_year,T_l_t(1,1:num_year)/T_l_t(1,1),'-b',2000:1999+num_year,T_l_t(2,1:num_year)/T_l_t(2,1),'-r',2000:1999+num_year,T_l_t(3,1:num_year)/T_l_t(3,1),'-g',2000:1999+num_year,T_l_t(4,1:num_year)/T_l_t(4,1),'-c',2000:1999+num_year,T_l_t(5,1:num_year)/T_l_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+num_year,T_h_t(1,1:num_year)/T_h_t(1,1),'-b',2000:1999+num_year,T_h_t(2,1:num_year)/T_h_t(2,1),'-r',2000:1999+num_year,T_h_t(3,1:num_year)/T_h_t(3,1),'-g',2000:1999+num_year,T_h_t(4,1:num_year)/T_h_t(4,1),'-c',2000:1999+num_year,T_h_t(5,1:num_year)/T_h_t(5,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% 
% sgtitle('Figure. Productivity by sector and country relative to 2000')
% saveas(gcf, 'plots_transition/final_2/Decompose_TFP_turn_on_rental.png');
% close(gcf)





