clear;
clc;


% -----------------------------------------------------------------------
% Compute transition path from initial steady state to new steady state
% -----------------------------------------------------------------------

set(0, 'DefaultFigureVisible', 'off');
mkdir plots_transition
mkdir plots_transition\final_2



%% Paramaters

% -----------------------------------------------------------------------
% Set control parameters
% -----------------------------------------------------------------------
Time_varying_IO_table = 0; % if 0, default IO table year is 2001
                           % if 1, time-varying IO table
                           % if 2, invariant IO table using the avg between 2000-2014
Time_varying_sec_invest = 0; % if 0, invariant sector investment share, default is avg across time
                             % if 1, time-varying sector investment share
Same_capital_share = 0; % if 1, capital share is the same across sectors and countries
Same_va_share = 0;      % if 1, value added share isthe same across sectors and countries 
Same_int_demand = 0;    % if 1, same intermediate demand across sectors and countries
Sector_specific_share = 0; 

Last_data_year = 2011;
Load_transition_whole = 1; % 0 means load whole transition from 1995 to ss (used 0 when ran first time)
Save_transition = 0; % 1 means save transition but used it when ran first time now dont need it

% -----------------------------------------------------------------------
% Set parameter values, and load data moments for calibration
% For four country setting, think about country 1 for LAC, country 2 for
% CHN, country 3 for US, country 4 for IND, and country 5 for ROW
% MINE: 1 is AUS, 2 is BRA, 3 is CHN, 4 is DEU, 5 is IND, 6 is JPN, 7 is
% MEX, 8 is USA
% ------------------------------------------------------------del-----------


I=9;                    % num of countries
theta=2;                % dispersion of productivity
eta=2;                  % elsaticity of subs b/w varieties
beta=0.96;              % annual discount factor
sigma=1;                % Intertemporal elasticity of substitution
epsilon_c=1;            % Substitution for consumption good
epsilon_x=1;            % Substitution for investment good

% country specific
del=0.06*ones(I,1);     % annual capital stock depreciation rate
lambda=0.55*ones(I,1);  % adjustment cost elasticity


% declare variables to load (a is g now, and m is l, h is h)
% alpha is my beta_K, nu is va/go, mu is IO_GH/total_IO_G
num_year = Last_data_year-1995+1;
T = num_year;

exr = zeros(I,num_year);

alpha_g_t = zeros(I,num_year);
alpha_l_t = zeros(I,num_year);
alpha_h_t = zeros(I,num_year);

nu_g_t = zeros(I,num_year);
nu_l_t = zeros(I,num_year);
nu_h_t = zeros(I,num_year);

va_g_t = zeros(I,num_year);
va_l_t = zeros(I,num_year);
va_h_t = zeros(I,num_year);

mu_gg_t = zeros(I,num_year);
mu_gl_t = zeros(I,num_year);
mu_gh_t = zeros(I,num_year);

mu_lg_t = zeros(I,num_year);
mu_ll_t = zeros(I,num_year);
mu_lh_t = zeros(I,num_year);

mu_hg_t = zeros(I,num_year);
mu_hl_t = zeros(I,num_year);
mu_hh_t = zeros(I,num_year);


govaii_data = readtable('wiot_go_va_ii.csv');
io_data = readtable('IO_MX_sectoral.csv');
betaK_data = readtable('betaKj.csv');
ppi_data = readtable('final_price.csv');
emp_data = readtable('emp_sectoral.csv');
cap_data = readtable('k_usd.csv');
pc_px_inv_data = readtable('pc_px_inv.csv');
exp_data = readtable('X_i_sectoral.csv');
imp_data = readtable('M_i_sectoral.csv');
exr_data = readtable('er_i.csv');


% time-invariant IO linkage
count = 1;
for i=1:I
for year = 1:T
        
        exr(i,year) = exr_data.er_i(count);

        alpha_g_t(i,year) = betaK_data.betaKG_i(count);
        alpha_l_t(i,year) = betaK_data.betaKL_i(count);
        alpha_h_t(i,year) = betaK_data.betaKH_i(count);
        
        nu_g_t(i,year) = govaii_data.va_G_usd_i(count)/govaii_data.Y_G_usd_i(count);
        nu_l_t(i,year) = govaii_data.va_L_usd_i(count)/govaii_data.Y_L_usd_i(count);
        nu_h_t(i,year) = govaii_data.va_H_usd_i(count)/govaii_data.Y_H_usd_i(count);
        
        va_g_t(i,year) = govaii_data.va_G_usd_i(count)/(govaii_data.va_G_usd_i(count)+govaii_data.va_L_usd_i(count)+govaii_data.va_H_usd_i(count));
        va_l_t(i,year) = govaii_data.va_L_usd_i(count)/(govaii_data.va_G_usd_i(count)+govaii_data.va_L_usd_i(count)+govaii_data.va_H_usd_i(count));
        va_h_t(i,year) = govaii_data.va_H_usd_i(count)/(govaii_data.va_G_usd_i(count)+govaii_data.va_L_usd_i(count)+govaii_data.va_H_usd_i(count));
        %va_s_t(i,year) = num(9,4)/( num(9,1)+num(9,2)+num(9,3)+num(9,4) );
        
        % mu_ij_t means inputs from j to i (share of good j in total
        % intermediates in sector i
        mu_gg_t(i,year) = io_data.IO_GG(count)/(io_data.IO_GG(count)+io_data.IO_GL(count)+io_data.IO_GH(count));
        mu_gl_t(i,year) = io_data.IO_GL(count)/(io_data.IO_GG(count)+io_data.IO_GL(count)+io_data.IO_GH(count));
        mu_gh_t(i,year) = io_data.IO_GH(count)/(io_data.IO_GG(count)+io_data.IO_GL(count)+io_data.IO_GH(count));
        
        mu_lg_t(i,year) = io_data.IO_LG(count)/(io_data.IO_LG(count)+io_data.IO_LL(count)+io_data.IO_LH(count));
        mu_ll_t(i,year) = io_data.IO_LL(count)/(io_data.IO_LG(count)+io_data.IO_LL(count)+io_data.IO_LH(count));
        mu_lh_t(i,year) = io_data.IO_LH(count)/(io_data.IO_LG(count)+io_data.IO_LL(count)+io_data.IO_LH(count));
        
        mu_hg_t(i,year) = io_data.IO_HG(count)/(io_data.IO_HG(count)+io_data.IO_HL(count)+io_data.IO_HH(count));
        mu_hl_t(i,year) = io_data.IO_HL(count)/(io_data.IO_HG(count)+io_data.IO_HL(count)+io_data.IO_HH(count));
        mu_hh_t(i,year) = io_data.IO_HH(count)/(io_data.IO_HG(count)+io_data.IO_HL(count)+io_data.IO_HH(count));
        
        count = count + 1;
end
end


if Time_varying_IO_table == 0
    
    IO_year = 2; %1 is 1995
    alpha_g_t = repmat(alpha_g_t(:,IO_year),1,num_year);
    alpha_l_t = repmat(alpha_l_t(:,IO_year),1,num_year);
    alpha_h_t = repmat(alpha_h_t(:,IO_year),1,num_year);
    
    nu_g_t = repmat(nu_g_t(:,IO_year),1,num_year);
    nu_l_t = repmat(nu_l_t(:,IO_year),1,num_year);
    nu_h_t = repmat(nu_h_t(:,IO_year),1,num_year);
    
    mu_gg_t = repmat(mu_gg_t(:,IO_year),1,num_year);
    mu_gl_t = repmat(mu_gl_t(:,IO_year),1,num_year);
    mu_gh_t = repmat(mu_gh_t(:,IO_year),1,num_year);
    mu_lg_t = repmat(mu_lg_t(:,IO_year),1,num_year);
    mu_ll_t = repmat(mu_ll_t(:,IO_year),1,num_year);
    mu_lh_t = repmat(mu_lh_t(:,IO_year),1,num_year);
    mu_hg_t = repmat(mu_hg_t(:,IO_year),1,num_year);
    mu_hl_t = repmat(mu_hl_t(:,IO_year),1,num_year);
    mu_hh_t = repmat(mu_hh_t(:,IO_year),1,num_year);
    
elseif  Time_varying_IO_table == 2

    alpha_g_t = repmat(mean(alpha_g_t,2),1,num_year);
    alpha_l_t = repmat(mean(alpha_l_t,2),1,num_year);
    alpha_h_t = repmat(mean(alpha_h_t,2),1,num_year);
    
    nu_g_t = repmat(mean(nu_g_t,2),1,num_year);
    nu_l_t = repmat(mean(nu_l_t,2),1,num_year);
    nu_h_t = repmat(mean(nu_h_t,2),1,num_year);
    
    mu_gg_t = repmat(mean(mu_gg_t,2),1,num_year);
    mu_gl_t = repmat(mean(mu_gl_t,2),1,num_year);
    mu_gh_t = repmat(mean(mu_gh_t,2),1,num_year);
    mu_lg_t = repmat(mean(mu_lg_t,2),1,num_year);
    mu_ll_t = repmat(mean(mu_ll_t,2),1,num_year);
    mu_lh_t = repmat(mean(mu_lh_t,2),1,num_year);
    mu_hg_t = repmat(mean(mu_hg_t,2),1,num_year);
    mu_hl_t = repmat(mean(mu_hl_t,2),1,num_year);
    mu_hh_t = repmat(mean(mu_hh_t,2),1,num_year);        
end

if Same_capital_share == 1
    
    ca_sh = 0.4;
    alpha_g_t = alpha_g_t*0+ca_sh;
    alpha_l_t = alpha_l_t*0+ca_sh;
    alpha_h_t = alpha_h_t*0+ca_sh;
end

if Same_va_share == 1
    
    va_sh = 0.7;
    nu_g_t = nu_g_t*0+va_sh;
    nu_l_t = nu_l_t*0+va_sh;
    nu_h_t = nu_h_t*0+va_sh;
end


if Same_int_demand == 1
    
   int_sh = 0.33;
   mu_gg_t = mu_gg_t*0+int_sh;
   mu_gl_t = mu_gg_t*0+int_sh;
   mu_gh_t = mu_gg_t*0+0.34;
   mu_lg_t = mu_gg_t*0+int_sh;
   mu_ll_t = mu_gg_t*0+int_sh;
   mu_lh_t = mu_gg_t*0+0.34;
   mu_hg_t = mu_gg_t*0+int_sh;
   mu_hl_t = mu_gg_t*0+int_sh;
   mu_hh_t = mu_gg_t*0+0.34;
end    

if Sector_specific_share == 1
    
    alpha_g_t = repmat(mean(alpha_g_t),I,1);
    alpha_l_t = repmat(mean(alpha_l_t),I,1);
    alpha_h_t = repmat(mean(alpha_h_t),I,1);
    
    nu_g_t = repmat(mean(nu_g_t),I,1);
    nu_l_t = repmat(mean(nu_l_t),I,1);
    nu_h_t = repmat(mean(nu_h_t),I,1);
    
    mu_gg_t = repmat(mean(mu_gg_t),I,1);
    mu_gl_t = repmat(mean(mu_gl_t),I,1);
    mu_gh_t = repmat(mean(mu_gh_t),I,1);
    mu_lg_t = repmat(mean(mu_lg_t),I,1);
    mu_ll_t = repmat(mean(mu_ll_t),I,1);
    mu_lh_t = repmat(mean(mu_lh_t),I,1);
    mu_hg_t = repmat(mean(mu_hg_t),I,1);
    mu_hl_t = repmat(mean(mu_hl_t),I,1);
    mu_hh_t = repmat(mean(mu_hh_t),I,1);

    %adjust = 0.021; % prevent negative final demand, and consumption share
    %mu_sa_t = mu_sa_t-adjust;
    %mu_ss_t = mu_ss_t+adjust;
    
end

% time-varying price and allocation
PY_g_t = zeros(I,num_year);
PY_l_t = zeros(I,num_year);
PY_h_t = zeros(I,num_year);

P_g_t = zeros(I,num_year);
P_l_t = zeros(I,num_year);
P_h_t = zeros(I,num_year);

l_g_t = zeros(I,num_year);
l_l_t = zeros(I,num_year);
l_h_t = zeros(I,num_year);

k_g_t = zeros(I,num_year);
k_l_t = zeros(I,num_year);
k_h_t = zeros(I,num_year);

PC_g_t_raw = zeros(I,num_year);
PC_l_t_raw = zeros(I,num_year);
PC_h_t_raw = zeros(I,num_year);

PX_g_t_raw = zeros(I,num_year);
PX_l_t_raw = zeros(I,num_year);
PX_h_t_raw = zeros(I,num_year);

EX_g_t = zeros(I,num_year);
EX_l_t = zeros(I,num_year);
EX_h_t = zeros(I,num_year);

IM_g_t = zeros(I,num_year);
IM_l_t = zeros(I,num_year);
IM_h_t = zeros(I,num_year);


count = 1;
for i=1:I
    for year=1:num_year
        % sector nominal output
        PY_g_t(i,year) = govaii_data.Y_G_usd_i(count);
        PY_l_t(i,year) = govaii_data.Y_L_usd_i(count);
        PY_h_t(i,year) = govaii_data.Y_H_usd_i(count);
    
        
        %this is refined data (p_sea_t/p_sea_2005)xp_ggdc_2005
        P_g_t(i,year) = ppi_data.p_G(count);
        P_l_t(i,year) = ppi_data.p_L(count);
        P_h_t(i,year) = ppi_data.p_H(count);
        
        % sector labor
        l_g_t(i,year) = emp_data.empG(count);
        l_l_t(i,year) = emp_data.empL(count);
        l_h_t(i,year) = emp_data.empH(count);

        l_g_t(i,year) = emp_data.empG(count)*1e3; % change from per thousand to each worker
        l_l_t(i,year) = emp_data.empL(count)*1e3;
        l_h_t(i,year) = emp_data.empH(count)*1e3;
        
        % sector capital
        k_g_t(i,year) = cap_data.k_G_usd(count);
        k_l_t(i,year) = cap_data.k_L_usd(count);
        k_h_t(i,year) = cap_data.k_H_usd(count);
        
        % sector final consumption expenditure
        PC_g_t_raw(i,year) = pc_px_inv_data.PC_G(count);
        PC_l_t_raw(i,year) = pc_px_inv_data.PC_L(count);
        PC_h_t_raw(i,year) = pc_px_inv_data.PC_H(count);
        
        % sector final investment expenditure (GFCF+inventory)
        PX_g_t_raw(i,year) = pc_px_inv_data.PX_G(count)+pc_px_inv_data.INV_G(count);
        PX_l_t_raw(i,year) = pc_px_inv_data.PX_L(count)+pc_px_inv_data.INV_L(count);
        PX_h_t_raw(i,year) = pc_px_inv_data.PX_H(count)+pc_px_inv_data.INV_H(count);     
        
        % sector export
        EX_g_t(i,year) = exp_data.X_i_G(count);
        EX_l_t(i,year) = exp_data.X_i_L(count);
        EX_h_t(i,year) = exp_data.X_i_H(count);
        
        % sector import
        IM_g_t(i,year) = imp_data.M_i_G(count);
        IM_l_t(i,year) = imp_data.M_i_L(count);
        IM_h_t(i,year) = imp_data.M_i_H(count); 
        count = count + 1;
    end
end

% time-varying bilteral trade
bt_data=readtable('bilateral_trade.csv');
bt_t=zeros(3,I,I,num_year); % index for destination (importer) first

count = 1;
for i=1:I % index for importer
    for j=1:I % index for exporter
        for year=1:num_year
            for sec=1:3 % index for sector
                bt_t(sec,i,j,year) = bt_data.bt_all(count);% all bt
                %bt_t(sec,i,j,year) = num((sec-1)*8+(i-1)*128+(j-1)*32+year,5)-num((sec-1)*8+(i-1)*128+(j-1)*32+year,6); % int bt
                count = count + 1;
            end
        end
    end
end


bt_g_t=squeeze(bt_t(1,:,:,:));
bt_l_t=squeeze(bt_t(2,:,:,:));
bt_h_t=squeeze(bt_t(3,:,:,:));

%bt_h_t(bt_h_t<1) = 1;
 
btl_g_t = bt_g_t; % bilateral trade level
btl_l_t = bt_l_t; % bilateral trade level
btl_h_t = bt_h_t; % bilateral trade level

bt_g_t = bt_g_t./repmat(sum(bt_g_t,2),1,I,1); % convert to trade share
bt_l_t = bt_l_t./repmat(sum(bt_l_t,2),1,I,1);
bt_h_t = bt_h_t./repmat(sum(bt_h_t,2),1,I,1);

% sector investment share to produce the investment good
PX_t_raw = PX_g_t_raw+PX_l_t_raw+PX_h_t_raw;
omega_xg_t_raw = PX_g_t_raw./PX_t_raw;
omega_xl_t_raw = PX_l_t_raw./PX_t_raw;
omega_xh_t_raw = PX_h_t_raw./PX_t_raw;

PC_t_raw = PC_g_t_raw+PC_l_t_raw+PC_h_t_raw;
omega_cg_t_raw = PC_g_t_raw./PC_t_raw;
omega_cl_t_raw = PC_l_t_raw./PC_t_raw;
omega_ch_t_raw = PC_h_t_raw./PC_t_raw;


omega_xg_t = omega_xg_t_raw; % match total and sector investment share and level
omega_xl_t = omega_xl_t_raw;
omega_xh_t = omega_xh_t_raw;

PX_t = PX_t_raw; % match total investment, and sector consumption share
PX_g_t = PX_t.*omega_xg_t;
PX_l_t = PX_t.*omega_xl_t;
PX_h_t = PX_t.*omega_xh_t;


% sector labor and capital share from the data (capital share is from nominal capital)
lh_g_t = l_g_t./(l_g_t+l_l_t+l_h_t);
lh_l_t = l_l_t./(l_g_t+l_l_t+l_h_t);
lh_h_t = l_h_t./(l_g_t+l_l_t+l_h_t);

kh_g_t = k_g_t./(k_g_t+k_l_t+k_h_t);
kh_l_t = k_l_t./(k_g_t+k_l_t+k_h_t);
kh_h_t = k_h_t./(k_g_t+k_l_t+k_h_t);


% real capital stock from WIOD 1995-2011
%[num,txt,raw]=xlsread('Capital_notes','SEA');
real_k_data = readtable("real_k_sectoral.csv");
real_K_g_t = zeros(3,2011-1995+1);
real_K_l_t = zeros(3,2011-1995+1);
real_K_h_t = zeros(3,2011-1995+1);

real_K_g_t_2 = zeros(3,2011-1995+1);
real_K_l_t_2 = zeros(3,2011-1995+1);
real_K_h_t_2 = zeros(3,2011-1995+1);

count = 1;
for i=1:I
    for year=1:2011-1995+1
        % sector real capital stock 1995 price
        real_K_g_t(i,year) = real_k_data.real_k_G(count);
        real_K_l_t(i,year) = real_k_data.real_k_L(count);
        real_K_h_t(i,year) = real_k_data.real_k_H(count);
        count = count + 1;
    end
end 
real_K_t_1 = real_K_g_t+real_K_l_t+real_K_h_t;  

real_K_t = real_K_t_1;
real_K_t = real_K_t(:,1:num_year); % year 1995 to 2014


figure('Position', get(0, 'Screensize'))
subplot(2,1,1)
plot(1995:2011,real_K_t(1,:)/real_K_t(5,1),'-b',1995:2011,real_K_t(2,:)/real_K_t(5,1),'-r',1995:2011,real_K_t(3,:)/real_K_t(5,1),'-g',1995:2011,real_K_t(4,:)/real_K_t(5,1),'-c',1995:2011,real_K_t(5,:)/real_K_t(5,1),'-k',...
    1995:2011,real_K_t(7,:)/real_K_t(5,1),'-m',1995:2011,real_K_t(9,:)/real_K_t(5,1),'-y','Linewidth',2.5)
title('K GFCF relative to IND 1995')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','Location','best')
grid on

subplot(2,1,2)
plot(1995:2011,real_K_t(1,:)/real_K_t(1,1),'-b',1995:2011,real_K_t(2,:)/real_K_t(2,1),'-r',1995:2011,real_K_t(3,:)/real_K_t(3,1),'-g',1995:2011,real_K_t(4,:)/real_K_t(4,1),'-c',1995:2011,real_K_t(5,:)/real_K_t(5,1),'-k',...
    1995:2011,real_K_t(7,:)/real_K_t(7,1),'-m',1995:2011,real_K_t(9,:)/real_K_t(9,1),'-y','Linewidth',2.5)
title('K GFCF relative to 1995')
grid on

sgtitle('Figure 1. Real capital stock 1995-2011 USD')
saveas(gcf, 'plots_transition/01_Real_capital_stock.png');
close(gcf)


%% Back out the shocks in the level

% back out model consistent final demand F=PC+PX, using the IO table
F_g_t = zeros(I,num_year);
F_l_t = zeros(I,num_year);
F_h_t = zeros(I,num_year);

PQ_g_t = zeros(I,num_year);
PQ_l_t = zeros(I,num_year);
PQ_h_t = zeros(I,num_year);

for year=1:num_year
    
    pi_g = squeeze(bt_g_t(:,:,year));
    pi_l = squeeze(bt_l_t(:,:,year));
    pi_h = squeeze(bt_h_t(:,:,year));
    
    [F_g_t(:,year),F_l_t(:,year),F_h_t(:,year),PQ_g_t(:,year),PQ_l_t(:,year),PQ_h_t(:,year)]...
        = back_out_final_demand(PY_g_t(:,year),PY_l_t(:,year),PY_h_t(:,year),...
        pi_g,pi_l,pi_h,I,nu_g_t(:,year),nu_l_t(:,year),nu_h_t(:,year),...
        mu_gg_t(:,year),mu_gl_t(:,year),mu_gh_t(:,year),...
        mu_lg_t(:,year),mu_ll_t(:,year),mu_lh_t(:,year),...
        mu_hg_t(:,year),mu_hl_t(:,year),mu_hh_t(:,year));
    
end

F_t = F_g_t+F_l_t+F_h_t;
F_g_t
F_l_t
F_h_t

% --- Feasibility: ensure PX_t <= F_t elementwise (keeps PC_t >= 0)
scale_fac = min(1, F_t ./ max(PX_t, eps));  % eps avoids divide-by-zero
PX_g_t = PX_g_t .* scale_fac;
PX_l_t = PX_l_t .* scale_fac;
PX_h_t = PX_h_t .* scale_fac;
PX_t   = PX_g_t + PX_l_t + PX_h_t;          % refresh PX_t after scaling


% check the consistency for model implied sector absorption
PQ_g_t_raw = PY_g_t+IM_g_t-EX_g_t;
PQ_l_t_raw = PY_l_t+IM_l_t-EX_l_t;
PQ_h_t_raw = PY_h_t+IM_h_t-EX_h_t;


PQ_g_t_raw-PQ_g_t % PY matched, pai matched, then PQ is matched
PQ_l_t_raw-PQ_l_t
PQ_h_t_raw-PQ_h_t


% back out model consistent NX,GDP, WGDP(numeraire) and WPC 
NX_g_t = PY_g_t-PQ_g_t;
NX_l_t = PY_l_t-PQ_l_t;
NX_h_t = PY_h_t-PQ_h_t;
NX_t = NX_g_t+NX_l_t+NX_h_t;

GDP_t =  PY_g_t.*nu_g_t+PY_l_t.*nu_l_t+PY_h_t.*nu_h_t;
WGDP_t = sum(GDP_t);


% --- Author-style consumption construction: PC_total from F-PX_total,
%     then allocate across sectors using *consumption shares from data*.
PC_t = max(F_t - PX_t, 0)  % guarantees nonneg total consumption

% Start from raw consumption shares; clip/renormalize row-wise for safety.
omega_cg_t = omega_cg_t_raw;
omega_cl_t = omega_cl_t_raw;
omega_ch_t = omega_ch_t_raw;

% Allocate sectoral consumption from PC_total
PC_g_t = PC_t .* omega_cg_t;
PC_l_t = PC_t .* omega_cl_t;
PC_h_t = PC_t .* omega_ch_t;

PC_g_t


check_valid_omega = [min(omega_cg_t(:)),min(omega_cl_t(:)),min(omega_ch_t(:))]
check_valid_omega = [min(omega_xg_t(:)),min(omega_xl_t(:)),min(omega_xh_t(:))]

P_c_t = (P_g_t./omega_cg_t).^omega_cg_t.*...
    (P_l_t./omega_cl_t).^omega_cl_t.*...
    (P_h_t./omega_ch_t).^omega_ch_t;

P_x_t = (P_g_t./omega_xg_t).^omega_xg_t.*...
    (P_l_t./omega_xl_t).^omega_xl_t.*...
    (P_h_t./omega_xh_t).^omega_xh_t;


% Price by sector and country relative to 1995
figure('Position', get(0, 'Screensize'))
subplot(2,2,1)
plot(1995:2011,P_g_t(1,:)/P_g_t(1,1),'-b',1995:2011,P_g_t(2,:)/P_g_t(2,1),'-r',1995:2011,P_g_t(3,:)/P_g_t(3,1),'-g',1995:2011,P_g_t(4,:)/P_g_t(4,1),'-c',1995:2011,P_g_t(5,:)/P_g_t(5,1),'-k',...
    1995:2011,P_g_t(7,:)/P_g_t(7,1),'-m',1995:2011,P_g_t(9,:)/P_g_t(9,1),'-y','Linewidth',2.5)
title('Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','Location','best')
grid on

subplot(2,2,2)
plot(1995:2011,P_l_t(1,:)/P_l_t(1,1),'-b',1995:2011,P_l_t(2,:)/P_l_t(2,1),'-r',1995:2011,P_l_t(3,:)/P_l_t(3,1),'-g',1995:2011,P_l_t(4,:)/P_l_t(4,1),'-c',1995:2011,P_l_t(5,:)/P_l_t(5,1),'-k',...
    1995:2011,P_l_t(7,:)/P_l_t(7,1),'-m',1995:2011,P_l_t(9,:)/P_l_t(9,1),'-y','Linewidth',2.5)
title('Low-Skill Service')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','Location','best')
grid on

subplot(2,2,3)
plot(1995:2011,P_h_t(1,:)/P_h_t(1,1),'-b',1995:2011,P_h_t(2,:)/P_h_t(2,1),'-r',1995:2011,P_h_t(3,:)/P_h_t(3,1),'-g',1995:2011,P_h_t(4,:)/P_h_t(4,1),'-c',1995:2011,P_h_t(5,:)/P_h_t(5,1),'-k',...
    1995:2011,P_h_t(7,:)/P_h_t(7,1),'-m',1995:2011,P_h_t(9,:)/P_h_t(9,1),'-y','Linewidth',2.5)
title('High-Skill Service')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','Location','best')
grid on

sgtitle('Figure 2. Prices by sector and country relative to 1995 (raw data)')
saveas(gcf, 'plots_transition/02_Prices_raw.png');
close(gcf)

% exr for BRA, CHN, IND
figure('Position', get(0, 'Screensize'))
subplot(2,2,1)
plot(1995:2011,exr(2,:),'-r','Linewidth',2.5)
title('BRA')
grid on

subplot(2,2,2)
plot(1995:2011,exr(3,:),'-g','Linewidth',2.5)
title('CHN')
grid on

subplot(2,2,3)
plot(1995:2011,exr(5,:),'-k','Linewidth',2.5)
title('IND')
grid on

subplot(2,2,4)
plot(1995:2011,exr(7,:),'-m','Linewidth',2.5)
title('MEX')
grid on

sgtitle('Figure 3. Exchange rate used in WIOD')
saveas(gcf, 'plots_transition/03_exr.png');
close(gcf)

% pin down the price level, based on restricting RGDP over investment good
% equals the real GDP moments, that is: (GDP_US_2005/WGDP_2005)/P_c_US_2005 = GDP_US_2005 
% or P_c_US_2005 = 1/WGDP_2005, GDP_US_2005 = 1.3321*1e7 million
numeraire = WGDP_t;
conversion = 1e1/P_c_t(9,11)./repmat(numeraire,I,1); % captures the GDP price of US in 2005
%conversion = 1e9/P_c_t(9,11)./repmat(numeraire,I,1); % check whether the level matters

P_g_t = P_g_t.*conversion;
P_l_t = P_l_t.*conversion;
P_h_t = P_h_t.*conversion;

P_c_t = P_c_t.*conversion;
P_x_t = P_x_t.*conversion;

% convert the relevant nominal variables (except for sector prices) in terms of the numeraire WGDP
PY_g_t_raw = PY_g_t;
PY_l_t_raw = PY_l_t;
PY_h_t_raw = PY_h_t;

PY_g_t = PY_g_t./repmat(numeraire,I,1);
PY_l_t = PY_l_t./repmat(numeraire,I,1);
PY_h_t = PY_h_t./repmat(numeraire,I,1);

PQ_g_t = PQ_g_t./repmat(numeraire,I,1);
PQ_l_t = PQ_l_t./repmat(numeraire,I,1);
PQ_h_t = PQ_h_t./repmat(numeraire,I,1);

NX_g_t = NX_g_t./repmat(numeraire,I,1);
NX_l_t = NX_l_t./repmat(numeraire,I,1);
NX_h_t = NX_h_t./repmat(numeraire,I,1);

F_g_t = F_g_t./repmat(numeraire,I,1);
F_l_t = F_l_t./repmat(numeraire,I,1);
F_h_t = F_h_t./repmat(numeraire,I,1);

F_t = F_t./repmat(numeraire,I,1);

GDP_t = GDP_t./repmat(numeraire,I,1);
NX_t = NX_t./repmat(numeraire,I,1);

PC_t = PC_t./repmat(numeraire,I,1);
PC_g_t = PC_g_t./repmat(numeraire,I,1);
PC_l_t = PC_l_t./repmat(numeraire,I,1);
PC_h_t = PC_h_t./repmat(numeraire,I,1);

PX_t = PX_t./repmat(numeraire,I,1);
PX_g_t = PX_g_t./repmat(numeraire,I,1);
PX_l_t = PX_l_t./repmat(numeraire,I,1);
PX_h_t = PX_h_t./repmat(numeraire,I,1);

rho_t = PX_t./GDP_t;

GDP_t(9,11)/P_c_t(9,11)

% back out the bilateral trade cost
% spike for trade cost for import from BRA is affected by decreasing trade share on BRA service, seems to relativ with BRA depreciation
% so how does WIOD construct the bt level in US, directly in US or convert from local currency
% price also increases relative to BRA in service, due to BRA depreciation
% and lower price. Maybe it is fine, as trade cost level is high in service
tao_g_t = bt_g_t*0;
tao_l_t = bt_l_t*0;
tao_h_t = bt_h_t*0;

for year=1:num_year
    for i=1:I % index for importer
        for j=1:I % index for exporter
            tao_g_t(i,j,year) = (bt_g_t(i,j,year)/bt_g_t(j,j,year))^(-1/theta)*P_g_t(i,year)/P_g_t(j,year);
            tao_l_t(i,j,year) = (bt_l_t(i,j,year)/bt_l_t(j,j,year))^(-1/theta)*P_l_t(i,year)/P_l_t(j,year);
            tao_h_t(i,j,year) = (bt_h_t(i,j,year)/bt_h_t(j,j,year))^(-1/theta)*P_h_t(i,year)/P_h_t(j,year);         
        end
    end
end

[min(tao_g_t(:)) min(tao_l_t(:)) min(tao_h_t(:))]

% trade cost for Indian
figure('Position', get(0, 'Screensize'))
subplot(2,3,1)
plot(1995:2011,squeeze(tao_g_t(1,5,:))/squeeze(tao_g_t(1,5,1)),'-b',1995:2011,squeeze(tao_g_t(2,5,:))/squeeze(tao_g_t(2,5,1)),'-r',1995:2011,squeeze(tao_g_t(3,5,:))/squeeze(tao_g_t(3,5,1)),'-g',1995:2011,squeeze(tao_g_t(4,5,:))/squeeze(tao_g_t(4,5,1)),'-c',...
    1995:2011,squeeze(tao_g_t(7,5,:))/squeeze(tao_g_t(7,5,1)),'-m',1995:2011,squeeze(tao_g_t(9,5,:))/squeeze(tao_g_t(9,5,1)),'-y',1995:2011,squeeze(tao_g_t(8,5,:))/squeeze(tao_g_t(8,5,1)),'-k',1995:2011,ones(1,num_year),'--black','Linewidth',2.5)
title('Goods exports')
legend('To AUS','To BRA','To CHN','To DEU','To MEX','To USA','To ROW','value of 1','location','best')
grid on

subplot(2,3,2)
plot(1995:2011,squeeze(tao_l_t(1,5,:))/squeeze(tao_l_t(1,5,1)),'-b',1995:2011,squeeze(tao_l_t(2,5,:))/squeeze(tao_l_t(2,5,1)),'-r',1995:2011,squeeze(tao_l_t(3,5,:))/squeeze(tao_l_t(3,5,1)),'-g',1995:2011,squeeze(tao_l_t(4,5,:))/squeeze(tao_l_t(4,5,1)),'-c',...
    1995:2011,squeeze(tao_l_t(7,5,:))/squeeze(tao_l_t(7,5,1)),'-m',1995:2011,squeeze(tao_l_t(9,5,:))/squeeze(tao_l_t(9,5,1)),'-y',1995:2011,squeeze(tao_l_t(8,5,:))/squeeze(tao_l_t(8,5,1)),'-k',1995:2011,ones(1,num_year),'--black','Linewidth',2.5)
title('Low-skill service exports')
grid on

subplot(2,3,3)
plot(1995:2011,squeeze(tao_h_t(1,5,:))/squeeze(tao_h_t(1,5,1)),'-b',1995:2011,squeeze(tao_h_t(2,5,:))/squeeze(tao_h_t(2,5,1)),'-r',1995:2011,squeeze(tao_h_t(3,5,:))/squeeze(tao_h_t(3,5,1)),'-g',1995:2011,squeeze(tao_h_t(4,5,:))/squeeze(tao_h_t(4,5,1)),'-c',...
    1995:2011,squeeze(tao_h_t(7,5,:))/squeeze(tao_h_t(7,5,1)),'-m',1995:2011,squeeze(tao_h_t(9,5,:))/squeeze(tao_h_t(9,5,1)),'-y',1995:2011,squeeze(tao_h_t(8,5,:))/squeeze(tao_h_t(8,5,1)),'-k',1995:2011,ones(1,num_year),'--black','Linewidth',2.5)
title('High-skill service exports')
grid on

subplot(2,3,4)
plot(1995:2011,squeeze(tao_g_t(5,1,:))/squeeze(tao_g_t(5,1,1)),'-b',1995:2011,squeeze(tao_g_t(5,2,:))/squeeze(tao_g_t(5,2,1)),'-r',1995:2011,squeeze(tao_g_t(5,3,:))/squeeze(tao_g_t(5,3,1)),'-g',1995:2011,squeeze(tao_g_t(5,4,:))/squeeze(tao_g_t(5,4,1)),'-c',...
    1995:2011,squeeze(tao_g_t(5,7,:))/squeeze(tao_g_t(5,7,1)),'-m',1995:2011,squeeze(tao_g_t(5,9,:))/squeeze(tao_g_t(5,9,1)),'-y',1995:2011,squeeze(tao_g_t(5,8,:))/squeeze(tao_g_t(5,8,1)),'-k',1995:2011,ones(1,num_year),'--black','Linewidth',2.5)
title('Goods imports')
legend('From AUS','From BRA','From CHN','From DEU','From MEX','From USA','From ROW','value of 1','location','best')
grid on

subplot(2,3,5)
plot(1995:2011,squeeze(tao_l_t(5,1,:))/squeeze(tao_l_t(5,1,1)),'-b',1995:2011,squeeze(tao_l_t(5,2,:))/squeeze(tao_l_t(5,2,1)),'-r',1995:2011,squeeze(tao_l_t(5,3,:))/squeeze(tao_l_t(5,3,1)),'-g',1995:2011,squeeze(tao_l_t(5,4,:))/squeeze(tao_l_t(5,4,1)),'-c',...
    1995:2011,squeeze(tao_l_t(5,7,:))/squeeze(tao_l_t(5,7,1)),'-m',1995:2011,squeeze(tao_l_t(5,9,:))/squeeze(tao_l_t(5,9,1)),'-y',1995:2011,squeeze(tao_l_t(5,8,:))/squeeze(tao_l_t(5,8,1)),'-k',1995:2011,ones(1,num_year),'--black','Linewidth',2.5)
title('Low-skill service imports')
grid on

subplot(2,3,6)
plot(1995:2011,squeeze(tao_h_t(5,1,:))/squeeze(tao_h_t(5,1,1)),'-b',1995:2011,squeeze(tao_h_t(5,2,:))/squeeze(tao_h_t(5,2,1)),'-r',1995:2011,squeeze(tao_h_t(5,3,:))/squeeze(tao_h_t(5,3,1)),'-g',1995:2011,squeeze(tao_h_t(5,4,:))/squeeze(tao_h_t(5,4,1)),'-c',...
    1995:2011,squeeze(tao_h_t(5,7,:))/squeeze(tao_h_t(5,7,1)),'-m',1995:2011,squeeze(tao_h_t(5,9,:))/squeeze(tao_h_t(5,9,1)),'-y',1995:2011,squeeze(tao_h_t(5,8,:))/squeeze(tao_h_t(5,8,1)),'-k',1995:2011,ones(1,num_year),'--black','Linewidth',2.5)
title('High-skill service imports')
grid on

%sgtitle('Figure. Bilateral trade costs for Indian') % also reflect exchange rate change
sgtitle('Bilateral trade costs for India relative to 1995') % also reflect exchange rate change
saveas(gcf, 'plots_transition/04_Trade_cost_Indian.png');
%saveas(gcf, 'plots_transition/Trade_cost_Indian_trade_share.png');
%saveas(gcf, 'plots_transition/Trade_cost_Indian_trade_price.png');
close(gcf)

tao_g_t(tao_g_t<1)=1;
tao_l_t(tao_l_t<1)=1;
tao_h_t(tao_h_t<1)=1;

% tao_g_t(tao_g_t==Inf)=5000;
% tao_l_t(tao_l_t==Inf)=5000;
% tao_h_t(tao_h_t==Inf)=5000;


[min(tao_g_t(:)) min(tao_l_t(:)) min(tao_h_t(:))]
[max(tao_g_t(:)) max(tao_l_t(:)) max(tao_h_t(:))]

% labor quantity 
l_t = l_g_t+l_l_t+l_h_t; % in thousand people
%l_t(3,:) = l_t(3,:)*2.5;

% use fi to denote share of income sent to global portfolio
fi_t = zeros(I,num_year);
diff_t = fi_t*0;
for i=1:num_year
    [fi_t(:,i),diff_t(:,i)] = global_portfolio(I,GDP_t(:,i),NX_t(:,i),l_t(:,i));
end    

portfolio_t = sum(fi_t.*GDP_t);
%transfer_t = portfolio_t.*l_t./sum(l_t);
transfer_t = portfolio_t.*GDP_t./sum(GDP_t);
%transfer_t = portfolio_t/I;
NX_portfolio_t = fi_t.*GDP_t-transfer_t;

figure('Position', get(0, 'Screensize')) 
subplot(2,1,1) % IND relative low GDP, but high pop, easy to have negative net export, which result in rou=1 for early periods
plot(1995:2011,(NX_portfolio_t(1,:)./GDP_t(1,:))/(NX_portfolio_t(1,1)./GDP_t(1,1)),'-b',1995:2011,(NX_portfolio_t(2,:)./GDP_t(2,:))/(NX_portfolio_t(2,1)./GDP_t(2,1)),'-r',1995:2011,(NX_portfolio_t(3,:)./GDP_t(3,:))/(NX_portfolio_t(3,1)./GDP_t(3,1)),'-g',1995:2011,(NX_portfolio_t(4,:)./GDP_t(4,:))/(NX_portfolio_t(4,1)./GDP_t(4,1)),'-c',1995:2011,(NX_portfolio_t(5,:)./GDP_t(5,:))/(NX_portfolio_t(5,1)./GDP_t(5,1)),'-k',1995:2011,(NX_portfolio_t(7,:)./GDP_t(7,:))/(NX_portfolio_t(7,1)./GDP_t(7,1)),'-m',1995:2011,(NX_portfolio_t(9,:)./GDP_t(9,:))/(NX_portfolio_t(9,1)./GDP_t(9,1)),'-y',...
    1995:2011,(NX_t(1,:)./GDP_t(1,:))/(NX_t(1,1)./GDP_t(1,1)),'--b',1995:2011,(NX_t(2,:)./GDP_t(2,:))/(NX_t(2,1)./GDP_t(2,1)),'--r',1995:2011,(NX_t(3,:)./GDP_t(3,:))/(NX_t(3,1)./GDP_t(3,1)),'--g',1995:2011,(NX_t(4,:)./GDP_t(4,:))/(NX_t(4,1)./GDP_t(4,1)),'--c',1995:2011,(NX_t(5,:)./GDP_t(5,:))/(NX_t(5,1)./GDP_t(5,1)),'--k',1995:2011,(NX_t(7,:)./GDP_t(7,:))/(NX_t(7,1)./GDP_t(7,1)),'--m',1995:2011,(NX_t(9,:)./GDP_t(9,:))/(NX_t(9,1)./GDP_t(9,1)),'--y',...
    'Linewidth',2.5)
title('Net export over GDP')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','Raw Data','location','best')
grid on

subplot(2,1,2)
plot(1995:2011,fi_t(1,:)/fi_t(1,1),'-b',1995:2011,fi_t(2,:)/fi_t(2,1),'-r',1995:2011,fi_t(3,:)/fi_t(3,1),'-g',1995:2011,fi_t(4,:)/fi_t(4,1),'-c',1995:2011,fi_t(5,:)/fi_t(5,1),'-k',1995:2011,fi_t(7,:)/fi_t(7,1),'-m',1995:2011,fi_t(9,:)/fi_t(9,1),'-y','Linewidth',2.5)
title('Portfolio investment share')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','northeast') 
grid on

sgtitle('Figure. Static global portfolio relative to 1995')
saveas(gcf, 'plots_transition\05_Target_NX_over_GDP.png');
close(gcf)


% Employment
figure('Position', get(0, 'Screensize'))
subplot(2,1,1)
plot(1995:2011,l_t(1,:),'-b',1995:2011,l_t(2,:),'-r',1995:2011,l_t(3,:),'-g',1995:2011,l_t(4,:),'-c',1995:2011,l_t(5,:),'-k',1995:2011,l_t(7,:),'-m',1995:2011,l_t(9,:),'-y','Linewidth',2.5)
title('Employment (in thousands of people)')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','northeast') 
grid on

subplot(2,1,2)
plot(1995:2011,l_t(1,:)/l_t(1,1),'-b',1995:2011,l_t(2,:)/l_t(2,1),'-r',1995:2011,l_t(3,:)/l_t(3,1),'-g',1995:2011,l_t(4,:)/l_t(4,1),'-c',1995:2011,l_t(5,:)/l_t(5,1),'-k',1995:2011,l_t(7,:)/l_t(7,1),'-m',1995:2011,l_t(9,:)/l_t(9,1),'-y','Linewidth',2.5)
title('Employment relative to 1995')
grid on

sgtitle('Figure. Employment')
saveas(gcf, 'plots_transition/06_Employment.png');
close(gcf)



%% Given the investment efficiency, back out the capital stocks, and sector productivity

A_x_guess = zeros(I,1);

A_x_guess(1) = 1.5;  % AUS
A_x_guess(2) = 0.24;  % BRA
A_x_guess(3) = 3.5;  % CHN
A_x_guess(4) = 0.75;  % DEU
A_x_guess(5) = 1.3;  % IND
A_x_guess(6) = 10.5;  % JPN
A_x_guess(7) = 0.62;  % MEX
A_x_guess(8) = 3.9;  % ROW
A_x_guess(9) = 0.51;  % USA


k_t = real_K_t; % directly input capital quantity from WIOD 1995 to 2014
X_t = PX_t./P_x_t;

K_t = repmat(k_t(:,1),1,num_year+1); % input K1995
K_t(:,2) = (1-del).*K_t(:,1)+A_x_guess.*( X_t(:,1).^lambda.*K_t(:,1).^(1-lambda) ); 
r_t(:,1) = ( alpha_g_t(:,1).*PY_g_t(:,1).*nu_g_t(:,1)+alpha_l_t(:,1).*PY_l_t(:,1).*nu_l_t(:,1)+...
    alpha_h_t(:,1).*PY_h_t(:,1).*nu_h_t(:,1) )./K_t(:,1);
for i=1:num_year-1
    
    r_t(:,i+1) = ( alpha_g_t(:,i+1).*PY_g_t(:,i+1).*nu_g_t(:,i+1)+alpha_l_t(:,i+1).*PY_l_t(:,i+1).*nu_l_t(:,i+1)+...
        alpha_h_t(:,i+1).*PY_h_t(:,i+1).*nu_h_t(:,i+1) )./K_t(:,i+1);
        
    LHS = PC_t(:,i+1)./PC_t(:,i).*PX_t(:,i)./( K_t(:,i+1)-(1-del).*K_t(:,i) );
    temp = LHS/beta./lambda-r_t(:,i+1)-(1-lambda)./lambda.*PX_t(:,i+1)./K_t(:,i+1);
    K_t(:,i+2) = (1-del)./lambda.*PX_t(:,i+1)./temp+(1-del).*K_t(:,i+1);
end
Ax_t = ( K_t(:,2:num_year+1)-(1-del).*K_t(:,1:num_year) )./(X_t.^lambda.*K_t(:,1:num_year).^(1-lambda));
P_k_t = 1./lambda./Ax_t.*(X_t./K_t(:,1:num_year)).^(1-lambda).*P_x_t;
PK_t = P_k_t.*K_t(:,1:num_year);

L_t = l_t; % per thousand people
w_t = ( (1-alpha_g_t).*PY_g_t.*nu_g_t+(1-alpha_l_t).*PY_l_t.*nu_l_t+...
    (1-alpha_h_t).*PY_h_t.*nu_h_t )./L_t;
w_t_raw = ( (1-alpha_g_t).*PY_g_t_raw.*nu_g_t+(1-alpha_l_t).*PY_l_t_raw.*nu_l_t+...
    (1-alpha_h_t).*PY_h_t_raw.*nu_h_t )./L_t;
GDP_t_raw =  PY_g_t_raw.*nu_g_t+PY_l_t_raw.*nu_l_t+PY_h_t_raw.*nu_h_t;
Ax_t_raw = ( k_t(:,2:num_year)-(1-del).*k_t(:,1:num_year-1) )./(X_t(:,1:num_year-1).^lambda.*k_t(:,1:num_year-1).^(1-lambda));
%Ax_t_PWT_raw = ( real_K_t_simul_PWT(:,2:8)-(1-del).*real_K_t_simul_PWT(:,1:8-1) )./(X_t(:,1:8-1).^lambda.*real_K_t_simul_PWT(:,1:8-1).^(1-lambda));


disp('seq of capital stock')
K_t
disp('seq of Ax')
Ax_t



% undertand why countries differ in Ax
figure('Position', get(0, 'Screensize'))
subplot(2,3,1)
plot(1995:2011,K_t(1,1:num_year),'-b',1995:2011,K_t(2,1:num_year),'-r',1995:2011,K_t(3,1:num_year),'-g',1995:2011,K_t(4,1:num_year),'-c',1995:2011,K_t(5,1:num_year),'-k',1995:2011,K_t(7,1:num_year),'-m',1995:2011,K_t(9,1:num_year),'-y','Linewidth',2.5)
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
title('Capital stock K')
grid on

subplot(2,3,2)
plot(1995:2011,K_t(1,1:num_year)/K_t(1,1),'-b',1995:2011,K_t(2,1:num_year)/K_t(2,1),'-r',1995:2011,K_t(3,1:num_year)/K_t(3,1),'-g',1995:2011,K_t(4,1:num_year)/K_t(4,1),'-c',1995:2011,K_t(5,1:num_year)/K_t(5,1),'-k',1995:2011,K_t(7,1:num_year)/K_t(7,1),'-m',1995:2011,K_t(9,1:num_year)/K_t(9,1),'-y','Linewidth',2.5)
title('Capital stock K over K1995')
grid on

subplot(2,3,3)
plot(1995:2011,P_k_t(1,1:num_year),'-b',1995:2011,P_k_t(2,1:num_year),'-r',1995:2011,P_k_t(3,1:num_year),'-g',1995:2011,P_k_t(4,1:num_year),'-c',1995:2011,P_k_t(5,1:num_year),'-k',1995:2011,P_k_t(7,1:num_year),'-m',1995:2011,P_k_t(9,1:num_year),'-y','Linewidth',2.5)
title('Price of new capital in WGDP')
grid on

subplot(2,3,4)
plot(1995:2011,Ax_t(1,1:num_year),'-b',1995:2011,Ax_t(2,1:num_year),'-r',1995:2011,Ax_t(3,1:num_year),'-g',1995:2011,Ax_t(4,1:num_year),'-c',1995:2011,Ax_t(5,1:num_year),'-k',1995:2011,Ax_t(7,1:num_year),'-m',1995:2011,Ax_t(9,1:num_year),'-y','Linewidth',2.5)
title('Investment efficiency')
grid on

subplot(2,3,5)
plot(1995:2011,Ax_t(1,1:num_year)/Ax_t(1,1),'-b',1995:2011,Ax_t(2,1:num_year)/Ax_t(2,1),'-r',1995:2011,Ax_t(3,1:num_year)/Ax_t(3,1),'-g',1995:2011,Ax_t(4,1:num_year)/Ax_t(4,1),'-c',1995:2011,Ax_t(5,1:num_year)/Ax_t(5,1),'-k',1995:2011,Ax_t(7,1:num_year)/Ax_t(7,1),'-m',1995:2011,Ax_t(9,1:num_year)/Ax_t(9,1),'-y','Linewidth',2.5)
title('Investment efficiency A over A1995')
grid on

saveas(gcf, 'plots_transition\07_Capital_target_data.png');
close(gcf)


% sector final demand
PC_g_t_raw = PC_g_t_raw./repmat(numeraire,I,1);
PC_l_t_raw = PC_l_t_raw./repmat(numeraire,I,1);
PC_h_t_raw = PC_h_t_raw./repmat(numeraire,I,1);

PX_g_t_raw = PX_g_t_raw./repmat(numeraire,I,1);
PX_l_t_raw = PX_l_t_raw./repmat(numeraire,I,1);
PX_h_t_raw = PX_h_t_raw./repmat(numeraire,I,1);

PC_t_raw = PC_g_t_raw+PC_l_t_raw+PC_h_t_raw;
PX_t_raw = PX_g_t_raw+PX_l_t_raw+PX_h_t_raw;
F_t_raw = PC_t_raw+PX_t_raw;
F_g_t_raw = PC_g_t_raw+PX_g_t_raw;
F_l_t_raw = PC_l_t_raw+PX_l_t_raw;
F_h_t_raw = PC_h_t_raw+PX_h_t_raw;

g_t = PC_t./F_t;
g_t_raw = PC_t_raw./F_t_raw;

figure('Position', get(0, 'Screensize'))
subplot(2,3,1)
plot(1995:2011,omega_cg_t(1,:),'-b',1995:2011,omega_cg_t(2,:),'-r',1995:2011,omega_cg_t(3,:),'-g',1995:2011,omega_cg_t(4,:),'-c',1995:2011,omega_cg_t(5,:),'-k',1995:2011,omega_cg_t(7,:),'-m',1995:2011,omega_cg_t(9,:),'-y', ...
    1995:2011,omega_cg_t_raw(1,:),'--b',1995:2011,omega_cg_t_raw(2,:),'--r',1995:2011,omega_cg_t_raw(3,:),'--g',1995:2011,omega_cg_t_raw(4,:),'--c',1995:2011,omega_cg_t_raw(5,:),'--k',1995:2011,omega_cg_t_raw(7,:),'--m',1995:2011,omega_cg_t_raw(9,:),'--y','Linewidth',2.5)
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','northwest') 
title('Cons share Goods')
grid on

subplot(2,3,2)
plot(1995:2011,omega_cl_t(1,:),'-b',1995:2011,omega_cl_t(2,:),'-r',1995:2011,omega_cl_t(3,:),'-g',1995:2011,omega_cl_t(4,:),'-c',1995:2011,omega_cl_t(5,:),'-k',1995:2011,omega_cl_t(7,:),'-m',1995:2011,omega_cl_t(9,:),'-y', ...
    1995:2011,omega_cl_t_raw(1,:),'--b',1995:2011,omega_cl_t_raw(2,:),'--r',1995:2011,omega_cl_t_raw(3,:),'--g',1995:2011,omega_cl_t_raw(4,:),'--c',1995:2011,omega_cl_t_raw(5,:),'--k',1995:2011,omega_cl_t_raw(7,:),'--m',1995:2011,omega_cl_t_raw(9,:),'--y','Linewidth',2.5)
title('Cons share Low-skill services')
grid on

subplot(2,3,3)
plot(1995:2011,omega_ch_t(1,:),'-b',1995:2011,omega_ch_t(2,:),'-r',1995:2011,omega_ch_t(3,:),'-g',1995:2011,omega_ch_t(4,:),'-c',1995:2011,omega_ch_t(5,:),'-k',1995:2011,omega_ch_t(7,:),'-m',1995:2011,omega_ch_t(9,:),'-y', ...
    1995:2011,omega_ch_t_raw(1,:),'--b',1995:2011,omega_ch_t_raw(2,:),'--r',1995:2011,omega_ch_t_raw(3,:),'--g',1995:2011,omega_ch_t_raw(4,:),'--c',1995:2011,omega_ch_t_raw(5,:),'--k',1995:2011,omega_ch_t_raw(7,:),'--m',1995:2011,omega_ch_t_raw(9,:),'--y','Linewidth',2.5)
title('Cons share High-skill services')
grid on

subplot(2,3,4)
plot(1995:2011,omega_xg_t(1,:),'-b',1995:2011,omega_xg_t(2,:),'-r',1995:2011,omega_xg_t(3,:),'-g',1995:2011,omega_xg_t(4,:),'-c',1995:2011,omega_xg_t(5,:),'-k',1995:2011,omega_xg_t(7,:),'-m',1995:2011,omega_xg_t(9,:),'-y', ...
    1995:2011,omega_xg_t_raw(1,:),'--b',1995:2011,omega_xg_t_raw(2,:),'--r',1995:2011,omega_xg_t_raw(3,:),'--g',1995:2011,omega_xg_t_raw(4,:),'--c',1995:2011,omega_xg_t_raw(5,:),'--k',1995:2011,omega_xg_t_raw(7,:),'--m',1995:2011,omega_xg_t_raw(9,:),'--y','Linewidth',2.5)
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','northwest') 
title('Invest Goods')
grid on

subplot(2,3,5)
plot(1995:2011,omega_xl_t(1,:),'-b',1995:2011,omega_xl_t(2,:),'-r',1995:2011,omega_xl_t(3,:),'-g',1995:2011,omega_xl_t(4,:),'-c',1995:2011,omega_xl_t(5,:),'-k',1995:2011,omega_xl_t(7,:),'-m',1995:2011,omega_xl_t(9,:),'-y', ...
    1995:2011,omega_xl_t_raw(1,:),'--b',1995:2011,omega_xl_t_raw(2,:),'--r',1995:2011,omega_xl_t_raw(3,:),'--g',1995:2011,omega_xl_t_raw(4,:),'--c',1995:2011,omega_xl_t_raw(5,:),'--k',1995:2011,omega_xl_t_raw(7,:),'--m',1995:2011,omega_xl_t_raw(9,:),'--y','Linewidth',2.5)
title('Invest share Low-Skill Services')
grid on

subplot(2,3,6)
plot(1995:2011,omega_xh_t(1,:),'-b',1995:2011,omega_xh_t(2,:),'-r',1995:2011,omega_xh_t(3,:),'-g',1995:2011,omega_xh_t(4,:),'-c',1995:2011,omega_xh_t(5,:),'-k',1995:2011,omega_xh_t(7,:),'-m',1995:2011,omega_xh_t(9,:),'-y', ...
    1995:2011,omega_xh_t_raw(1,:),'--b',1995:2011,omega_xh_t_raw(2,:),'--r',1995:2011,omega_xh_t_raw(3,:),'--g',1995:2011,omega_xh_t_raw(4,:),'--c',1995:2011,omega_xh_t_raw(5,:),'--k',1995:2011,omega_xh_t_raw(7,:),'--m',1995:2011,omega_xh_t_raw(9,:),'--y','Linewidth',2.5)
title('Invest share High-tech Mfg')
grid on

sgtitle('Figure. Targeted moments vs raw data')
saveas(gcf, 'plots_transition\08_Target_final_demand_share.png');
close(gcf)

% back out productivity
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
diag(bt_h_t(:,:,2))
u_h_t(:,2)
P_h_t(:,2)

for i=1:num_year
    T_g_t(:,i) = (gam.^theta).*diag(bt_g_t(:,:,i)).*( u_g_t(:,i)./P_g_t(:,i) ).^theta;
    T_l_t(:,i) = (gam.^theta).*diag(bt_l_t(:,:,i)).*( u_l_t(:,i)./P_l_t(:,i) ).^theta;
    T_h_t(:,i) = (gam.^theta).*diag(bt_h_t(:,:,i)).*( u_h_t(:,i)./P_h_t(:,i) ).^theta;
end

T_g_t
T_l_t
T_h_t

% decompose TFP
figure('Position', get(0, 'Screensize'))
subplot(2,2,1)
plot(1995:2011,T_g_t(1,1:num_year)/T_g_t(1,1),'-b',1995:2011,T_g_t(2,1:num_year)/T_g_t(2,1),'-r',1995:2011,T_g_t(3,1:num_year)/T_g_t(3,1),'-g',1995:2011,T_g_t(4,1:num_year)/T_g_t(4,1),'-c',1995:2011,T_g_t(5,1:num_year)/T_g_t(5,1),'-k',1995:2011,T_g_t(7,1:num_year)/T_g_t(7,1),'-m',1995:2011,T_g_t(9,1:num_year)/T_g_t(9,1),'-y','Linewidth',2.5)
title('Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
grid on

subplot(2,2,2)
plot(1995:2011,T_l_t(1,1:num_year)/T_l_t(1,1),'-b',1995:2011,T_l_t(2,1:num_year)/T_l_t(2,1),'-r',1995:2011,T_l_t(3,1:num_year)/T_l_t(3,1),'-g',1995:2011,T_l_t(4,1:num_year)/T_l_t(4,1),'-c',1995:2011,T_l_t(5,1:num_year)/T_l_t(5,1),'-k',1995:2011,T_l_t(7,1:num_year)/T_l_t(7,1),'-m',1995:2011,T_l_t(9,1:num_year)/T_l_t(9,1),'-y','Linewidth',2.5)
title('Low-skill services')
grid on

subplot(2,2,3)
plot(1995:2011,T_g_t(1,1:num_year)/T_g_t(1,1),'-b',1995:2011,T_g_t(2,1:num_year)/T_g_t(2,1),'-r',1995:2011,T_g_t(3,1:num_year)/T_g_t(3,1),'-g',1995:2011,T_g_t(4,1:num_year)/T_g_t(4,1),'-c',1995:2011,T_g_t(5,1:num_year)/T_g_t(5,1),'-k',1995:2011,T_g_t(7,1:num_year)/T_g_t(7,1),'-m',1995:2011,T_g_t(9,1:num_year)/T_g_t(9,1),'-y','Linewidth',2.5)
title('High-skill services')
grid on

sgtitle('Figure. Productivity by sector and country relative to 1995')
saveas(gcf, 'plots_transition/final_2/01_Decompose_TFP.png');
close(gcf)

TFP_decompose(I,num_year,theta,eta,r_t,w_t,P_g_t,P_l_t,P_h_t,bt_g_t,bt_l_t,bt_h_t,...
    alpha_g_t,alpha_l_t,alpha_h_t,nu_g_t,nu_l_t,nu_h_t,...
    mu_gg_t,mu_gl_t,mu_gh_t,mu_lg_t,mu_ll_t,mu_lh_t,...
    mu_hg_t,mu_hl_t,mu_hh_t);


%Price by sector and country relative to 2000
figure('Position', get(0, 'Screensize'))
subplot(2,3,1)
plot(1995:2011,P_g_t(1,:)/P_g_t(1,1),'-b',1995:2011,P_g_t(2,:)/P_g_t(2,1),'-r',1995:2011,P_g_t(3,:)/P_g_t(3,1),'-g',1995:2011,P_g_t(4,:)/P_g_t(4,1),'-c',1995:2011,P_g_t(5,:)/P_g_t(5,1),'-k',1995:2011,P_g_t(7,:)/P_g_t(7,1),'-m',1995:2011,P_g_t(9,:)/P_g_t(9,1),'-y','Linewidth',2.5)
title('Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
grid on

subplot(2,3,2)
plot(1995:2011,P_l_t(1,:)/P_l_t(1,1),'-b',1995:2011,P_l_t(2,:)/P_l_t(2,1),'-r',1995:2011,P_l_t(3,:)/P_l_t(3,1),'-g',1995:2011,P_l_t(4,:)/P_l_t(4,1),'-c',1995:2011,P_l_t(5,:)/P_l_t(5,1),'-k',1995:2011,P_l_t(7,:)/P_l_t(7,1),'-m',1995:2011,P_l_t(9,:)/P_l_t(9,1),'-y','Linewidth',2.5)
title('Low-skill services')
grid on

subplot(2,3,3)
plot(1995:2011,P_h_t(1,:)/P_h_t(1,1),'-b',1995:2011,P_h_t(2,:)/P_h_t(2,1),'-r',1995:2011,P_h_t(3,:)/P_h_t(3,1),'-g',1995:2011,P_h_t(4,:)/P_h_t(4,1),'-c',1995:2011,P_h_t(5,:)/P_h_t(5,1),'-k',1995:2011,P_h_t(7,:)/P_h_t(7,1),'-m',1995:2011,P_h_t(9,:)/P_h_t(9,1),'-y','Linewidth',2.5)
title('High-skill services')
grid on

subplot(2,3,4)
plot(1995:2011,w_t(1,:)/w_t(1,1),'-b',1995:2011,w_t(2,:)/w_t(2,1),'-r',1995:2011,w_t(3,:)/w_t(3,1),'-g',1995:2011,w_t(4,:)/w_t(4,1),'-c',1995:2011,w_t(5,:)/w_t(5,1),'-k',1995:2011,w_t(7,:)/w_t(7,1),'-m',1995:2011,w_t(9,:)/w_t(9,1),'-y','Linewidth',2.5)
title('Wage')
grid on

subplot(2,3,5)
plot(1995:2011,r_t(1,:)/r_t(1,1),'-b',1995:2011,r_t(2,:)/r_t(2,1),'-r',1995:2011,r_t(3,:)/r_t(3,1),'-g',1995:2011,r_t(4,:)/r_t(4,1),'-c',1995:2011,r_t(5,:)/r_t(5,1),'-k',1995:2011,r_t(7,:)/r_t(7,1),'-m',1995:2011,r_t(9,:)/r_t(9,1),'-y','Linewidth',2.5)
title('Rental rate')
grid on

sgtitle('Figure. Prices by sector and country relative to 1995')
saveas(gcf, 'plots_transition/09_Prices.png');
close(gcf)

%% Given the investment efficiency, solve the steady-state


A_x2=Ax_t(:,num_year);
T_g2=T_g_t(:,num_year);
T_l2=T_l_t(:,num_year);
T_h2=T_h_t(:,num_year);

L2 = L_t(:,num_year);

d_g2 = tao_g_t(:,:,num_year);
d_l2 = tao_l_t(:,:,num_year);
d_h2 = tao_h_t(:,:,num_year);

alpha_g2=alpha_g_t(:,num_year);
alpha_l2=alpha_l_t(:,num_year);
alpha_h2=alpha_h_t(:,num_year);

nu_g2=nu_g_t(:,num_year);
nu_l2=nu_l_t(:,num_year);
nu_h2=nu_h_t(:,num_year);

mu_gg2 = mu_gg_t(:,num_year);
mu_gl2 = mu_gl_t(:,num_year);
mu_gh2 = mu_gh_t(:,num_year);

mu_lg2 = mu_lg_t(:,num_year);
mu_ll2 = mu_ll_t(:,num_year);
mu_lh2 = mu_lh_t(:,num_year);

mu_hg2 = mu_hg_t(:,num_year);
mu_hl2 = mu_hl_t(:,num_year);
mu_hh2 = mu_hh_t(:,num_year);


fi2 = fi_t(:,num_year);
omega_cg2 = omega_cg_t(:,num_year);
omega_cl2 = omega_cl_t(:,num_year);
omega_ch2 = omega_ch_t(:,num_year);

omega_xg2 = omega_xg_t(:,num_year);
omega_xl2 = omega_xl_t(:,num_year);
omega_xh2 = omega_xh_t(:,num_year);

omega_bar = zeros(I,1); % not used

diary(['plots_transition\Results.txt']);

disp(' ');
disp('---------------------------------------------------------------');
disp('Terminal steady state');

mu_sum_g = mu_gg2 + mu_gl2 + mu_gh2
mu_sum_l = mu_lg2 + mu_ll2 + mu_lh2
mu_sum_h = mu_hg2 + mu_hl2 + mu_hh2

nu_g2
nu_l2
nu_h2

alpha_g2
alpha_l2
alpha_h2

omega_c_sum = omega_cg2 + omega_cl2 + omega_ch2
omega_x_sum = omega_xg2 + omega_xl2 + omega_xh2

tic;
% Solve for equilibrium wages
w0 = w_t(:,num_year);
w2 = SS_planner(w0,A_x2,T_g2,T_l2,T_h2,d_g2,d_l2,d_h2,L2,I,...
    theta,eta,del,lambda,epsilon_c,epsilon_x,alpha_g2,alpha_l2,alpha_h2,nu_g2,nu_l2,nu_h2,...
    beta,mu_gg2,mu_gl2,mu_gh2,mu_lg2,mu_ll2,mu_lh2,mu_hg2,mu_hl2,mu_hh2,...
    fi2,omega_bar,omega_cg2,omega_cl2,omega_ch2,omega_xg2,omega_xl2,omega_xh2);

% Solve for remaining allocations
[w2_check,r2,q2,P_g2,P_l2,P_h2,P_c2,P_x2,C2,X2,K2,A2,F2,...
    K_g2,K_l2,K_h2,L_g2,L_l2,L_h2,Y_g2,Y_l2,Y_h2,pi_g2,pi_l2,pi_h2] = ...
    SS_alloc_planner(w2,A_x2,T_g2,T_l2,T_h2,d_g2,d_l2,d_h2,L2,I,...
    theta,eta,del,lambda,epsilon_c,epsilon_x,alpha_g2,alpha_l2,alpha_h2,nu_g2,nu_l2,nu_h2,...
    beta,mu_gg2,mu_gl2,mu_gh2,mu_lg2,mu_ll2,mu_lh2,mu_hg2,mu_hl2,mu_hh2,...
    fi2,omega_bar,omega_cg2,omega_cl2,omega_ch2,omega_xg2,omega_xl2,omega_xh2);

rho2 = (P_x2.*X2)./(w2.*L2+r2.*K2);
R2 = lambda.*( r2+(1-lambda)./lambda.*P_x2.*X2./K2+...
    (1-del)./lambda./A_x2.*P_x2.*(X2./K2).^(1-lambda) )./...
    ( P_x2./A_x2.*(X2./K2).^(1-lambda) )-1; % steady-state real interest rate
GDP2 = w2.*L2+r2.*K2;
PC2 = P_c2.*C2;
PX2 = P_x2.*X2;
WGDP2 = sum(GDP2);
WPC2 = sum(PC2);
WPX2 = sum(PX2);

disp(['Data at ' num2str(Last_data_year) ' vs steady-state'])
disp('wage')
[w_t(:,num_year) w2]
disp('labor')
[l_t(:,num_year) L2]
disp('rental rate')
[r_t(:,num_year) r2]
disp('capital stock')
[K_t(:,1) K_t(:,num_year) K2]

% check moments between data and model
disp('ss trade imbalance in the model')
[A2 F2]

%% Given the investment efficiency, Solve the whole transition from 2000 to the steady-state


% solve the whole transition from 2000
TT = 150*0+300*0+700;
A_x_tt = [Ax_t  repmat(Ax_t(:,num_year),1,(TT-num_year))];
T_g_tt = [T_g_t repmat(T_g_t(:,num_year),1,(TT-num_year))];
T_l_tt = [T_l_t repmat(T_l_t(:,num_year),1,(TT-num_year))];
T_h_tt = [T_h_t repmat(T_h_t(:,num_year),1,(TT-num_year))];

L_tt = [L_t repmat(L_t(:,num_year),1,(TT-num_year))];
omega_cg_tt = [omega_cg_t repmat(omega_cg_t(:,num_year),1,(TT-num_year))];
omega_cl_tt = [omega_cl_t repmat(omega_cl_t(:,num_year),1,(TT-num_year))];
omega_ch_tt = [omega_ch_t repmat(omega_ch_t(:,num_year),1,(TT-num_year))];

omega_xg_tt = [omega_xg_t repmat(omega_xg_t(:,num_year),1,(TT-num_year))];
omega_xl_tt = [omega_xl_t repmat(omega_xl_t(:,num_year),1,(TT-num_year))];
omega_xh_tt = [omega_xh_t repmat(omega_xh_t(:,num_year),1,(TT-num_year))];

fi_tt = [fi_t repmat(fi_t(:,num_year),1,(TT-num_year))];

alpha_g_tt = [alpha_g_t repmat(alpha_g_t(:,num_year),1,(TT-num_year))];
alpha_l_tt = [alpha_l_t repmat(alpha_l_t(:,num_year),1,(TT-num_year))];
alpha_h_tt = [alpha_h_t repmat(alpha_h_t(:,num_year),1,(TT-num_year))];

nu_g_tt = [nu_g_t repmat(nu_g_t(:,num_year),1,(TT-num_year))];
nu_l_tt = [nu_l_t repmat(nu_l_t(:,num_year),1,(TT-num_year))];
nu_h_tt = [nu_h_t repmat(nu_h_t(:,num_year),1,(TT-num_year))];

mu_gg_tt = [mu_gg_t repmat(mu_gg_t(:,num_year),1,(TT-num_year))];
mu_gl_tt = [mu_gl_t repmat(mu_gl_t(:,num_year),1,(TT-num_year))];
mu_gh_tt = [mu_gh_t repmat(mu_gh_t(:,num_year),1,(TT-num_year))];

mu_lg_tt = [mu_lg_t repmat(mu_lg_t(:,num_year),1,(TT-num_year))];
mu_ll_tt = [mu_ll_t repmat(mu_ll_t(:,num_year),1,(TT-num_year))];
mu_lh_tt = [mu_lh_t repmat(mu_lh_t(:,num_year),1,(TT-num_year))];

mu_hg_tt = [mu_hg_t repmat(mu_hg_t(:,num_year),1,(TT-num_year))];
mu_hl_tt = [mu_hl_t repmat(mu_hl_t(:,num_year),1,(TT-num_year))];
mu_hh_tt = [mu_hh_t repmat(mu_hh_t(:,num_year),1,(TT-num_year))];


d_g_tt = zeros(I,I,TT);
d_l_tt = zeros(I,I,TT);
d_h_tt = zeros(I,I,TT);

for tt=1:TT
    if tt<=num_year
        d_g_tt(:,:,tt) = tao_g_t(:,:,tt);
        d_l_tt(:,:,tt) = tao_l_t(:,:,tt);
        d_h_tt(:,:,tt) = tao_h_t(:,:,tt);
       
    else
        d_g_tt(:,:,tt) = tao_g_t(:,:,num_year);
        d_l_tt(:,:,tt) = tao_l_t(:,:,num_year);
        d_h_tt(:,:,tt) = tao_h_t(:,:,num_year);
        
    end
end

K_1 = K_t(:,1);
A_1 = NX_t(:,1);


disp(' ');
disp('---------------------------------------------------------------');
disp(['Computing whole transition from 1995 to ' num2str(Last_data_year) ' until final steady-state']);


% Initial guess for wage, world interest rate, and nominal investment rate.
speed = 0.75;

tic;
if Load_transition_whole == 0
    
    w0(:,1) = w_t(:,1);
    for tt=1:TT-1
        w0(:,tt+1) = w2 + speed*(w0(:,tt)-w2);
    end
    rho0 = repmat(rho2,1,TT); 
else
    
    rho0 = textread('plots_transition/whole_transition/rho_tt.txt');
    w0 = textread('plots_transition/whole_transition/w_tt.txt');
end

[rho_tt,w_tt,Zr_tt] = comp_inv_labor_planner(rho0,w0,A_x_tt,T_g_tt,T_l_tt,T_h_tt,d_g_tt,d_l_tt,d_h_tt,L_tt,K_1,A_1,I,TT,...
    theta,eta,del,lambda,epsilon_c,epsilon_x,alpha_g_tt,alpha_l_tt,alpha_h_tt,nu_g_tt,nu_l_tt,nu_h_tt,...
    beta,sigma,r2,P_c2,P_x2,C2,X2,K2,A_x2,mu_gg_tt,mu_gl_tt,mu_gh_tt,mu_lg_tt,mu_ll_tt,mu_lh_tt,...
    mu_hg_tt,mu_hl_tt,mu_hh_tt,...
    fi_tt,omega_bar,omega_cg_tt,omega_cl_tt,omega_ch_tt,omega_xg_tt,omega_xl_tt,omega_xh_tt);


% Solve for remaining allocations as a function of wages
[w_check_tt,r_tt,P_g_tt,P_l_tt,P_h_tt,P_c_tt,P_x_tt,C_tt,X_tt,K_tt,A_tt,F_tt, ...
    K_g_tt,K_l_tt,K_h_tt,L_g_tt,L_l_tt,L_h_tt,Y_g_tt,Y_l_tt,Y_h_tt,pi_g_tt,pi_l_tt,pi_h_tt] = ...
    alloc_planner(w_tt,A_x_tt,T_g_tt,T_l_tt,T_h_tt,d_g_tt,d_l_tt,d_h_tt,L_tt,K_1,A_1,I,TT,...
    theta,eta,del,lambda,epsilon_c,epsilon_x,alpha_g_tt,alpha_l_tt,alpha_h_tt,nu_g_tt,nu_l_tt,nu_h_tt,...
    rho_tt,mu_gg_tt,mu_gl_tt,mu_gh_tt,mu_lg_tt,mu_ll_tt,mu_lh_tt,mu_hg_tt,mu_hl_tt,mu_hh_tt,...
    fi_tt,omega_bar,omega_cg_tt,omega_cl_tt,omega_ch_tt,omega_xg_tt,omega_xl_tt,omega_xh_tt);

% real interest rate for return in the next period R(t+1)
R_tt = lambda.*( [r_tt(:,2:TT),r2]+(1-lambda)./lambda.*[P_x_tt(:,2:TT),P_x2].*[X_tt(:,2:TT),X2]./[K_tt(:,2:TT),K2]+...
    (1-del)./lambda./[A_x_tt(:,2:TT),A_x2].*[P_x_tt(:,2:TT),P_x2].*( [X_tt(:,2:TT),X2]./[K_tt(:,2:TT),K2] ).^(1-lambda) )./...
    ( P_x_tt(:,1:TT)./A_x_tt(:,1:TT).*( X_tt(:,1:TT)./K_tt(:,1:TT) ).^(1-lambda).*[P_c_tt(:,2:TT),P_c2]./P_c_tt(:,1:TT) )-1;
GDP_tt = w_tt.*L_tt+r_tt.*K_tt(:,1:TT);
PC_tt = P_c_tt.*C_tt;
PX_tt = P_x_tt.*X_tt;
P_k_tt = 1./lambda./A_x_tt.*(X_tt./K_tt(:,1:TT)).^(1-lambda).*P_x_tt;
WGDP_tt = sum(GDP_tt);
WPC_tt = sum(PC_tt);
WPX_tt = sum(PX_tt);

% sector output and absorption
PY_g_tt = P_g_tt.*Y_g_tt;
PY_l_tt = P_l_tt.*Y_l_tt;
PY_h_tt = P_h_tt.*Y_h_tt;

PQ_g_tt = mu_gg_tt.*(1-nu_g_tt).*PY_g_tt + mu_lg_tt.*(1-nu_l_tt).*PY_l_tt + mu_hg_tt.*(1-nu_h_tt).*PY_h_tt +...
    omega_cg_tt.*PC_tt+omega_xg_tt.*PX_tt;
PQ_l_tt = mu_gl_tt.*(1-nu_g_tt).*PY_g_tt + mu_ll_tt.*(1-nu_l_tt).*PY_l_tt + mu_hl_tt.*(1-nu_h_tt).*PY_h_tt +...
    omega_cl_tt.*PC_tt+omega_xl_tt.*PX_tt;
PQ_h_tt = mu_gh_tt.*(1-nu_g_tt).*PY_g_tt + mu_lh_tt.*(1-nu_l_tt).*PY_l_tt + mu_hh_tt.*(1-nu_h_tt).*PY_h_tt +...
    omega_ch_tt.*PC_tt+omega_xh_tt.*PX_tt;


INTM_g_tt = mu_gg_tt.*(1-nu_g_tt).*PY_g_tt + mu_lg_tt.*(1-nu_l_tt).*PY_l_tt + mu_hg_tt.*(1-nu_h_tt).*PY_h_tt ;
INTM_l_tt = mu_gl_tt.*(1-nu_g_tt).*PY_g_tt + mu_ll_tt.*(1-nu_l_tt).*PY_l_tt + mu_hl_tt.*(1-nu_h_tt).*PY_h_tt ;
INTM_h_tt = mu_gh_tt.*(1-nu_g_tt).*PY_g_tt + mu_lh_tt.*(1-nu_l_tt).*PY_l_tt + mu_hh_tt.*(1-nu_h_tt).*PY_h_tt ;

NX_g_tt = PY_g_tt-PQ_g_tt;
NX_l_tt = PY_l_tt-PQ_l_tt;
NX_h_tt = PY_h_tt-PQ_h_tt;

va_g_tt = nu_g_tt.*PY_g_tt./GDP_tt;
va_l_tt = nu_l_tt.*PY_l_tt./GDP_tt;
va_h_tt = nu_h_tt.*PY_h_tt./GDP_tt;


u_g_tt = (r_tt./(alpha_g_tt.*nu_g_tt)).^(alpha_g_tt.*nu_g_tt) .* (w_tt./((1-alpha_g_tt).*nu_g_tt)).^((1-alpha_g_tt).*nu_g_tt) ...
    .* (P_g_tt./( mu_gg_tt.*(1-nu_g_tt) )).^( mu_gg_tt.*(1-nu_g_tt) )...
    .* (P_l_tt./( mu_gl_tt.*(1-nu_g_tt) )).^( mu_gl_tt.*(1-nu_g_tt) )...
    .* (P_h_tt./( mu_gh_tt.*(1-nu_g_tt) )).^( mu_gh_tt.*(1-nu_g_tt) );

u_l_tt = (r_tt./(alpha_l_tt.*nu_l_tt)).^(alpha_l_tt.*nu_l_tt) .* (w_tt./((1-alpha_l_tt).*nu_l_tt)).^((1-alpha_l_tt).*nu_l_tt) ...
    .* (P_g_tt./( mu_lg_tt.*(1-nu_l_tt) )).^( mu_lg_tt.*(1-nu_l_tt) )...
    .* (P_l_tt./( mu_ll_tt.*(1-nu_l_tt) )).^( mu_ll_tt.*(1-nu_l_tt) )...
    .* (P_h_tt./( mu_lh_tt.*(1-nu_l_tt) )).^( mu_lh_tt.*(1-nu_l_tt) );

u_h_tt = (r_tt./(alpha_h_tt.*nu_h_tt)).^(alpha_h_tt.*nu_h_tt) .* (w_tt./((1-alpha_h_tt).*nu_h_tt)).^((1-alpha_h_tt).*nu_h_tt) ...
    .* (P_g_tt./( mu_hg_tt.*(1-nu_h_tt) )).^( mu_hg_tt.*(1-nu_h_tt) )...
    .* (P_l_tt./( mu_hl_tt.*(1-nu_h_tt) )).^( mu_hl_tt.*(1-nu_h_tt) )...
    .* (P_h_tt./( mu_hh_tt.*(1-nu_h_tt) )).^( mu_hh_tt.*(1-nu_h_tt) );

TFP_g_tt = u_g_tt./P_g_tt;
TFP_l_tt = u_l_tt./P_l_tt;
TFP_h_tt = u_h_tt./P_h_tt;

RHS = beta.*lambda.*( r_tt(:,2:num_year)+(1-lambda)./lambda.*P_x_tt(:,2:num_year).*X_tt(:,2:num_year)./K_tt(:,2:num_year)+...
    (1-del)./lambda.*P_x_tt(:,2:num_year).*X_tt(:,2:num_year)./( K_tt(:,3:num_year+1)-(1-del).*K_tt(:,2:num_year) ));
check_Euler = RHS./( P_x_tt(:,1:num_year-1).*X_tt(:,1:num_year-1)./( K_tt(:,2:num_year)-(1-del).*K_tt(:,1:num_year-1) ) )-1;

Zr_K = lambda.*beta.*( [r_tt(:,2:TT),r2]+(1-lambda)./lambda.*[P_x_tt(:,2:TT),P_x2].*[X_tt(:,2:TT),X2]./[K_tt(:,2:TT),K2]+...
    (1-del)./lambda.*[P_x_tt(:,2:TT),P_x2].*[X_tt(:,2:TT),X2]./( [K_tt(:,3:TT),K2,K2]-(1-del).*[K_tt(:,2:TT),K2] ) )./...
    ( P_x_tt(:,1:TT).*X_tt(:,1:TT)./( [K_tt(:,2:TT),K2]-(1-del).*K_tt(:,1:TT) ) )./...
    ( [P_c_tt(:,2:TT),P_c2].*[C_tt(:,2:TT),C2]./( P_c_tt(:,1:TT).*C_tt(:,1:TT) ) )-1;

max(abs(Zr_tt(:)))
mean(w_check_tt(:)./w_tt(:))

[rho2 rho_tt(:,TT)]

disp('Euler eq')
Zr_tt(:,1:num_year-1)
Zr_tt(:,1:num_year-1)-check_Euler


% investment expenditure: data vs model
disp(['Data vs model at the last data period ' num2str(Last_data_year)])
disp('GDP')
[GDP_t(:,num_year) GDP_tt(:,num_year)]
disp('PC')
[PC_t(:,num_year) PC_tt(:,num_year)]
disp('PX')
[PX_t(:,num_year) PX_tt(:,num_year)]
disp('NX')
[NX_t(:,num_year) A_tt(:,num_year)]
disp('PC+PX+NX')
[PC_t(:,num_year)+PX_t(:,num_year)+NX_t(:,num_year) PC_tt(:,num_year)+PX_tt(:,num_year)+A_tt(:,num_year)]
disp('rou')
[rho_t(:,num_year) rho_tt(:,num_year)]

disp('PY_g')
[PY_g_t(:,num_year) P_g_tt(:,num_year).*Y_g_tt(:,num_year)]
disp('PY_l')
[PY_l_t(:,num_year) P_l_tt(:,num_year).*Y_l_tt(:,num_year)]
disp('PY_h')
[PY_h_t(:,num_year) P_h_tt(:,num_year).*Y_h_tt(:,num_year)]

disp('F_g')
[F_g_t(:,num_year) PC_tt(:,num_year).*omega_cg_tt(:,num_year)+PX_tt(:,num_year).*omega_xg_tt(:,num_year)]
disp('F_l')
[F_l_t(:,num_year) PC_tt(:,num_year).*omega_cl_tt(:,num_year)+PX_tt(:,num_year).*omega_xl_tt(:,num_year)]
disp('F_h')
[F_h_t(:,num_year) PC_tt(:,num_year).*omega_ch_tt(:,num_year)+PX_tt(:,num_year).*omega_xh_tt(:,num_year)]

disp('w')
[w_t(:,num_year) w_tt(:,num_year)]
disp('r')
[r_t(:,num_year) r_tt(:,num_year)]
disp('P_g')
[P_g_t(:,num_year) P_g_tt(:,num_year)]
disp('P_l')
[P_l_t(:,num_year) P_l_tt(:,num_year)]
disp('P_h')
[P_h_t(:,num_year) P_h_tt(:,num_year)]


disp('Pai_g')
[bt_g_t(:,:,num_year) zeros(I,2) pi_g_tt(:,:,num_year)]
disp('Pai_l')
[bt_l_t(:,:,num_year) zeros(I,2) pi_l_tt(:,:,num_year)]
disp('Pai_h')
[bt_h_t(:,:,num_year) zeros(I,2) pi_h_tt(:,:,num_year)]


disp(['Calibrated Ax at ' num2str(Last_data_year)])
A_x_tt(:,num_year)

disp('w 1995')
[w_t(:,1) w_tt(:,1)]

disp('PX 1995')
[PX_t(:,1) PX_tt(:,1)]

disp(['PX' num2str(Last_data_year)])
[PX_t(:,num_year) PX_tt(:,num_year)]

disp(['K' num2str(Last_data_year)])
[K_t(:,num_year) K_tt(:,num_year)]

disp(['K' num2str(Last_data_year+1)])
[K_t(:,num_year+1) K_tt(:,num_year+1)]

disp(['K1995 ' 'K' num2str(Last_data_year)  ' Kss'])
[K_t(:,1) K_t(:,num_year) K2]

diary off;


%% check the fit of the guessed Ax, based on implied capital and investment
figure('Position', get(0, 'Screensize'))
subplot(2,3,1)
plot(1995:2011,K_tt(1,1:num_year),'-b',1995:2011,K_tt(2,1:num_year),'-r',1995:2011,K_tt(3,1:num_year),'-g',1995:2011,K_tt(4,1:num_year),'-c',1995:2011,K_tt(5,1:num_year),'-k',1995:2011,K_tt(7,1:num_year),'-m',1995:2011,K_tt(9,1:num_year),'-y', ...
    1995:2011,K_t(1,1:num_year),'--b',1995:2011,K_t(2,1:num_year),'--r',1995:2011,K_t(3,1:num_year),'--g',1995:2011,K_t(4,1:num_year),'--c',1995:2011,K_t(5,1:num_year),'--k',1995:2011,K_t(7,1:num_year),'--m',1995:2011,K_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
title('Capital stock K')
grid on

subplot(2,3,2)
plot(1995:2011,PX_tt(1,1:num_year),'-b',1995:2011,PX_tt(2,1:num_year),'-r',1995:2011,PX_tt(3,1:num_year),'-g',1995:2011,PX_tt(4,1:num_year),'-c',1995:2011,PX_tt(5,1:num_year),'-k',1995:2011,PX_tt(7,1:num_year),'-m',1995:2011,PX_tt(9,1:num_year),'-y', ...
    1995:2011,PX_t(1,1:num_year),'--b',1995:2011,PX_t(2,1:num_year),'--r',1995:2011,PX_t(3,1:num_year),'--g',1995:2011,PX_t(4,1:num_year),'--c',1995:2011,PX_t(5,1:num_year),'--k',1995:2011,PX_t(7,1:num_year),'--m',1995:2011,PX_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Investment PX')
grid on

subplot(2,3,3)
plot(1995:2011,P_x_tt(1,1:num_year),'-b',1995:2011,P_x_tt(2,1:num_year),'-r',1995:2011,P_x_tt(3,1:num_year),'-g',1995:2011,P_x_tt(4,1:num_year),'-c',1995:2011,P_x_tt(5,1:num_year),'-k',1995:2011,P_x_tt(7,1:num_year),'-m',1995:2011,P_x_tt(9,1:num_year),'-y', ...
    1995:2011,P_x_t(1,1:num_year),'--b',1995:2011,P_x_t(2,1:num_year),'--r',1995:2011,P_x_t(3,1:num_year),'--g',1995:2011,P_x_t(4,1:num_year),'--c',1995:2011,P_x_t(5,1:num_year),'-k',1995:2011,P_x_t(7,1:num_year),'--m',1995:2011,P_x_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Investment price P_x')
grid on

subplot(2,3,4)
plot(1995:2011,X_tt(1,1:num_year),'-b',1995:2011,X_tt(2,1:num_year),'-r',1995:2011,X_tt(3,1:num_year),'-g',1995:2011,X_tt(4,1:num_year),'-c',1995:2011,X_tt(5,1:num_year),'-k',1995:2011,X_tt(7,1:num_year),'-m',1995:2011,X_tt(9,1:num_year),'-y', ...
    1995:2011,X_t(1,1:num_year),'--b',1995:2011,X_t(2,1:num_year),'--r',1995:2011,X_t(3,1:num_year),'--g',1995:2011,X_t(4,1:num_year),'--c',1995:2011,X_t(5,1:num_year),'--k',1995:2011,X_t(7,1:num_year),'--m',1995:2011,X_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Investment quantity X')
grid on

subplot(2,3,5)
plot(1995:2011,X_tt(1,1:num_year)./K_tt(1,1:num_year),'-b',1995:2011,X_tt(2,1:num_year)./K_tt(2,1:num_year),'-r',1995:2011,X_tt(3,1:num_year)./K_tt(3,1:num_year),'-g',1995:2011,X_tt(4,1:num_year)./K_tt(4,1:num_year),'-c',1995:2011,X_tt(5,1:num_year)./K_tt(5,1:num_year),'-k',1995:2011,X_tt(7,1:num_year)./K_tt(7,1:num_year),'-m',1995:2011,X_tt(9,1:num_year)./K_tt(9,1:num_year),'-y', ...
    1995:2011,X_t(1,1:num_year)./K_t(1,1:num_year),'--b',1995:2011,X_t(2,1:num_year)./K_t(2,1:num_year),'--r',1995:2011,X_t(3,1:num_year)./K_t(3,1:num_year),'--g',1995:2011,X_t(4,1:num_year)./K_t(4,1:num_year),'--c',1995:2011,X_t(5,1:num_year)./K_t(5,1:num_year),'--k',1995:2011,X_t(7,1:num_year)./K_t(7,1:num_year),'--m',1995:2011,X_t(9,1:num_year)./K_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Investment over capital X/K')
grid on

subplot(2,3,6)
plot(1995:2011,K_tt(1,1:num_year)./K_tt(1,1),'-b',1995:2011,K_tt(2,1:num_year)./K_tt(2,1),'-r',1995:2011,K_tt(3,1:num_year)./K_tt(3,1),'-g',1995:2011,K_tt(4,1:num_year)./K_tt(4,1),'-c',1995:2011,K_tt(5,1:num_year)./K_tt(5,1),'-k',1995:2011,K_tt(7,1:num_year)./K_tt(7,1),'-m',1995:2011,K_tt(9,1:num_year)./K_tt(9,1),'-y', ...
    1995:2011,K_t(1,1:num_year)./K_t(1,1),'--b',1995:2011,K_t(2,1:num_year)./K_t(2,1),'--r',1995:2011,K_t(3,1:num_year)./K_t(3,1),'--g',1995:2011,K_t(4,1:num_year)./K_t(4,1),'--c',1995:2011,K_t(5,1:num_year)./K_t(5,1),'--k',1995:2011,K_t(7,1:num_year)./K_t(7,1),'--m',1995:2011,K_t(9,1:num_year)./K_t(9,1),'--y', ...
    'Linewidth',2.5)
title('Capital growth Kt/K1995')
grid on

saveas(gcf, 'plots_transition\10_Check_fit_Ax.png');
close(gcf)



%% Plot the whole transition to the steady-state
Tss = 50;

figure('Position', get(0, 'Screensize'))
subplot(2,3,1)
plot(1:TT+Tss,[w_tt(1,:),w2(1)*ones(1,Tss)],'-b',1:TT+Tss,[w_tt(2,:),w2(2)*ones(1,Tss)],'-r',1:TT+Tss,[w_tt(3,:),w2(3)*ones(1,Tss)],'-g',1:TT+Tss,[w_tt(4,:),w2(4)*ones(1,Tss)],'-c',1:TT+Tss,[w_tt(5,:),w2(5)*ones(1,Tss)],'-k',1:TT+Tss,[w_tt(7,:),w2(7)*ones(1,Tss)],'-m',1:TT+Tss,[w_tt(9,:),w2(9)*ones(1,Tss)],'-y', ...
    'Linewidth',2.5)
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
title('Wage')
grid on

subplot(2,3,2)
plot(1:TT+Tss,[r_tt(1,:),r2(1)*ones(1,Tss)],'-b',1:TT+Tss,[r_tt(2,:),r2(2)*ones(1,Tss)],'-r',1:TT+Tss,[r_tt(3,:),r2(3)*ones(1,Tss)],'-g',1:TT+Tss,[r_tt(4,:),r2(4)*ones(1,Tss)],'-c',1:TT+Tss,[r_tt(5,:),r2(5)*ones(1,Tss)],'-k',1:TT+Tss,[r_tt(7,:),r2(7)*ones(1,Tss)],'-m',1:TT+Tss,[r_tt(9,:),r2(9)*ones(1,Tss)],'-y', ...
    'Linewidth',2.5)
title('Rental rate')
%legend('country 1','country 2','country 3','country 4','location','best')
grid on

subplot(2,3,3)
plot(1:TT+Tss,[K_tt(1,:),K2(1)*ones(1,Tss-1)],'-b',1:TT+Tss,[K_tt(2,:),K2(2)*ones(1,Tss-1)],'-r',1:TT+Tss,[K_tt(3,:),K2(3)*ones(1,Tss-1)],'-g',1:TT+Tss,[K_tt(4,:),K2(4)*ones(1,Tss-1)],'-c',1:TT+Tss,[K_tt(5,:),K2(5)*ones(1,Tss-1)],'-k',1:TT+Tss,[K_tt(7,:),K2(7)*ones(1,Tss-1)],'-m',1:TT+Tss,[K_tt(9,:),K2(9)*ones(1,Tss-1)],'-y', ...
    'Linewidth',2.5)
title('Capital accumulation')
%legend('country 1','country 2','country 3','country 4','location','best')
grid on

subplot(2,3,4)
plot(1:TT+Tss,[K_tt(1,:)./[L_tt(1,:) L2(1)],K2(1)/L2(1)*ones(1,Tss-1)],'-b',1:TT+Tss,[K_tt(2,:)./[L_tt(2,:) L2(2)],K2(2)/L2(2)*ones(1,Tss-1)],'-r',1:TT+Tss,[K_tt(3,:)./[L_tt(3,:) L2(3)],K2(3)/L2(3)*ones(1,Tss-1)],'-g',1:TT+Tss,[K_tt(4,:)./[L_tt(4,:) L2(4)],K2(4)/L2(4)*ones(1,Tss-1)],'-c',1:TT+Tss,[K_tt(5,:)./[L_tt(5,:) L2(5)],K2(5)/L2(5)*ones(1,Tss-1)],'-k',1:TT+Tss,[K_tt(7,:)./[L_tt(7,:) L2(7)],K2(7)/L2(7)*ones(1,Tss-1)],'-m',1:TT+Tss,[K_tt(9,:)./[L_tt(9,:) L2(9)],K2(9)/L2(9)*ones(1,Tss-1)],'-y', ...
    'Linewidth',2.5)
title('Capital per worker')
%legend('country 1','country 2','country 3','country 4','location','best')
grid on

subplot(2,3,5)
plot(1:TT+Tss,[C_tt(1,:),C2(1)*ones(1,Tss)],'-b',1:TT+Tss,[C_tt(2,:),C2(2)*ones(1,Tss)],'-r',1:TT+Tss,[C_tt(3,:),C2(3)*ones(1,Tss)],'-g',1:TT+Tss,[C_tt(4,:),C2(4)*ones(1,Tss)],'-c',1:TT+Tss,[C_tt(5,:),C2(5)*ones(1,Tss)],'-k',1:TT+Tss,[C_tt(7,:),C2(7)*ones(1,Tss)],'-m',1:TT+Tss,[C_tt(9,:),C2(9)*ones(1,Tss)],'-y', ...
    'Linewidth',2.5)
title('Consumption')
%legend('country 1','country 2','country 3','country 4','location','best')
grid on

subplot(2,3,6)
plot(1:TT+Tss,[X_tt(1,:),X2(1)*ones(1,Tss)],'-b',1:TT+Tss,[X_tt(2,:),X2(2)*ones(1,Tss)],'-r',1:TT+Tss,[X_tt(3,:),X2(3)*ones(1,Tss)],'-g',1:TT+Tss,[X_tt(4,:),X2(4)*ones(1,Tss)],'-c',1:TT+Tss,[X_tt(5,:),X2(5)*ones(1,Tss)],'-k',1:TT+Tss,[X_tt(7,:),X2(7)*ones(1,Tss)],'-m',1:TT+Tss,[X_tt(9,:),X2(9)*ones(1,Tss)],'-y', ...
    'Linewidth',2.5)
title('Investment')
%legend('country 1','country 2','country 3','country 4','location','best')
grid on

saveas(gcf, 'plots_transition/11_Aggregate.png');


figure('Position', get(0, 'Screensize'))
subplot(3,3,1)
plot(1:TT+Tss,[K_tt(1,:),K2(1)*ones(1,Tss-1)],'-b','Linewidth',2.5)
title('AUS')
grid on

subplot(3,3,2)
plot(1:TT+Tss,[K_tt(2,:),K2(2)*ones(1,Tss-1)],'-r','Linewidth',2.5)
title('BRA')
grid on

subplot(3,3,3)
plot(1:TT+Tss,[K_tt(3,:),K2(3)*ones(1,Tss-1)],'-g','Linewidth',2.5)
title('CHN')
grid on

subplot(3,3,4)
plot(1:TT+Tss,[K_tt(4,:),K2(4)*ones(1,Tss-1)],'-c','Linewidth',2.5)
title('DEU')
grid on

subplot(3,3,5)
plot(1:TT+Tss,[K_tt(5,:),K2(5)*ones(1,Tss-1)],'-k','Linewidth',2.5)
title('IND')
grid on

subplot(3,3,6)
plot(1:TT+Tss,[K_tt(6,:),K2(6)*ones(1,Tss-1)],'-r','Linewidth',2.5)
title('JPN')
grid on

subplot(3,3,7)
plot(1:TT+Tss,[K_tt(7,:),K2(7)*ones(1,Tss-1)],'-m','Linewidth',2.5)
title('MEX')
grid on

subplot(3,3,8)
plot(1:TT+Tss,[K_tt(9,:),K2(9)*ones(1,Tss-1)],'-y','Linewidth',2.5)
title('USA')
grid on

subplot(3,3,9)
plot(1:TT+Tss,[K_tt(8,:),K2(8)*ones(1,Tss-1)],'-c','Linewidth',2.5)
title('ROW')
grid on

saveas(gcf, 'plots_transition/12_Capital_dynamics_by_country.png');


figure('Position', get(0, 'Screensize'))
subplot(2,4,1)
plot(1:TT+Tss,[rho_tt(1,:),rho2(1)*ones(1,Tss)],'-b',1:TT+Tss,[rho_tt(2,:),rho2(2)*ones(1,Tss)],'-r',1:TT+Tss,[rho_tt(3,:),rho2(3)*ones(1,Tss)],'-g',1:TT+Tss,[rho_tt(4,:),rho2(4)*ones(1,Tss)],'-c',1:TT+Tss,[rho_tt(5,:),rho2(5)*ones(1,Tss)],'-k',1:TT+Tss,[rho_tt(7,:),rho2(7)*ones(1,Tss)],'-m',1:TT+Tss,[rho_tt(9,:),rho2(9)*ones(1,Tss)],'-y', ...
    'Linewidth',2.5)
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
title('Investment rate')

subplot(2,4,2)
plot(1:TT,Zr_tt(1,:),'-b',1:TT,Zr_tt(2,:),'-r',1:TT,Zr_tt(3,:),'-g',1:TT,Zr_tt(4,:),'-c',1:TT,Zr_tt(5,:),'-k',1:TT,Zr_tt(7,:),'-m',1:TT,Zr_tt(9,:),'-y', ...
    'Linewidth',2.5)
title('Diff in Euler eq')
%legend('country 1','country 2','country 3','country 4','location','best')

subplot(2,4,3)
plot(1:TT+Tss,[WGDP_tt,WGDP2*ones(1,Tss)],'-b',1:TT+Tss,[WPC_tt,WPC2*ones(1,Tss)],'-r',1:TT+Tss,[WPX_tt,WPX2*ones(1,Tss)],'-g','Linewidth',2.5)
title('Numeraire')
legend('World GDP','World PC','World PX','location','best')

subplot(2,4,4)
plot(1:TT+Tss,[-A_tt,-A2*ones(1,Tss)],'-b',1:TT+Tss,[F_tt,F2*ones(1,Tss)],'--r','Linewidth',2.5)
title('Net transfer inflow vs Trade deficit')

subplot(2,4,5)
plot(1:TT+Tss,[P_c_tt(1,:).*C_tt(1,:),P_c2(1)*C2(1)*ones(1,Tss)],'-b',1:TT+Tss,[P_c_tt(2,:).*C_tt(2,:),P_c2(2)*C2(2)*ones(1,Tss)],'-r',1:TT+Tss,[P_c_tt(3,:).*C_tt(3,:),P_c2(3)*C2(3)*ones(1,Tss)],'-g',1:TT+Tss,[P_c_tt(4,:).*C_tt(4,:),P_c2(4)*C2(4)*ones(1,Tss)],'-c',1:TT+Tss,[P_c_tt(5,:).*C_tt(5,:),P_c2(5)*C2(5)*ones(1,Tss)],'-k',1:TT+Tss,[P_c_tt(7,:).*C_tt(7,:),P_c2(7)*C2(7)*ones(1,Tss)],'-m',1:TT+Tss,[P_c_tt(9,:).*C_tt(9,:),P_c2(9)*C2(9)*ones(1,Tss)],'-y', ...
    'Linewidth',2.5)
title('Consumption expenditure (PC)')
%legend('country 1','country 2','country 3','country 4','location','best')

subplot(2,4,6)
plot(1:TT+Tss,[P_x_tt(1,:).*X_tt(1,:),P_x2(1)*X2(1)*ones(1,Tss)],'-b',1:TT+Tss,[P_x_tt(2,:).*X_tt(2,:),P_x2(2)*X2(2)*ones(1,Tss)],'-r',1:TT+Tss,[P_x_tt(3,:).*X_tt(3,:),P_x2(3)*X2(3)*ones(1,Tss)],'-g',1:TT+Tss,[P_x_tt(4,:).*X_tt(4,:),P_x2(4)*X2(4)*ones(1,Tss)],'-c',1:TT+Tss,[P_x_tt(5,:).*X_tt(5,:),P_x2(5)*X2(5)*ones(1,Tss)],'-k',1:TT+Tss,[P_x_tt(7,:).*X_tt(7,:),P_x2(7)*X2(7)*ones(1,Tss)],'-m',1:TT+Tss,[P_x_tt(9,:).*X_tt(9,:),P_x2(9)*X2(9)*ones(1,Tss)],'-y', ...
    'Linewidth',2.5)
title('Investment expenditure (PX)')
%legend('country 1','country 2','country 3','country 4','location','best')

subplot(2,4,7)
plot(1:TT,rho_tt,'-b',1:TT,repmat(rho2,1,TT),'--r','Linewidth',2.5)
title('rho: solution vs initial guess')

subplot(2,4,8)
plot(1:TT,Zr_K(1,:),'-b',1:TT,Zr_K(2,:),'-r',1:TT,Zr_K(3,:),'-g',1:TT,Zr_K(4,:),'-c',1:TT,Zr_K(5,:),'-k',1:TT,Zr_K(7,:),'-m',1:TT,Zr_K(9,:),'-y', ...
    'Linewidth',2.5)
title('Diff in Euler eq K')
%legend('country 1','country 2','country 3','country 4','location','best')

saveas(gcf, 'plots_transition/13_rho.png');

 
% Sector capital and labor share
figure('Position', get(0, 'Screensize'))
subplot(2,4,1)
plot(1:TT+Tss,[K_g_tt(1,:)./K_tt(1,1:TT),K_g2(1)/K2(1)*ones(1,Tss)],'-b',1:TT+Tss,[K_l_tt(1,:)./K_tt(1,1:TT),K_l2(1)/K2(1)*ones(1,Tss)],'-r',1:TT+Tss,[K_h_tt(1,:)./K_tt(1,1:TT),K_h2(1)/K2(1)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Capital share AUS')
legend('G','L','H','location','best')

subplot(2,4,2)
plot(1:TT+Tss,[K_g_tt(2,:)./K_tt(2,1:TT),K_g2(2)/K2(2)*ones(1,Tss)],'-b',1:TT+Tss,[K_l_tt(2,:)./K_tt(2,1:TT),K_l2(2)/K2(2)*ones(1,Tss)],'-r',1:TT+Tss,[K_h_tt(2,:)./K_tt(2,1:TT),K_h2(2)/K2(2)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Capital share BRA')
legend('G','L','H','location','best')

subplot(2,4,3)
plot(1:TT+Tss,[K_g_tt(3,:)./K_tt(3,1:TT),K_g2(3)/K2(3)*ones(1,Tss)],'-b',1:TT+Tss,[K_l_tt(3,:)./K_tt(3,1:TT),K_l2(3)/K2(3)*ones(1,Tss)],'-r',1:TT+Tss,[K_h_tt(3,:)./K_tt(3,1:TT),K_h2(3)/K2(3)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Capital share CHN')
legend('G','L','H','location','best')

subplot(2,4,4)
plot(1:TT+Tss,[K_g_tt(4,:)./K_tt(4,1:TT),K_g2(4)/K2(4)*ones(1,Tss)],'-b',1:TT+Tss,[K_l_tt(4,:)./K_tt(4,1:TT),K_l2(4)/K2(4)*ones(1,Tss)],'-r',1:TT+Tss,[K_h_tt(4,:)./K_tt(4,1:TT),K_h2(4)/K2(4)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Capital share DEU')
legend('G','L','H','location','best')

subplot(2,4,5)
plot(1:TT+Tss,[K_g_tt(5,:)./K_tt(5,1:TT),K_g2(5)/K2(5)*ones(1,Tss)],'-b',1:TT+Tss,[K_l_tt(5,:)./K_tt(5,1:TT),K_l2(5)/K2(5)*ones(1,Tss)],'-r',1:TT+Tss,[K_h_tt(5,:)./K_tt(5,1:TT),K_h2(5)/K2(5)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Capital share IND')
legend('G','L','H','location','best')

subplot(2,4,6)
plot(1:TT+Tss,[K_g_tt(7,:)./K_tt(7,1:TT),K_g2(7)/K2(7)*ones(1,Tss)],'-b',1:TT+Tss,[K_l_tt(7,:)./K_tt(7,1:TT),K_l2(7)/K2(7)*ones(1,Tss)],'-r',1:TT+Tss,[K_h_tt(7,:)./K_tt(7,1:TT),K_h2(7)/K2(7)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Capital share MEX')
legend('G','L','H','location','best')

subplot(2,4,7)
plot(1:TT+Tss,[K_g_tt(9,:)./K_tt(9,1:TT),K_g2(9)/K2(9)*ones(1,Tss)],'-b',1:TT+Tss,[K_l_tt(9,:)./K_tt(9,1:TT),K_l2(9)/K2(9)*ones(1,Tss)],'-r',1:TT+Tss,[K_h_tt(9,:)./K_tt(9,1:TT),K_h2(9)/K2(9)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Capital share USA')
legend('G','L','H','location','best')

subplot(2,4,8)
plot(1:TT+Tss,[K_g_tt(8,:)./K_tt(8,1:TT),K_g2(8)/K2(8)*ones(1,Tss)],'-b',1:TT+Tss,[K_l_tt(8,:)./K_tt(8,1:TT),K_l2(8)/K2(8)*ones(1,Tss)],'-r',1:TT+Tss,[K_h_tt(8,:)./K_tt(8,1:TT),K_h2(8)/K2(8)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Capital share ROW')
legend('G','L','H','location','best')

saveas(gcf, 'plots_transition/14_Sector_Capital_share.png');

subplot(2,4,1)
plot(1:TT+Tss,[L_g_tt(1,:)./L_tt(1,:),L_g2(1)/L2(1)*ones(1,Tss)],'-b',1:TT+Tss,[L_l_tt(1,:)./L_tt(1,:),L_l2(1)/L2(1)*ones(1,Tss)],'-r',1:TT+Tss,[L_h_tt(1,:)./L_tt(1,:),L_h2(1)/L2(1)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Labor share AUS')
legend('G','L','H','location','best')

subplot(2,4,2)
plot(1:TT+Tss,[L_g_tt(2,:)./L_tt(2,:),L_g2(2)/L2(2)*ones(1,Tss)],'-b',1:TT+Tss,[L_l_tt(2,:)./L_tt(2,:),L_l2(2)/L2(2)*ones(1,Tss)],'-r',1:TT+Tss,[L_h_tt(2,:)./L_tt(2,:),L_h2(2)/L2(2)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Labor share BRA')
legend('G','L','H','location','best')

subplot(2,4,3)
plot(1:TT+Tss,[L_g_tt(3,:)./L_tt(3,:),L_g2(3)/L2(3)*ones(1,Tss)],'-b',1:TT+Tss,[L_l_tt(3,:)./L_tt(3,:),L_l2(3)/L2(3)*ones(1,Tss)],'-r',1:TT+Tss,[L_h_tt(3,:)./L_tt(3,:),L_h2(3)/L2(3)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Labor share CHN')
legend('G','L','H','location','best')

subplot(2,4,4)
plot(1:TT+Tss,[L_g_tt(4,:)./L_tt(4,:),L_g2(4)/L2(4)*ones(1,Tss)],'-b',1:TT+Tss,[L_l_tt(4,:)./L_tt(4,:),L_l2(4)/L2(4)*ones(1,Tss)],'-r',1:TT+Tss,[L_h_tt(4,:)./L_tt(4,:),L_h2(4)/L2(4)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Labor share DEU')
legend('G','L','H','location','best')

subplot(2,4,5)
plot(1:TT+Tss,[L_g_tt(5,:)./L_tt(5,:),L_g2(5)/L2(5)*ones(1,Tss)],'-b',1:TT+Tss,[L_l_tt(5,:)./L_tt(5,:),L_l2(5)/L2(5)*ones(1,Tss)],'-r',1:TT+Tss,[L_h_tt(5,:)./L_tt(5,:),L_h2(5)/L2(5)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Labor share IND')
legend('G','L','H','location','best')

subplot(2,4,6)
plot(1:TT+Tss,[L_g_tt(7,:)./L_tt(7,:),L_g2(7)/L2(7)*ones(1,Tss)],'-b',1:TT+Tss,[L_l_tt(7,:)./L_tt(7,:),L_l2(7)/L2(7)*ones(1,Tss)],'-r',1:TT+Tss,[L_h_tt(7,:)./L_tt(7,:),L_h2(7)/L2(7)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Labor share MEX')
legend('G','L','H','location','best')

subplot(2,4,7)
plot(1:TT+Tss,[L_g_tt(9,:)./L_tt(9,:),L_g2(9)/L2(9)*ones(1,Tss)],'-b',1:TT+Tss,[L_l_tt(9,:)./L_tt(9,:),L_l2(9)/L2(9)*ones(1,Tss)],'-r',1:TT+Tss,[L_h_tt(9,:)./L_tt(9,:),L_h2(9)/L2(9)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Labor share USA')
legend('G','L','H','location','best')

subplot(2,4,8)
plot(1:TT+Tss,[L_g_tt(8,:)./L_tt(8,:),L_g2(8)/L2(8)*ones(1,Tss)],'-b',1:TT+Tss,[L_l_tt(8,:)./L_tt(8,:),L_l2(8)/L2(8)*ones(1,Tss)],'-r',1:TT+Tss,[L_h_tt(8,:)./L_tt(8,:),L_h2(8)/L2(8)*ones(1,Tss)],'-g','Linewidth',2.5)
title('Labor share ROW')
legend('G','L','H','location','best')

saveas(gcf, 'plots_transition/15_Sector_Labor_share.png');


%% Compare solutions and targets 
mkdir plots_transition/whole_transition

if Save_transition == 1 | Load_transition_whole == 0
save('plots_transition/whole_transition/rho_tt.txt','rho_tt','-ASCII');
save('plots_transition/whole_transition/w_tt.txt','w_tt','-ASCII');
save('plots_transition/whole_transition/Zr_tt.txt','Zr_tt','-ASCII');
save('plots_transition/whole_transition/del.txt','del','-ASCII');
save('plots_transition/whole_transition/A_x_guess.txt','A_x_guess','-ASCII');
end 

ye = 2030-1995+1; % end year for plot


%Checking A_x_guess
figure('Position', get(0, 'Screensize'))
subplot(3,3,1)
plot(1995:1994+ye,PX_tt(1,1:ye),'-b',...
    1995:1995+num_year-1,PX_t(1,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment PX AUS')

subplot(3,3,2)
plot(1995:1994+ye,PX_tt(2,1:ye),'-b',...
    1995:1995+num_year-1,PX_t(2,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment PX BRA')

subplot(3,3,3)
plot(1995:1994+ye,PX_tt(3,1:ye),'-b',...
    1995:1995+num_year-1,PX_t(3,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment PX CHN')

subplot(3,3,4)
plot(1995:1994+ye,PX_tt(4,1:ye),'-b',...
    1995:1995+num_year-1,PX_t(4,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment PX DEU')

subplot(3,3,5)
plot(1995:1994+ye,PX_tt(5,1:ye),'-b',...
    1995:1995+num_year-1,PX_t(5,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment PX IND')

subplot(3,3,6)
plot(1995:1994+ye,PX_tt(6,1:ye),'-b',...
    1995:1995+num_year-1,PX_t(6,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment PX JPN')

subplot(3,3,7)
plot(1995:1994+ye,PX_tt(7,1:ye),'-b',...
    1995:1995+num_year-1,PX_t(7,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment PX MEX')

subplot(3,3,8)
plot(1995:1994+ye,PX_tt(8,1:ye),'-b',...
    1995:1995+num_year-1,PX_t(8,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment PX row')

subplot(3,3,9)
plot(1995:1994+ye,PX_tt(9,1:ye),'-b',...
    1995:1995+num_year-1,PX_t(9,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment PX USA')

sgtitle('Figure. Model fit for aggregate allocation')
saveas(gcf, 'plots_transition/1601_Fit_agg_investmentPX.png');
close(gcf)


figure('Position', get(0, 'Screensize'))
subplot(3,3,1)
plot(1995:1994+ye,rho_tt(1,1:ye),'-b',...
    1995:1995+num_year-1,rho_t(1,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment rate AUS')

subplot(3,3,2)
plot(1995:1994+ye,rho_tt(2,1:ye),'-b',...
    1995:1995+num_year-1,rho_t(2,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment rate BRA')

subplot(3,3,3)
plot(1995:1994+ye,rho_tt(3,1:ye),'-b',...
    1995:1995+num_year-1,rho_t(3,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment rate CHN')

subplot(3,3,4)
plot(1995:1994+ye,rho_tt(4,1:ye),'-b',...
    1995:1995+num_year-1,rho_t(4,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment rate DEU')

subplot(3,3,5)
plot(1995:1994+ye,rho_tt(5,1:ye),'-b',...
    1995:1995+num_year-1,rho_t(5,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment rate IND')

subplot(3,3,6)
plot(1995:1994+ye,rho_tt(6,1:ye),'-b',...
    1995:1995+num_year-1,rho_t(6,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment rate JPN')

subplot(3,3,7)
plot(1995:1994+ye,rho_tt(7,1:ye),'-b',...
    1995:1995+num_year-1,rho_t(7,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment rate MEX')

subplot(3,3,8)
plot(1995:1994+ye,rho_tt(8,1:ye),'-b',...
    1995:1995+num_year-1,rho_t(8,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment rate row')

subplot(3,3,9)
plot(1995:1994+ye,rho_tt(9,1:ye),'-b',...
    1995:1995+num_year-1,rho_t(9,1:num_year),'--b', ...
    'Linewidth',2.5)
title('Investment rate USA')

sgtitle('Figure. Model fit for aggregate allocation')
saveas(gcf, 'plots_transition/1602_Fit_agg_rho.png');
close(gcf)


% Aggregate allocation
figure('Position', get(0, 'Screensize'))
subplot(2,4,1)
plot(1995:1994+ye,GDP_tt(1,1:ye),'-b',1995:1994+ye,GDP_tt(2,1:ye),'-r',1995:1994+ye,GDP_tt(3,1:ye),'-g',1995:1994+ye,GDP_tt(4,1:ye),'-c',1995:1994+ye,GDP_tt(5,1:ye),'-k',1995:1994+ye,GDP_tt(7,1:ye),'-m',1995:1994+ye,GDP_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,GDP_t(1,1:num_year),'--b',1995:1995+num_year-1,GDP_t(2,1:num_year),'--r',1995:1995+num_year-1,GDP_t(3,1:num_year),'--g',1995:1995+num_year-1,GDP_t(4,1:num_year),'--c',1995:1995+num_year-1,GDP_t(5,1:num_year),'--k',1995:1995+num_year-1,GDP_t(7,1:num_year),'--m',1995:1995+num_year-1,GDP_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
title('GDP = PC+PX+NX')

subplot(2,4,2)
plot(1995:1994+ye,PC_tt(1,1:ye),'-b',1995:1994+ye,PC_tt(2,1:ye),'-r',1995:1994+ye,PC_tt(3,1:ye),'-g',1995:1994+ye,PC_tt(4,1:ye),'-c',1995:1994+ye,PC_tt(5,1:ye),'-k',1995:1994+ye,PC_tt(7,1:ye),'-m',1995:1994+ye,PC_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,PC_t(1,1:num_year),'--b',1995:1995+num_year-1,PC_t(2,1:num_year),'--r',1995:1995+num_year-1,PC_t(3,1:num_year),'--g',1995:1995+num_year-1,PC_t(4,1:num_year),'--c',1995:1995+num_year-1,PC_t(5,1:num_year),'--k',1995:1995+num_year-1,PC_t(7,1:num_year),'--m',1995:1995+num_year-1,PC_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Consumption PC')

subplot(2,4,3)
plot(1995:1994+ye,PX_tt(1,1:ye),'-b',1995:1994+ye,PX_tt(2,1:ye),'-r',1995:1994+ye,PX_tt(3,1:ye),'-g',1995:1994+ye,PX_tt(4,1:ye),'-c',1995:1994+ye,PX_tt(5,1:ye),'-k',1995:1994+ye,PX_tt(7,1:ye),'-m',1995:1994+ye,PX_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,PX_t(1,1:num_year),'--b',1995:1995+num_year-1,PX_t(2,1:num_year),'--r',1995:1995+num_year-1,PX_t(3,1:num_year),'--g',1995:1995+num_year-1,PX_t(4,1:num_year),'--c',1995:1995+num_year-1,PX_t(5,1:num_year),'--k',1995:1995+num_year-1,PX_t(7,1:num_year),'--m',1995:1995+num_year-1,PX_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Investment PX')

subplot(2,4,4)
plot(1995:1994+ye,-F_tt(1,1:ye),'-b',1995:1994+ye,-F_tt(2,1:ye),'-r',1995:1994+ye,-F_tt(3,1:ye),'-g',1995:1994+ye,-F_tt(4,1:ye),'-c',1995:1994+ye,-F_tt(5,1:ye),'-k',1995:1994+ye,-F_tt(7,1:ye),'-m',1995:1994+ye,-F_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,NX_t(1,1:num_year),'--b',1995:1995+num_year-1,NX_t(2,1:num_year),'--r',1995:1995+num_year-1,NX_t(3,1:num_year),'--g',1995:1995+num_year-1,NX_t(4,1:num_year),'--c',1995:1995+num_year-1,NX_t(5,1:num_year),'--k',1995:1995+num_year-1,NX_t(7,1:num_year),'--m',1995:1995+num_year-1,NX_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Net export NX')

subplot(2,4,5)
plot(1995:1994+ye,K_tt(1,1:ye),'-b',1995:1994+ye,K_tt(2,1:ye),'-r',1995:1994+ye,K_tt(3,1:ye),'-g',1995:1994+ye,K_tt(4,1:ye),'-c',1995:1994+ye,K_tt(5,1:ye),'-k',1995:1994+ye,K_tt(7,1:ye),'-m',1995:1994+ye,K_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,K_t(1,1:num_year),'--b',1995:1995+num_year-1,K_t(2,1:num_year),'--r',1995:1995+num_year-1,K_t(3,1:num_year),'--g',1995:1995+num_year-1,K_t(4,1:num_year),'--c',1995:1995+num_year-1,K_t(5,1:num_year),'--k',1995:1995+num_year-1,K_t(7,1:num_year),'--m',1995:1995+num_year-1,K_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Capital stock')

subplot(2,4,6)
plot(1995:1994+ye,L_tt(1,1:ye),'-b',1995:1994+ye,L_tt(2,1:ye),'-r',1995:1994+ye,L_tt(3,1:ye),'-g',1995:1994+ye,L_tt(4,1:ye),'-c',1995:1994+ye,L_tt(5,1:ye),'-k',1995:1994+ye,L_tt(7,1:ye),'-m',1995:1994+ye,L_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,l_t(1,1:num_year),'--b',1995:1995+num_year-1,l_t(2,1:num_year),'--r',1995:1995+num_year-1,l_t(3,1:num_year),'--g',1995:1995+num_year-1,l_t(4,1:num_year),'--c',1995:1995+num_year-1,l_t(5,1:num_year),'--k',1995:1995+num_year-1,l_t(7,1:num_year),'--m',1995:1995+num_year-1,l_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Employment')

subplot(2,4,7)
plot(1995:1994+ye,rho_tt(1,1:ye),'-b',1995:1994+ye,rho_tt(2,1:ye),'-r',1995:1994+ye,rho_tt(3,1:ye),'-g',1995:1994+ye,rho_tt(4,1:ye),'-c',1995:1994+ye,rho_tt(5,1:ye),'-k',1995:1994+ye,rho_tt(7,1:ye),'-m',1995:1994+ye,rho_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,rho_t(1,1:num_year),'--b',1995:1995+num_year-1,rho_t(2,1:num_year),'--r',1995:1995+num_year-1,rho_t(3,1:num_year),'--g',1995:1995+num_year-1,rho_t(4,1:num_year),'--c',1995:1995+num_year-1,rho_t(5,1:num_year),'--k',1995:1995+num_year-1,rho_t(7,1:num_year),'--m',1995:1995+num_year-1,rho_t(9,1:num_year),'--y',...
    'Linewidth',2.5)
title('Investment rate')

subplot(2,4,8)
plot(1995:1994+ye,A_tt(1,1:ye),'-b',1995:1994+ye,A_tt(2,1:ye),'-r',1995:1994+ye,A_tt(3,1:ye),'-g',1995:1994+ye,A_tt(4,1:ye),'-c',1995:1994+ye,A_tt(5,1:ye),'-k',1995:1994+ye,A_tt(7,1:ye),'-m',1995:1994+ye,A_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,NX_t(1,1:num_year),'--b',1995:1995+num_year-1,NX_t(2,1:num_year),'--r',1995:1995+num_year-1,NX_t(3,1:num_year),'--g',1995:1995+num_year-1,NX_t(4,1:num_year),'--c',1995:1995+num_year-1,NX_t(5,1:num_year),'--k',1995:1995+num_year-1,NX_t(7,1:num_year),'--m',1995:1995+num_year-1,NX_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Net transfer outflow')

sgtitle('Figure. Model fit for aggregate allocation')
saveas(gcf, 'plots_transition/16_Fit_agg_allocation.png');
close(gcf)

% sector prices
figure('Position', get(0, 'Screensize'))
subplot(3,4,1)
plot(1995:1994+ye,P_g_tt(1,1:ye),'-b',1995:1994+ye,P_g_tt(2,1:ye),'-r',1995:1994+ye,P_g_tt(3,1:ye),'-g',1995:1994+ye,P_g_tt(4,1:ye),'-c',1995:1994+ye,P_g_tt(5,1:ye),'-k',1995:1994+ye,P_g_tt(7,1:ye),'-m',1995:1994+ye,P_g_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,P_g_t(1,1:num_year),'--b',1995:1995+num_year-1,P_g_t(2,1:num_year),'--r',1995:1995+num_year-1,P_g_t(3,1:num_year),'--g',1995:1995+num_year-1,P_g_t(4,1:num_year),'--c',1995:1995+num_year-1,P_g_t(5,1:num_year),'--k',1995:1995+num_year-1,P_g_t(7,1:num_year),'--m',1995:1995+num_year-1,P_g_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Sector price Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 

subplot(3,4,2)
plot(1995:1994+ye,P_l_tt(1,1:ye),'-b',1995:1994+ye,P_l_tt(2,1:ye),'-r',1995:1994+ye,P_l_tt(3,1:ye),'-g',1995:1994+ye,P_l_tt(4,1:ye),'-c',1995:1994+ye,P_l_tt(5,1:ye),'-k',1995:1994+ye,P_l_tt(7,1:ye),'-m',1995:1994+ye,P_l_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,P_l_t(1,1:num_year),'--b',1995:1995+num_year-1,P_l_t(2,1:num_year),'--r',1995:1995+num_year-1,P_l_t(3,1:num_year),'--g',1995:1995+num_year-1,P_l_t(4,1:num_year),'--c',1995:1995+num_year-1,P_l_t(5,1:num_year),'--k',1995:1995+num_year-1,P_l_t(7,1:num_year),'--m',1995:1995+num_year-1,P_l_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Sector price Low-Skill Services')

subplot(3,4,3)
plot(1995:1994+ye,P_h_tt(1,1:ye),'-b',1995:1994+ye,P_h_tt(2,1:ye),'-r',1995:1994+ye,P_h_tt(3,1:ye),'-g',1995:1994+ye,P_h_tt(4,1:ye),'-c',1995:1994+ye,P_h_tt(5,1:ye),'-k',1995:1994+ye,P_h_tt(7,1:ye),'-m',1995:1994+ye,P_h_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,P_h_t(1,1:num_year),'--b',1995:1995+num_year-1,P_h_t(2,1:num_year),'--r',1995:1995+num_year-1,P_h_t(3,1:num_year),'--g',1995:1995+num_year-1,P_h_t(4,1:num_year),'--c',1995:1995+num_year-1,P_h_t(5,1:num_year),'--k',1995:1995+num_year-1,P_h_t(7,1:num_year),'--m',1995:1995+num_year-1,P_h_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Sector price High-Skill Services')

subplot(3,4,4)
plot(1995:1994+ye,w_tt(1,1:ye),'-b',1995:1994+ye,w_tt(2,1:ye),'-r',1995:1994+ye,w_tt(3,1:ye),'-g',1995:1994+ye,w_tt(4,1:ye),'-c',1995:1994+ye,w_tt(5,1:ye),'-k',1995:1994+ye,w_tt(7,1:ye),'-m',1995:1994+ye,w_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,w_t(1,1:num_year),'--b',1995:1995+num_year-1,w_t(2,1:num_year),'--r',1995:1995+num_year-1,w_t(3,1:num_year),'--g',1995:1995+num_year-1,w_t(4,1:num_year),'--c',1995:1995+num_year-1,w_t(5,1:num_year),'--k',1995:1995+num_year-1,w_t(7,1:num_year),'--m',1995:1995+num_year-1,w_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Wage')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 

subplot(3,4,5)
plot(1995:1994+ye,r_tt(1,1:ye),'-b',1995:1994+ye,r_tt(2,1:ye),'-r',1995:1994+ye,r_tt(3,1:ye),'-g',1995:1994+ye,r_tt(4,1:ye),'-c',1995:1994+ye,r_tt(5,1:ye),'-k',1995:1994+ye,r_tt(7,1:ye),'-m',1995:1994+ye,r_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,r_t(1,1:num_year),'--b',1995:1995+num_year-1,r_t(2,1:num_year),'--r',1995:1995+num_year-1,r_t(3,1:num_year),'--g',1995:1995+num_year-1,r_t(4,1:num_year),'--c',1995:1995+num_year-1,r_t(5,1:num_year),'--k',1995:1995+num_year-1,r_t(7,1:num_year),'--m',1995:1995+num_year-1,r_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Rental rate')

subplot(3,4,6)
plot(1995:1994+ye,P_c_tt(1,1:ye),'-b',1995:1994+ye,P_c_tt(2,1:ye),'-r',1995:1994+ye,P_c_tt(3,1:ye),'-g',1995:1994+ye,P_c_tt(4,1:ye),'-c',1995:1994+ye,P_c_tt(5,1:ye),'-k',1995:1994+ye,P_c_tt(7,1:ye),'-m',1995:1994+ye,P_c_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,P_c_t(1,1:num_year),'--b',1995:1995+num_year-1,P_c_t(2,1:num_year),'--r',1995:1995+num_year-1,P_c_t(3,1:num_year),'--g',1995:1995+num_year-1,P_c_t(4,1:num_year),'--c',1995:1995+num_year-1,P_c_t(5,1:num_year),'--k',1995:1995+num_year-1,P_c_t(7,1:num_year),'--m',1995:1995+num_year-1,P_c_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Consumption price')

subplot(3,4,7)
plot(1995:1994+ye,P_x_tt(1,1:ye),'-b',1995:1994+ye,P_x_tt(2,1:ye),'-r',1995:1994+ye,P_x_tt(3,1:ye),'-g',1995:1994+ye,P_x_tt(4,1:ye),'-c',1995:1994+ye,P_x_tt(5,1:ye),'-k',1995:1994+ye,P_x_tt(7,1:ye),'-m',1995:1994+ye,P_x_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,P_x_t(1,1:num_year),'--b',1995:1995+num_year-1,P_x_t(2,1:num_year),'--r',1995:1995+num_year-1,P_x_t(3,1:num_year),'--g',1995:1995+num_year-1,P_x_t(4,1:num_year),'--c',1995:1995+num_year-1,P_x_t(5,1:num_year),'--k',1995:1995+num_year-1,P_x_t(7,1:num_year),'--m',1995:1995+num_year-1,P_x_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Investment price')

subplot(3,4,8)
plot(1995:1994+ye,u_g_tt(1,1:ye),'-b',1995:1994+ye,u_g_tt(2,1:ye),'-r',1995:1994+ye,u_g_tt(3,1:ye),'-g',1995:1994+ye,u_g_tt(4,1:ye),'-c',1995:1994+ye,u_g_tt(5,1:ye),'-k',1995:1994+ye,u_g_tt(7,1:ye),'-m',1995:1994+ye,u_g_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Unit cost Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 

subplot(3,4,9)
plot(1995:1994+ye,u_l_tt(1,1:ye),'-b',1995:1994+ye,u_l_tt(2,1:ye),'-r',1995:1994+ye,u_l_tt(3,1:ye),'-g',1995:1994+ye,u_l_tt(4,1:ye),'-c',1995:1994+ye,u_l_tt(5,1:ye),'-k',1995:1994+ye,u_l_tt(7,1:ye),'-m',1995:1994+ye,u_l_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Unit cost Low-Skill Services')

subplot(3,4,10)
plot(1995:1994+ye,u_h_tt(1,1:ye),'-b',1995:1994+ye,u_h_tt(2,1:ye),'-r',1995:1994+ye,u_h_tt(3,1:ye),'-g',1995:1994+ye,u_h_tt(4,1:ye),'-c',1995:1994+ye,u_h_tt(5,1:ye),'-k',1995:1994+ye,u_h_tt(7,1:ye),'-m',1995:1994+ye,u_h_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Unit cost High-Skill Services')

sgtitle('Figure. Model fit for prices')
saveas(gcf, 'plots_transition/17_Fit_prices.png');
close(gcf)

% sector allocation
figure('Position', get(0, 'Screensize'))
subplot(3,3,1)
plot(1995:1994+ye,PY_g_tt(1,1:ye),'-b',1995:1994+ye,PY_g_tt(2,1:ye),'-r',1995:1994+ye,PY_g_tt(3,1:ye),'-g',1995:1994+ye,PY_g_tt(4,1:ye),'-c',1995:1994+ye,PY_g_tt(5,1:ye),'-k',1995:1994+ye,PY_g_tt(7,1:ye),'-m',1995:1994+ye,PY_g_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,PY_g_t(1,1:num_year),'--b',1995:1995+num_year-1,PY_g_t(2,1:num_year),'--r',1995:1995+num_year-1,PY_g_t(3,1:num_year),'--g',1995:1995+num_year-1,PY_g_t(4,1:num_year),'--c',1995:1995+num_year-1,PY_g_t(5,1:num_year),'--k',1995:1995+num_year-1,PY_g_t(7,1:num_year),'--m',1995:1995+num_year-1,PY_g_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Sector output Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 

subplot(3,3,2)
plot(1995:1994+ye,PY_l_tt(1,1:ye),'-b',1995:1994+ye,PY_l_tt(2,1:ye),'-r',1995:1994+ye,PY_l_tt(3,1:ye),'-g',1995:1994+ye,PY_l_tt(4,1:ye),'-c',1995:1994+ye,PY_l_tt(5,1:ye),'-k',1995:1994+ye,PY_l_tt(7,1:ye),'-m',1995:1994+ye,PY_l_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,PY_l_t(1,1:num_year),'--b',1995:1995+num_year-1,PY_l_t(2,1:num_year),'--r',1995:1995+num_year-1,PY_l_t(3,1:num_year),'--g',1995:1995+num_year-1,PY_l_t(4,1:num_year),'--c',1995:1995+num_year-1,PY_l_t(5,1:num_year),'--k',1995:1995+num_year-1,PY_l_t(7,1:num_year),'--m',1995:1995+num_year-1,PY_l_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Sector output Low-Skill Services')

subplot(3,3,3)
plot(1995:1994+ye,PY_h_tt(1,1:ye),'-b',1995:1994+ye,PY_h_tt(2,1:ye),'-r',1995:1994+ye,PY_h_tt(3,1:ye),'-g',1995:1994+ye,PY_h_tt(4,1:ye),'-c',1995:1994+ye,PY_h_tt(5,1:ye),'-k',1995:1994+ye,PY_h_tt(7,1:ye),'-m',1995:1994+ye,PY_h_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,PY_h_t(1,1:num_year),'--b',1995:1995+num_year-1,PY_h_t(2,1:num_year),'--r',1995:1995+num_year-1,PY_h_t(3,1:num_year),'--g',1995:1995+num_year-1,PY_h_t(4,1:num_year),'--c',1995:1995+num_year-1,PY_h_t(5,1:num_year),'--k',1995:1995+num_year-1,PY_h_t(7,1:num_year),'--m',1995:1995+num_year-1,PY_h_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Sector output High-Skill Services')

subplot(3,3,4)
plot(1995:1994+ye,PQ_g_tt(1,1:ye),'-b',1995:1994+ye,PQ_g_tt(2,1:ye),'-r',1995:1994+ye,PQ_g_tt(3,1:ye),'-g',1995:1994+ye,PQ_g_tt(4,1:ye),'-c',1995:1994+ye,PQ_g_tt(5,1:ye),'-k',1995:1994+ye,PQ_g_tt(7,1:ye),'-m',1995:1994+ye,PQ_g_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,PQ_g_t(1,1:num_year),'--b',1995:1995+num_year-1,PQ_g_t(2,1:num_year),'--r',1995:1995+num_year-1,PQ_g_t(3,1:num_year),'--g',1995:1995+num_year-1,PQ_g_t(4,1:num_year),'--c',1995:1995+num_year-1,PQ_g_t(5,1:num_year),'--k',1995:1995+num_year-1,PQ_g_t(7,1:num_year),'--m',1995:1995+num_year-1,PQ_g_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Sector absorprtion Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 

subplot(3,3,5)
plot(1995:1994+ye,PQ_l_tt(1,1:ye),'-b',1995:1994+ye,PQ_l_tt(2,1:ye),'-r',1995:1994+ye,PQ_l_tt(3,1:ye),'-g',1995:1994+ye,PQ_l_tt(4,1:ye),'-c',1995:1994+ye,PQ_l_tt(5,1:ye),'-k',1995:1994+ye,PQ_l_tt(7,1:ye),'-m',1995:1994+ye,PQ_l_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,PQ_l_t(1,1:num_year),'--b',1995:1995+num_year-1,PQ_l_t(2,1:num_year),'--r',1995:1995+num_year-1,PQ_l_t(3,1:num_year),'--g',1995:1995+num_year-1,PQ_l_t(4,1:num_year),'--c',1995:1995+num_year-1,PQ_l_t(5,1:num_year),'--k',1995:1995+num_year-1,PQ_l_t(7,1:num_year),'--m',1995:1995+num_year-1,PQ_l_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Sector absorption Low-Skill Services')

subplot(3,3,6)
plot(1995:1994+ye,PQ_h_tt(1,1:ye),'-b',1995:1994+ye,PQ_h_tt(2,1:ye),'-r',1995:1994+ye,PQ_h_tt(3,1:ye),'-g',1995:1994+ye,PQ_h_tt(4,1:ye),'-c',1995:1994+ye,PQ_h_tt(5,1:ye),'-k',1995:1994+ye,PQ_h_tt(7,1:ye),'-m',1995:1994+ye,PQ_h_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,PQ_h_t(1,1:num_year),'--b',1995:1995+num_year-1,PQ_h_t(2,1:num_year),'--r',1995:1995+num_year-1,PQ_h_t(3,1:num_year),'--g',1995:1995+num_year-1,PQ_h_t(4,1:num_year),'--c',1995:1995+num_year-1,PQ_h_t(5,1:num_year),'--k',1995:1995+num_year-1,PQ_h_t(7,1:num_year),'--m',1995:1995+num_year-1,PQ_h_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Sector absorption High-Skill Services')

subplot(3,3,7)
plot(1995:1994+ye,NX_g_tt(1,1:ye),'-b',1995:1994+ye,NX_g_tt(2,1:ye),'-r',1995:1994+ye,NX_g_tt(3,1:ye),'-g',1995:1994+ye,NX_g_tt(4,1:ye),'-c',1995:1994+ye,NX_g_tt(5,1:ye),'-k',1995:1994+ye,NX_g_tt(7,1:ye),'-m',1995:1994+ye,NX_g_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,NX_g_t(1,1:num_year),'--b',1995:1995+num_year-1,NX_g_t(2,1:num_year),'--r',1995:1995+num_year-1,NX_g_t(3,1:num_year),'--g',1995:1995+num_year-1,NX_g_t(4,1:num_year),'--c',1995:1995+num_year-1,NX_g_t(5,1:num_year),'--k',1995:1995+num_year-1,NX_g_t(7,1:num_year),'--m',1995:1995+num_year-1,NX_g_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Net export Agriculture')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 

subplot(3,3,8)
plot(1995:1994+ye,NX_l_tt(1,1:ye),'-b',1995:1994+ye,NX_l_tt(2,1:ye),'-r',1995:1994+ye,NX_l_tt(3,1:ye),'-g',1995:1994+ye,NX_l_tt(4,1:ye),'-c',1995:1994+ye,NX_l_tt(5,1:ye),'-k',1995:1994+ye,NX_l_tt(7,1:ye),'-m',1995:1994+ye,NX_l_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,NX_l_t(1,1:num_year),'--b',1995:1995+num_year-1,NX_l_t(2,1:num_year),'--r',1995:1995+num_year-1,NX_l_t(3,1:num_year),'--g',1995:1995+num_year-1,NX_l_t(4,1:num_year),'--c',1995:1995+num_year-1,NX_l_t(5,1:num_year),'--k',1995:1995+num_year-1,NX_l_t(7,1:num_year),'--m',1995:1995+num_year-1,NX_l_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Net export Low-Skill Services')

subplot(3,3,9)
plot(1995:1994+ye,NX_h_tt(1,1:ye),'-b',1995:1994+ye,NX_h_tt(2,1:ye),'-r',1995:1994+ye,NX_h_tt(3,1:ye),'-g',1995:1994+ye,NX_h_tt(4,1:ye),'-c',1995:1994+ye,NX_h_tt(5,1:ye),'-k',1995:1994+ye,NX_h_tt(7,1:ye),'-m',1995:1994+ye,NX_h_tt(9,1:ye),'-y', ...
    1995:1995+num_year-1,NX_h_t(1,1:num_year),'--b',1995:1995+num_year-1,NX_h_t(2,1:num_year),'--r',1995:1995+num_year-1,NX_h_t(3,1:num_year),'--g',1995:1995+num_year-1,NX_h_t(4,1:num_year),'--c',1995:1995+num_year-1,NX_h_t(5,1:num_year),'--k',1995:1995+num_year-1,NX_h_t(7,1:num_year),'--m',1995:1995+num_year-1,NX_h_t(9,1:num_year),'--y', ...
    'Linewidth',2.5)
title('Net export High-Skill Services')

sgtitle('Figure. Model fit for sector allocation')
saveas(gcf, 'plots_transition/18_Fit_sec_allocation.png');
close(gcf)

% sector bilateral trade share IND
figure('Position', get(0, 'Screensize'))
subplot(3,3,1)
plot(1995:1994+ye,squeeze(pi_g_tt(5,1,1:ye)),'-b',1995:1994+ye,squeeze(pi_g_tt(5,2,1:ye)),'-r',1995:1994+ye,squeeze(pi_g_tt(5,3,1:ye)),'-g',1995:1994+ye,squeeze(pi_g_tt(5,4,1:ye)),'-c',1995:1994+ye,squeeze(pi_g_tt(5,7,1:ye)),'-m',1995:1994+ye,squeeze(pi_g_tt(5,9,1:ye)),'-y', ...
    1995:1995+num_year-1,squeeze(bt_g_t(5,1,:)),'--b',1995:1995+num_year-1,squeeze(bt_g_t(5,2,:)),'--r',1995:1995+num_year-1,squeeze(bt_g_t(5,3,:)),'--g',1995:1995+num_year-1,squeeze(bt_g_t(5,4,:)),'--c',1995:1995+num_year-1,squeeze(bt_g_t(5,7,:)),'--m',1995:1995+num_year-1,squeeze(bt_g_t(5,9,:)),'--y', ...
    'Linewidth',2.5)
title('Bilateral trade for IND in Goods')
legend('AUS','BRA','CHN','DEU','MEX','USA','location','best') 

subplot(3,3,2)
plot(1995:1994+ye,squeeze(pi_l_tt(5,1,1:ye)),'-b',1995:1994+ye,squeeze(pi_l_tt(5,2,1:ye)),'-r',1995:1994+ye,squeeze(pi_l_tt(5,3,1:ye)),'-g',1995:1994+ye,squeeze(pi_l_tt(5,4,1:ye)),'-c',1995:1994+ye,squeeze(pi_l_tt(5,7,1:ye)),'-m',1995:1994+ye,squeeze(pi_l_tt(5,9,1:ye)),'-y', ...
    1995:1995+num_year-1,squeeze(bt_l_t(5,1,:)),'--b',1995:1995+num_year-1,squeeze(bt_l_t(5,2,:)),'--r',1995:1995+num_year-1,squeeze(bt_l_t(5,3,:)),'--g',1995:1995+num_year-1,squeeze(bt_l_t(5,4,:)),'--c',1995:1995+num_year-1,squeeze(bt_l_t(5,7,:)),'--m',1995:1995+num_year-1,squeeze(bt_l_t(5,9,:)),'--y', ...
    'Linewidth',2.5)
title('Bilateral trade for IND in Low-Skill Services')

subplot(3,3,3)
plot(1995:1994+ye,squeeze(pi_h_tt(5,1,1:ye)),'-b',1995:1994+ye,squeeze(pi_h_tt(5,2,1:ye)),'-r',1995:1994+ye,squeeze(pi_h_tt(5,3,1:ye)),'-g',1995:1994+ye,squeeze(pi_h_tt(5,4,1:ye)),'-c',1995:1994+ye,squeeze(pi_h_tt(5,7,1:ye)),'-m',1995:1994+ye,squeeze(pi_h_tt(5,9,1:ye)),'-y', ...
    1995:1995+num_year-1,squeeze(bt_h_t(5,1,:)),'--b',1995:1995+num_year-1,squeeze(bt_h_t(5,2,:)),'--r',1995:1995+num_year-1,squeeze(bt_h_t(5,3,:)),'--g',1995:1995+num_year-1,squeeze(bt_h_t(5,4,:)),'--c',1995:1995+num_year-1,squeeze(bt_h_t(5,7,:)),'--m',1995:1995+num_year-1,squeeze(bt_h_t(5,9,:)),'--y', ...
    'Linewidth',2.5)
title('Bilateral trade for IND in High-Skill Services')

subplot(3,3,4)
plot(1995:1994+ye,squeeze(pi_g_tt(3,1,1:ye)),'-b',1995:1994+ye,squeeze(pi_g_tt(3,2,1:ye)),'-r',1995:1994+ye,squeeze(pi_g_tt(3,4,1:ye)),'-c',1995:1994+ye,squeeze(pi_g_tt(3,5,1:ye)),'-k',1995:1994+ye,squeeze(pi_g_tt(3,7,1:ye)),'-m',1995:1994+ye,squeeze(pi_g_tt(3,9,1:ye)),'-y', ...
    1995:1995+num_year-1,squeeze(bt_g_t(3,1,:)),'--b',1995:1995+num_year-1,squeeze(bt_g_t(3,2,:)),'--r',1995:1995+num_year-1,squeeze(bt_g_t(3,4,:)),'--c',1995:1995+num_year-1,squeeze(bt_g_t(3,5,:)),'--k',1995:1995+num_year-1,squeeze(bt_g_t(3,7,:)),'--m',1995:1995+num_year-1,squeeze(bt_g_t(3,9,:)),'--y', ...
    'Linewidth',2.5)
title('Bilateral trade for CHN in Goods')
legend('AUS','BRA','DEU','IND','MEX','USA','location','best') 

subplot(3,3,5)
plot(1995:1994+ye,squeeze(pi_l_tt(3,1,1:ye)),'-b',1995:1994+ye,squeeze(pi_l_tt(3,2,1:ye)),'-r',1995:1994+ye,squeeze(pi_l_tt(3,4,1:ye)),'-c',1995:1994+ye,squeeze(pi_l_tt(3,5,1:ye)),'-k',1995:1994+ye,squeeze(pi_l_tt(3,7,1:ye)),'-m',1995:1994+ye,squeeze(pi_l_tt(3,9,1:ye)),'-y', ...
    1995:1995+num_year-1,squeeze(bt_l_t(3,1,:)),'--b',1995:1995+num_year-1,squeeze(bt_l_t(3,2,:)),'--r',1995:1995+num_year-1,squeeze(bt_l_t(3,4,:)),'--c',1995:1995+num_year-1,squeeze(bt_l_t(3,5,:)),'--k',1995:1995+num_year-1,squeeze(bt_l_t(3,7,:)),'--m',1995:1995+num_year-1,squeeze(bt_l_t(3,9,:)),'--y', ...
    'Linewidth',2.5)
title('Bilateral trade for CHN in Low-Skill Services')

subplot(3,3,6)
plot(1995:1994+ye,squeeze(pi_h_tt(3,1,1:ye)),'-b',1995:1994+ye,squeeze(pi_h_tt(3,2,1:ye)),'-r',1995:1994+ye,squeeze(pi_h_tt(3,4,1:ye)),'-c',1995:1994+ye,squeeze(pi_h_tt(3,5,1:ye)),'-k',1995:1994+ye,squeeze(pi_h_tt(3,7,1:ye)),'-m',1995:1994+ye,squeeze(pi_h_tt(3,9,1:ye)),'-y', ...
    1995:1995+num_year-1,squeeze(bt_h_t(3,1,:)),'--b',1995:1995+num_year-1,squeeze(bt_h_t(3,2,:)),'--r',1995:1995+num_year-1,squeeze(bt_h_t(3,4,:)),'--c',1995:1995+num_year-1,squeeze(bt_h_t(3,5,:)),'--k',1995:1995+num_year-1,squeeze(bt_h_t(3,7,:)),'--m',1995:1995+num_year-1,squeeze(bt_h_t(3,9,:)),'--y', ...
    'Linewidth',2.5)
title('Bilateral trade for CHN in High-Skill Services')

subplot(3,3,7)
plot(1995:1994+ye,squeeze(pi_g_tt(9,1,1:ye)),'-b',1995:1994+ye,squeeze(pi_g_tt(9,2,1:ye)),'-r',1995:1994+ye,squeeze(pi_g_tt(9,3,1:ye)),'-g',1995:1994+ye,squeeze(pi_g_tt(9,4,1:ye)),'-c',1995:1994+ye,squeeze(pi_g_tt(9,5,1:ye)),'-k',1995:1994+ye,squeeze(pi_g_tt(9,7,1:ye)),'-m', ...
    1995:1995+num_year-1,squeeze(bt_g_t(9,1,:)),'--b',1995:1995+num_year-1,squeeze(bt_g_t(9,2,:)),'--r',1995:1995+num_year-1,squeeze(bt_g_t(9,3,:)),'--g',1995:1995+num_year-1,squeeze(bt_g_t(9,4,:)),'--c',1995:1995+num_year-1,squeeze(bt_g_t(9,5,:)),'--k',1995:1995+num_year-1,squeeze(bt_g_t(9,7,:)),'--m', ...
    'Linewidth',2.5)
title('Bilateral trade for USA in Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','location','best') 

subplot(3,3,8)
plot(1995:1994+ye,squeeze(pi_l_tt(9,1,1:ye)),'-b',1995:1994+ye,squeeze(pi_l_tt(9,2,1:ye)),'-r',1995:1994+ye,squeeze(pi_l_tt(9,3,1:ye)),'-g',1995:1994+ye,squeeze(pi_l_tt(9,4,1:ye)),'-c',1995:1994+ye,squeeze(pi_l_tt(9,5,1:ye)),'-k',1995:1994+ye,squeeze(pi_l_tt(9,7,1:ye)),'-m', ...
    1995:1995+num_year-1,squeeze(bt_l_t(9,1,:)),'--b',1995:1995+num_year-1,squeeze(bt_l_t(9,2,:)),'--r',1995:1995+num_year-1,squeeze(bt_l_t(9,3,:)),'--g',1995:1995+num_year-1,squeeze(bt_l_t(9,4,:)),'--c',1995:1995+num_year-1,squeeze(bt_l_t(9,5,:)),'--k',1995:1995+num_year-1,squeeze(bt_l_t(9,7,:)),'--m', ...
    'Linewidth',2.5)
title('Bilateral trade for USA in Low-Skill Services')

subplot(3,3,9)
plot(1995:1994+ye,squeeze(pi_h_tt(9,1,1:ye)),'-b',1995:1994+ye,squeeze(pi_h_tt(9,2,1:ye)),'-r',1995:1994+ye,squeeze(pi_h_tt(9,3,1:ye)),'-g',1995:1994+ye,squeeze(pi_h_tt(9,4,1:ye)),'-c',1995:1994+ye,squeeze(pi_h_tt(9,5,1:ye)),'-k',1995:1994+ye,squeeze(pi_h_tt(9,7,1:ye)),'-m', ...
    1995:1995+num_year-1,squeeze(bt_h_t(9,1,:)),'--b',1995:1995+num_year-1,squeeze(bt_h_t(9,2,:)),'--r',1995:1995+num_year-1,squeeze(bt_h_t(9,3,:)),'--g',1995:1995+num_year-1,squeeze(bt_h_t(9,4,:)),'--c',1995:1995+num_year-1,squeeze(bt_h_t(9,5,:)),'--k',1995:1995+num_year-1,squeeze(bt_h_t(9,7,:)),'--m', ...
    'Linewidth',2.5)
title('Bilateral trade for USA in High-Skill Services')

sgtitle('Figure. Model fit for bilateral trade share')
saveas(gcf, 'plots_transition/19_Fit_trade_share.png');
close(gcf)


% demand shocks
figure('Position', get(0, 'Screensize'))
subplot(4,3,1)
plot(1995:1994+ye,T_g_tt(1,1:ye),'-b',1995:1994+ye,T_g_tt(2,1:ye),'-r',1995:1994+ye,T_g_tt(3,1:ye),'-g',1995:1994+ye,T_g_tt(4,1:ye),'-c',1995:1994+ye,T_g_tt(5,1:ye),'-k',1995:1994+ye,T_g_tt(7,1:ye),'-m',1995:1994+ye,T_g_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Productivity Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 

subplot(4,3,2)
plot(1995:1994+ye,T_l_tt(1,1:ye),'-b',1995:1994+ye,T_l_tt(2,1:ye),'-r',1995:1994+ye,T_l_tt(3,1:ye),'-g',1995:1994+ye,T_l_tt(4,1:ye),'-c',1995:1994+ye,T_l_tt(5,1:ye),'-k',1995:1994+ye,T_l_tt(7,1:ye),'-m',1995:1994+ye,T_l_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Productivity Low-Skill Services')

subplot(4,3,3)
plot(1995:1994+ye,T_h_tt(1,1:ye),'-b',1995:1994+ye,T_h_tt(2,1:ye),'-r',1995:1994+ye,T_h_tt(3,1:ye),'-g',1995:1994+ye,T_h_tt(4,1:ye),'-c',1995:1994+ye,T_h_tt(5,1:ye),'-k',1995:1994+ye,T_h_tt(7,1:ye),'-m',1995:1994+ye,T_h_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Productivity High-Skill Services')

subplot(4,3,4)
plot(1995:1994+ye,omega_cg_tt(1,1:ye),'-b',1995:1994+ye,omega_cg_tt(2,1:ye),'-r',1995:1994+ye,omega_cg_tt(3,1:ye),'-g',1995:1994+ye,omega_cg_tt(4,1:ye),'-c',1995:1994+ye,omega_cg_tt(5,1:ye),'-k',1995:1994+ye,omega_cg_tt(7,1:ye),'-m',1995:1994+ye,omega_cg_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Consumption share Goods')

subplot(4,3,5)
plot(1995:1994+ye,omega_cl_tt(1,1:ye),'-b',1995:1994+ye,omega_cl_tt(2,1:ye),'-r',1995:1994+ye,omega_cl_tt(3,1:ye),'-g',1995:1994+ye,omega_cl_tt(4,1:ye),'-c',1995:1994+ye,omega_cl_tt(5,1:ye),'-k',1995:1994+ye,omega_cl_tt(7,1:ye),'-m',1995:1994+ye,omega_cl_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Consumption share Low-Skill Services')

subplot(4,3,6)
plot(1995:1994+ye,omega_ch_tt(1,1:ye),'-b',1995:1994+ye,omega_ch_tt(2,1:ye),'-r',1995:1994+ye,omega_ch_tt(3,1:ye),'-g',1995:1994+ye,omega_ch_tt(4,1:ye),'-c',1995:1994+ye,omega_ch_tt(5,1:ye),'-k',1995:1994+ye,omega_ch_tt(7,1:ye),'-m',1995:1994+ye,omega_ch_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Consumption share High-Skill Services')
 
subplot(4,3,7)
plot(1995:1994+ye,omega_xg_tt(1,1:ye),'-b',1995:1994+ye,omega_xg_tt(2,1:ye),'-r',1995:1994+ye,omega_xg_tt(3,1:ye),'-g',1995:1994+ye,omega_xg_tt(4,1:ye),'-c',1995:1994+ye,omega_xg_tt(5,1:ye),'-k',1995:1994+ye,omega_xg_tt(7,1:ye),'-m',1995:1994+ye,omega_xg_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Investment share Goods')

subplot(4,3,8)
plot(1995:1994+ye,omega_xl_tt(1,1:ye),'-b',1995:1994+ye,omega_xl_tt(2,1:ye),'-r',1995:1994+ye,omega_xl_tt(3,1:ye),'-g',1995:1994+ye,omega_xl_tt(4,1:ye),'-c',1995:1994+ye,omega_xl_tt(5,1:ye),'-k',1995:1994+ye,omega_xl_tt(7,1:ye),'-m',1995:1994+ye,omega_xl_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Investment share Low-Skill Services')

subplot(4,3,9)
plot(1995:1994+ye,omega_xh_tt(1,1:ye),'-b',1995:1994+ye,omega_xh_tt(2,1:ye),'-r',1995:1994+ye,omega_xh_tt(3,1:ye),'-g',1995:1994+ye,omega_xh_tt(4,1:ye),'-c',1995:1994+ye,omega_xh_tt(5,1:ye),'-k',1995:1994+ye,omega_xh_tt(7,1:ye),'-m',1995:1994+ye,omega_xh_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Investment share High-Skill Services')

subplot(4,3,10)
plot(1995:1994+ye,PC_tt(1,1:ye),'-b',1995:1994+ye,PC_tt(2,1:ye),'-r',1995:1994+ye,PC_tt(3,1:ye),'-g',1995:1994+ye,PC_tt(4,1:ye),'-c',1995:1994+ye,PC_tt(5,1:ye),'-k',1995:1994+ye,PC_tt(7,1:ye),'-m',1995:1994+ye,PC_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Consumption PC as share of WPC')

subplot(4,3,11)
plot(1995:1994+ye,A_x_tt(1,1:ye),'-b',1995:1994+ye,A_x_tt(2,1:ye),'-r',1995:1994+ye,A_x_tt(3,1:ye),'-g',1995:1994+ye,A_x_tt(4,1:ye),'-c',1995:1994+ye,A_x_tt(5,1:ye),'-k',1995:1994+ye,A_x_tt(7,1:ye),'-m',1995:1994+ye,A_x_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Investment efficiency')

subplot(4,3,12)
plot(1995:1994+ye,L_tt(1,1:ye),'-b',1995:1994+ye,L_tt(2,1:ye),'-r',1995:1994+ye,L_tt(3,1:ye),'-g',1995:1994+ye,L_tt(4,1:ye),'-c',1995:1994+ye,L_tt(5,1:ye),'-k',1995:1994+ye,L_tt(7,1:ye),'-m',1995:1994+ye,L_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Employment')

sgtitle('Figure. Aggregate and sector demand shocks')
saveas(gcf, 'plots_transition/20_Shock_demand.png');

% relative shocks
figure('Position', get(0, 'Screensize'))
subplot(4,3,1)
plot(1995:1994+ye,T_g_tt(1,2:ye+1)./T_g_tt(1,1:ye),'-b',1995:1994+ye,T_g_tt(2,2:ye+1)./T_g_tt(2,1:ye),'-r',1995:1994+ye,T_g_tt(3,2:ye+1)./T_g_tt(3,1:ye),'-g',1995:1994+ye,T_g_tt(4,2:ye+1)./T_g_tt(4,1:ye),'-c',1995:1994+ye,T_g_tt(5,2:ye+1)./T_g_tt(5,1:ye),'-k',1995:1994+ye,T_g_tt(7,2:ye+1)./T_g_tt(7,1:ye),'-m',1995:1994+ye,T_g_tt(9,2:ye+1)./T_g_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Productivity Goods (time)')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','northeast') 

subplot(4,3,2)
plot(1995:1994+ye,T_l_tt(1,2:ye+1)./T_l_tt(1,1:ye),'-b',1995:1994+ye,T_l_tt(2,2:ye+1)./T_l_tt(2,1:ye),'-r',1995:1994+ye,T_l_tt(3,2:ye+1)./T_l_tt(3,1:ye),'-g',1995:1994+ye,T_l_tt(4,2:ye+1)./T_l_tt(4,1:ye),'-c',1995:1994+ye,T_l_tt(5,2:ye+1)./T_l_tt(5,1:ye),'-k',1995:1994+ye,T_l_tt(7,2:ye+1)./T_l_tt(7,1:ye),'-m',1995:1994+ye,T_l_tt(9,2:ye+1)./T_l_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Productivity Low-Skill Services')

subplot(4,3,3)
plot(1995:1994+ye,T_h_tt(1,2:ye+1)./T_h_tt(1,1:ye),'-b',1995:1994+ye,T_h_tt(2,2:ye+1)./T_h_tt(2,1:ye),'-r',1995:1994+ye,T_h_tt(3,2:ye+1)./T_h_tt(3,1:ye),'-g',1995:1994+ye,T_h_tt(4,2:ye+1)./T_h_tt(4,1:ye),'-c',1995:1994+ye,T_h_tt(5,2:ye+1)./T_h_tt(5,1:ye),'-k',1995:1994+ye,T_h_tt(7,2:ye+1)./T_h_tt(7,1:ye),'-m',1995:1994+ye,T_h_tt(9,2:ye+1)./T_h_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Productivity High-Skill Services')

subplot(4,3,4)
plot(1995:1994+ye,T_g_tt(1,1:ye)./T_g_tt(5,1:ye),'-b',1995:1994+ye,T_g_tt(2,1:ye)./T_g_tt(5,1:ye),'-r',1995:1994+ye,T_g_tt(3,1:ye)./T_g_tt(5,1:ye),'-g',1995:1994+ye,T_g_tt(4,1:ye)./T_g_tt(5,1:ye),'-c',1995:1994+ye,T_g_tt(5,1:ye)./T_g_tt(5,1:ye),'-k',1995:1994+ye,T_g_tt(7,1:ye)./T_g_tt(7,1:ye),'-m',1995:1994+ye,T_g_tt(9,1:ye)./T_g_tt(5,1:ye),'-y', ...
    'Linewidth',2.5)
title('Productivity Goods (country IND)')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','northeast') 

subplot(4,3,5)
plot(1995:1994+ye,T_l_tt(1,1:ye)./T_l_tt(5,1:ye),'-b',1995:1994+ye,T_l_tt(2,1:ye)./T_l_tt(5,1:ye),'-r',1995:1994+ye,T_l_tt(3,1:ye)./T_l_tt(5,1:ye),'-g',1995:1994+ye,T_l_tt(4,1:ye)./T_l_tt(5,1:ye),'-c',1995:1994+ye,T_l_tt(5,1:ye)./T_l_tt(5,1:ye),'-k',1995:1994+ye,T_l_tt(7,1:ye)./T_l_tt(7,1:ye),'-m',1995:1994+ye,T_l_tt(9,1:ye)./T_l_tt(5,1:ye),'-y', ...
    'Linewidth',2.5)
title('Productivity Low-Skill Services')

subplot(4,3,6)
plot(1995:1994+ye,T_h_tt(1,1:ye)./T_h_tt(5,1:ye),'-b',1995:1994+ye,T_h_tt(2,1:ye)./T_h_tt(5,1:ye),'-r',1995:1994+ye,T_h_tt(3,1:ye)./T_h_tt(5,1:ye),'-g',1995:1994+ye,T_h_tt(4,1:ye)./T_h_tt(5,1:ye),'-c',1995:1994+ye,T_h_tt(5,1:ye)./T_h_tt(5,1:ye),'-k',1995:1994+ye,T_h_tt(7,1:ye)./T_h_tt(7,1:ye),'-m',1995:1994+ye,T_h_tt(9,1:ye)./T_h_tt(5,1:ye),'-y', ...
    'Linewidth',2.5)
title('Productivity High-Skill Services')

subplot(4,3,7)
plot(1995:1994+ye,TFP_g_tt(1,1:ye)./TFP_g_tt(5,1:ye),'-b',1995:1994+ye,TFP_g_tt(2,1:ye)./TFP_g_tt(5,1:ye),'-r',1995:1994+ye,TFP_g_tt(3,1:ye)./TFP_g_tt(5,1:ye),'-g',1995:1994+ye,TFP_g_tt(4,1:ye)./TFP_g_tt(5,1:ye),'-c',1995:1994+ye,TFP_g_tt(5,1:ye)./TFP_g_tt(5,1:ye),'-k',1995:1994+ye,TFP_g_tt(7,1:ye)./TFP_g_tt(7,1:ye),'-m',1995:1994+ye,TFP_g_tt(9,1:ye)./TFP_g_tt(5,1:ye),'-y', ...
    'Linewidth',2.5)
title('Measured TFP Goods (country IND)')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','northeast') 

subplot(4,3,8)
plot(1995:1994+ye,TFP_l_tt(1,1:ye)./TFP_l_tt(5,1:ye),'-b',1995:1994+ye,TFP_l_tt(2,1:ye)./TFP_l_tt(5,1:ye),'-r',1995:1994+ye,TFP_l_tt(3,1:ye)./TFP_l_tt(5,1:ye),'-g',1995:1994+ye,TFP_l_tt(4,1:ye)./TFP_l_tt(5,1:ye),'-c',1995:1994+ye,TFP_l_tt(5,1:ye)./TFP_l_tt(5,1:ye),'-k',1995:1994+ye,TFP_l_tt(7,1:ye)./TFP_l_tt(7,1:ye),'-m',1995:1994+ye,TFP_l_tt(9,1:ye)./TFP_l_tt(5,1:ye),'-y', ...
    'Linewidth',2.5)
title('Measured TFP Low-Skill Services')

subplot(4,3,9)
plot(1995:1994+ye,TFP_h_tt(1,1:ye)./TFP_h_tt(5,1:ye),'-b',1995:1994+ye,TFP_h_tt(2,1:ye)./TFP_h_tt(5,1:ye),'-r',1995:1994+ye,TFP_h_tt(3,1:ye)./TFP_h_tt(5,1:ye),'-g',1995:1994+ye,TFP_h_tt(4,1:ye)./TFP_h_tt(5,1:ye),'-c',1995:1994+ye,TFP_h_tt(5,1:ye)./TFP_h_tt(5,1:ye),'-k',1995:1994+ye,TFP_h_tt(7,1:ye)./TFP_h_tt(7,1:ye),'-m',1995:1994+ye,TFP_h_tt(9,1:ye)./TFP_h_tt(5,1:ye),'-y', ...
    'Linewidth',2.5)
title('Measured TFP High-Skill Services')

subplot(4,3,10)
plot(1995:1994+ye,A_x_tt(1,2:ye+1)./A_x_tt(1,1:ye),'-b',1995:1994+ye,A_x_tt(2,2:ye+1)./A_x_tt(2,1:ye),'-r',1995:1994+ye,A_x_tt(3,2:ye+1)./A_x_tt(3,1:ye),'-g',1995:1994+ye,A_x_tt(4,2:ye+1)./A_x_tt(4,1:ye),'-c',1995:1994+ye,A_x_tt(5,2:ye+1)./A_x_tt(5,1:ye),'-k',1995:1994+ye,A_x_tt(7,2:ye+1)./A_x_tt(7,1:ye),'-m',1995:1994+ye,A_x_tt(9,2:ye+1)./A_x_tt(9,1:ye),'-y', ...
    'Linewidth',2.5)
title('Investment efficiency (time)')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','northeast') 

subplot(4,3,11)
plot(1995:1994+ye,A_x_tt(1,1:ye)./A_x_tt(5,1:ye),'-b',1995:1994+ye,A_x_tt(2,1:ye)./A_x_tt(5,1:ye),'-r',1995:1994+ye,A_x_tt(3,1:ye)./A_x_tt(5,1:ye),'-g',1995:1994+ye,A_x_tt(4,1:ye)./A_x_tt(5,1:ye),'-c',1995:1994+ye,A_x_tt(5,1:ye)./A_x_tt(5,1:ye),'-k',1995:1994+ye,A_x_tt(7,1:ye)./A_x_tt(7,1:ye),'-m',1995:1994+ye,A_x_tt(9,1:ye)./A_x_tt(5,1:ye),'-y', ...
    'Linewidth',2.5)
title('Investment efficiency (country IND)')

sgtitle('Figure. Relative shocks across time and country')
saveas(gcf, 'plots_transition/21_Relative_shock.png');



% trade shock
figure('Position', get(0, 'Screensize'))
subplot(4,3,1)
plot(1995:1994+ye,squeeze(d_g_tt(5,1,1:ye)),'-b',1995:1994+ye,squeeze(d_g_tt(5,2,1:ye)),'-r',1995:1994+ye,squeeze(d_g_tt(5,3,1:ye)),'-g',1995:1994+ye,squeeze(d_g_tt(5,4,1:ye)),'-c',1995:1994+ye,squeeze(d_g_tt(5,7,1:ye)),'-m',1995:1994+ye,squeeze(d_g_tt(5,9,1:ye)),'-y', ...
    1995:1994+ye,ones(I,ye),'--black','Linewidth',2.5)
title('Cost shipped to IND in Goods')
legend('AUS','BRA','CHN','DEU','MEX','USA','Value of 1','location','best') 

subplot(4,3,2)
plot(1995:1994+ye,squeeze(d_l_tt(5,1,1:ye)),'-b',1995:1994+ye,squeeze(d_l_tt(5,2,1:ye)),'-r',1995:1994+ye,squeeze(d_l_tt(5,3,1:ye)),'-g',1995:1994+ye,squeeze(d_l_tt(5,4,1:ye)),'-c',1995:1994+ye,squeeze(d_l_tt(5,7,1:ye)),'-m',1995:1994+ye,squeeze(d_l_tt(5,9,1:ye)),'-y', ...
    1995:1994+ye,ones(I,ye),'--black','Linewidth',2.5)
title('Cost shipped to IND in Low-skill services')

subplot(4,3,3)
plot(1995:1994+ye,squeeze(d_h_tt(5,1,1:ye)),'-b',1995:1994+ye,squeeze(d_h_tt(5,2,1:ye)),'-r',1995:1994+ye,squeeze(d_h_tt(5,3,1:ye)),'-g',1995:1994+ye,squeeze(d_h_tt(5,4,1:ye)),'-c',1995:1994+ye,squeeze(d_h_tt(5,7,1:ye)),'-m',1995:1994+ye,squeeze(d_h_tt(5,9,1:ye)),'-y', ...
    1995:1994+ye,ones(I,ye),'--black','Linewidth',2.5)
title('Cost shipped to IND in High-skill services')

subplot(4,3,4)
plot(1995:1994+ye,squeeze(d_g_tt(3,1,1:ye)),'-b',1995:1994+ye,squeeze(d_g_tt(3,2,1:ye)),'-r',1995:1994+ye,squeeze(d_g_tt(3,5,1:ye)),'-k',1995:1994+ye,squeeze(d_g_tt(3,4,1:ye)),'-c',1995:1994+ye,squeeze(d_g_tt(3,7,1:ye)),'-m',1995:1994+ye,squeeze(d_g_tt(3,9,1:ye)),'-y', ...
    1995:1994+ye,ones(I,ye),'--black','Linewidth',2.5)
title('Cost shipped to CHN in Goods')
legend('AUS','BRA','IND','DEU','MEX','USA','Value of 1','location','best') 

subplot(4,3,5)
plot(1995:1994+ye,squeeze(d_l_tt(3,1,1:ye)),'-b',1995:1994+ye,squeeze(d_l_tt(3,2,1:ye)),'-r',1995:1994+ye,squeeze(d_l_tt(3,5,1:ye)),'-k',1995:1994+ye,squeeze(d_l_tt(3,4,1:ye)),'-c',1995:1994+ye,squeeze(d_l_tt(3,7,1:ye)),'-m',1995:1994+ye,squeeze(d_l_tt(3,9,1:ye)),'-y', ...
    1995:1994+ye,ones(I,ye),'--black','Linewidth',2.5)
title('Cost shipped to CHN in Low-skill services')

subplot(4,3,6)
plot(1995:1994+ye,squeeze(d_h_tt(3,1,1:ye)),'-b',1995:1994+ye,squeeze(d_h_tt(3,2,1:ye)),'-r',1995:1994+ye,squeeze(d_h_tt(3,5,1:ye)),'-k',1995:1994+ye,squeeze(d_h_tt(3,4,1:ye)),'-c',1995:1994+ye,squeeze(d_h_tt(3,7,1:ye)),'-m',1995:1994+ye,squeeze(d_h_tt(3,9,1:ye)),'-y', ...
    1995:1994+ye,ones(I,ye),'--black','Linewidth',2.5)
title('Cost shipped to CHN in High-skill services')

subplot(4,3,7)
plot(1995:1994+ye,squeeze(d_g_tt(9,1,1:ye)),'-b',1995:1994+ye,squeeze(d_g_tt(9,2,1:ye)),'-r',1995:1994+ye,squeeze(d_g_tt(9,3,1:ye)),'-g',1995:1994+ye,squeeze(d_g_tt(9,4,1:ye)),'-c',1995:1994+ye,squeeze(d_g_tt(9,7,1:ye)),'-m',1995:1994+ye,squeeze(d_g_tt(9,5,1:ye)),'-k', ...
    1995:1994+ye,ones(I,ye),'--black','Linewidth',2.5)
title('Cost shipped to USA in Goods')
legend('AUS','BRA','CHN','DEU','MEX','IND','Value of 1','location','best') 

subplot(4,3,8)
plot(1995:1994+ye,squeeze(d_l_tt(9,1,1:ye)),'-b',1995:1994+ye,squeeze(d_l_tt(9,2,1:ye)),'-r',1995:1994+ye,squeeze(d_l_tt(9,3,1:ye)),'-g',1995:1994+ye,squeeze(d_l_tt(9,4,1:ye)),'-c',1995:1994+ye,squeeze(d_l_tt(9,7,1:ye)),'-m',1995:1994+ye,squeeze(d_l_tt(9,5,1:ye)),'-k', ...
    1995:1994+ye,ones(I,ye),'--black','Linewidth',2.5)
title('Cost shipped to USA in Low-skill services')

subplot(4,3,9)
plot(1995:1994+ye,squeeze(d_h_tt(9,1,1:ye)),'-b',1995:1994+ye,squeeze(d_h_tt(9,2,1:ye)),'-r',1995:1994+ye,squeeze(d_h_tt(9,3,1:ye)),'-g',1995:1994+ye,squeeze(d_h_tt(9,4,1:ye)),'-c',1995:1994+ye,squeeze(d_h_tt(9,7,1:ye)),'-m',1995:1994+ye,squeeze(d_h_tt(9,5,1:ye)),'-k', ...
    1995:1994+ye,ones(I,ye),'--black','Linewidth',2.5)
title('Cost shipped to USA in High-skill services')

subplot(4,3,10)
plot(1995:1994+ye,squeeze(d_g_tt(2,1,1:ye)),'-b',1995:1994+ye,squeeze(d_g_tt(2,5,1:ye)),'-k',1995:1994+ye,squeeze(d_g_tt(2,3,1:ye)),'-g',1995:1994+ye,squeeze(d_g_tt(2,4,1:ye)),'-c',1995:1994+ye,squeeze(d_g_tt(2,7,1:ye)),'-m',1995:1994+ye,squeeze(d_g_tt(2,9,1:ye)),'-y', ...
    1995:1994+ye,ones(I,ye),'--black','Linewidth',2.5)
title('Cost shipped to BRA in Goods')
legend('AUS','IND','CHN','DEU','MEX','USA','Value of 1','location','best') 

subplot(4,3,11)
plot(1995:1994+ye,squeeze(d_l_tt(2,1,1:ye)),'-b',1995:1994+ye,squeeze(d_l_tt(2,5,1:ye)),'-k',1995:1994+ye,squeeze(d_l_tt(2,3,1:ye)),'-g',1995:1994+ye,squeeze(d_l_tt(2,4,1:ye)),'-c',1995:1994+ye,squeeze(d_l_tt(2,7,1:ye)),'-m',1995:1994+ye,squeeze(d_l_tt(2,9,1:ye)),'-y', ...
    1995:1994+ye,ones(I,ye),'--black','Linewidth',2.5)
title('Cost shipped to BRA in Low-skill services')

subplot(4,3,12)
plot(1995:1994+ye,squeeze(d_h_tt(2,1,1:ye)),'-b',1995:1994+ye,squeeze(d_h_tt(2,5,1:ye)),'-k',1995:1994+ye,squeeze(d_h_tt(2,3,1:ye)),'-g',1995:1994+ye,squeeze(d_h_tt(2,4,1:ye)),'-c',1995:1994+ye,squeeze(d_h_tt(2,7,1:ye)),'-m',1995:1994+ye,squeeze(d_h_tt(2,9,1:ye)),'-y', ...
    1995:1994+ye,ones(I,ye),'--black','Linewidth',2.5)
title('Cost shipped to BRA in High-skill services')

sgtitle('Figure. Bilateral trade cost')
saveas(gcf, 'plots_transition/22_Shock_trade.png');
close(gcf)


%% Changing comparative advantage and intermediate use

% share of intermediate use over absorption
intm_g_tt = INTM_g_tt./PQ_g_tt; % share within each country is the same, due to same trade share for int and final use
intm_l_tt = INTM_l_tt./PQ_l_tt;
intm_h_tt = INTM_h_tt./PQ_h_tt;

% export and import level
EX_g_tt = PQ_g_tt*0;
EX_l_tt = PQ_l_tt*0;
EX_h_tt = PQ_h_tt*0;

IM_g_tt = PQ_g_tt*0;
IM_l_tt = PQ_l_tt*0;
IM_h_tt = PQ_h_tt*0;

for i=1:I
    
    EX_g_tt(i,:) = squeeze(pi_g_tt(1,i,:))'.*PQ_g_tt(1,:)+squeeze(pi_g_tt(2,i,:))'.*PQ_g_tt(2,:)+...
        squeeze(pi_g_tt(3,i,:))'.*PQ_g_tt(3,:);
    EX_g_tt(i,:) = EX_g_tt(i,:)-squeeze(pi_g_tt(i,i,:))'.*PQ_g_tt(i,:);
    
    EX_l_tt(i,:) = squeeze(pi_l_tt(1,i,:))'.*PQ_l_tt(1,:)+squeeze(pi_l_tt(2,i,:))'.*PQ_l_tt(2,:)+...
        squeeze(pi_l_tt(3,i,:))'.*PQ_l_tt(3,:);
    EX_l_tt(i,:) = EX_l_tt(i,:)-squeeze(pi_l_tt(i,i,:))'.*PQ_l_tt(i,:);
    
    EX_h_tt(i,:) = squeeze(pi_h_tt(1,i,:))'.*PQ_h_tt(1,:)+squeeze(pi_h_tt(2,i,:))'.*PQ_h_tt(2,:)+...
        squeeze(pi_h_tt(3,i,:))'.*PQ_h_tt(3,:);
    EX_h_tt(i,:) = EX_h_tt(i,:)-squeeze(pi_h_tt(i,i,:))'.*PQ_h_tt(i,:);
    
    IM_g_tt(i,:) = PQ_g_tt(i,:)-squeeze(pi_g_tt(i,i,:))'.*PQ_g_tt(i,:);
    IM_l_tt(i,:) = PQ_l_tt(i,:)-squeeze(pi_l_tt(i,i,:))'.*PQ_l_tt(i,:);
    IM_h_tt(i,:) = PQ_h_tt(i,:)-squeeze(pi_h_tt(i,i,:))'.*PQ_h_tt(i,:);     
end    

EX_tt = EX_g_tt+EX_l_tt+EX_h_tt;
IM_tt = IM_g_tt+IM_l_tt+IM_h_tt;


% bilateral export and import level
EX_1_tt = PQ_g_tt*0;
EX_2_tt = PQ_l_tt*0;
EX_3_tt = PQ_h_tt*0;
EX_4_tt = PQ_g_tt*0;
EX_5_tt = PQ_g_tt*0;
EX_6_tt = PQ_g_tt*0;
EX_7_tt = PQ_g_tt*0;
EX_8_tt = PQ_g_tt*0;
EX_9_tt = PQ_g_tt*0;

IM_1_tt = PQ_g_tt*0;
IM_2_tt = PQ_l_tt*0;
IM_3_tt = PQ_h_tt*0;
IM_4_tt = PQ_g_tt*0;
IM_5_tt = PQ_g_tt*0;
IM_6_tt = PQ_g_tt*0;
IM_7_tt = PQ_l_tt*0;
IM_8_tt = PQ_h_tt*0;
IM_9_tt = PQ_g_tt*0;

for i=1:I
    
    EX_1_tt(i,:) = squeeze(pi_g_tt(i,1,:))'.*PQ_g_tt(i,:)+squeeze(pi_l_tt(i,1,:))'.*PQ_l_tt(i,:)+...
        squeeze(pi_h_tt(i,1,:))'.*PQ_h_tt(i,:);
    EX_2_tt(i,:) = squeeze(pi_g_tt(i,2,:))'.*PQ_g_tt(i,:)+squeeze(pi_l_tt(i,2,:))'.*PQ_l_tt(i,:)+...
        squeeze(pi_h_tt(i,2,:))'.*PQ_h_tt(i,:);
    EX_3_tt(i,:) = squeeze(pi_g_tt(i,3,:))'.*PQ_g_tt(i,:)+squeeze(pi_l_tt(i,3,:))'.*PQ_l_tt(i,:)+...
        squeeze(pi_h_tt(i,3,:))'.*PQ_h_tt(i,:);
    EX_4_tt(i,:) = squeeze(pi_g_tt(i,4,:))'.*PQ_g_tt(i,:)+squeeze(pi_l_tt(i,4,:))'.*PQ_l_tt(i,:)+...
        squeeze(pi_h_tt(i,4,:))'.*PQ_h_tt(i,:);
    EX_5_tt(i,:) = squeeze(pi_g_tt(i,5,:))'.*PQ_g_tt(i,:)+squeeze(pi_l_tt(i,5,:))'.*PQ_l_tt(i,:)+...
        squeeze(pi_h_tt(i,5,:))'.*PQ_h_tt(i,:);
    EX_6_tt(i,:) = squeeze(pi_g_tt(i,6,:))'.*PQ_g_tt(i,:)+squeeze(pi_l_tt(i,6,:))'.*PQ_l_tt(i,:)+...
        squeeze(pi_h_tt(i,6,:))'.*PQ_h_tt(i,:);
    EX_7_tt(i,:) = squeeze(pi_g_tt(i,7,:))'.*PQ_g_tt(i,:)+squeeze(pi_l_tt(i,7,:))'.*PQ_l_tt(i,:)+...
        squeeze(pi_h_tt(i,7,:))'.*PQ_h_tt(i,:);
    EX_8_tt(i,:) = squeeze(pi_g_tt(i,8,:))'.*PQ_g_tt(i,:)+squeeze(pi_l_tt(i,8,:))'.*PQ_l_tt(i,:)+...
        squeeze(pi_h_tt(i,8,:))'.*PQ_h_tt(i,:);
    EX_9_tt(i,:) = squeeze(pi_g_tt(i,9,:))'.*PQ_g_tt(i,:)+squeeze(pi_l_tt(i,9,:))'.*PQ_l_tt(i,:)+...
        squeeze(pi_h_tt(i,9,:))'.*PQ_h_tt(i,:);
    
    IM_1_tt(i,:) = squeeze(pi_g_tt(1,i,:))'.*PQ_g_tt(1,:)+squeeze(pi_l_tt(1,i,:))'.*PQ_l_tt(1,:)+...
        squeeze(pi_h_tt(1,i,:))'.*PQ_h_tt(1,:);
    IM_2_tt(i,:) = squeeze(pi_g_tt(2,i,:))'.*PQ_g_tt(2,:)+squeeze(pi_l_tt(2,i,:))'.*PQ_l_tt(2,:)+...
        squeeze(pi_h_tt(2,i,:))'.*PQ_h_tt(2,:);
    IM_3_tt(i,:) = squeeze(pi_g_tt(3,i,:))'.*PQ_g_tt(3,:)+squeeze(pi_l_tt(3,i,:))'.*PQ_l_tt(3,:)+...
        squeeze(pi_h_tt(3,i,:))'.*PQ_h_tt(3,:);
    IM_4_tt(i,:) = squeeze(pi_g_tt(4,i,:))'.*PQ_g_tt(4,:)+squeeze(pi_l_tt(4,i,:))'.*PQ_l_tt(4,:)+...
        squeeze(pi_h_tt(4,i,:))'.*PQ_h_tt(4,:);   
    IM_5_tt(i,:) = squeeze(pi_g_tt(5,i,:))'.*PQ_g_tt(5,:)+squeeze(pi_l_tt(5,i,:))'.*PQ_l_tt(5,:)+...
        squeeze(pi_h_tt(5,i,:))'.*PQ_h_tt(5,:);
    IM_6_tt(i,:) = squeeze(pi_g_tt(6,i,:))'.*PQ_g_tt(6,:)+squeeze(pi_l_tt(6,i,:))'.*PQ_l_tt(6,:)+...
        squeeze(pi_h_tt(6,i,:))'.*PQ_h_tt(6,:);
    IM_7_tt(i,:) = squeeze(pi_g_tt(7,i,:))'.*PQ_g_tt(7,:)+squeeze(pi_l_tt(7,i,:))'.*PQ_l_tt(7,:)+...
        squeeze(pi_h_tt(7,i,:))'.*PQ_h_tt(7,:);
    IM_8_tt(i,:) = squeeze(pi_g_tt(8,i,:))'.*PQ_g_tt(8,:)+squeeze(pi_l_tt(8,i,:))'.*PQ_l_tt(8,:)+...
        squeeze(pi_h_tt(8,i,:))'.*PQ_h_tt(8,:);
    IM_9_tt(i,:) = squeeze(pi_g_tt(9,i,:))'.*PQ_g_tt(9,:)+squeeze(pi_l_tt(9,i,:))'.*PQ_l_tt(9,:)+...
        squeeze(pi_h_tt(9,i,:))'.*PQ_h_tt(9,:);
end

% bilateral export and import share
ex_1_2_tt = EX_1_tt(2,:)./EX_tt(1,:);
ex_1_3_tt = EX_1_tt(3,:)./EX_tt(1,:);
ex_1_4_tt = EX_1_tt(4,:)./EX_tt(1,:);
ex_1_5_tt = EX_1_tt(5,:)./EX_tt(1,:);
ex_1_6_tt = EX_1_tt(6,:)./EX_tt(1,:);
ex_1_7_tt = EX_1_tt(7,:)./EX_tt(1,:);
ex_1_8_tt = EX_1_tt(8,:)./EX_tt(1,:);
ex_1_9_tt = EX_1_tt(9,:)./EX_tt(1,:);

ex_2_1_tt = EX_2_tt(1,:)./EX_tt(2,:);
ex_2_3_tt = EX_2_tt(3,:)./EX_tt(2,:);
ex_2_4_tt = EX_2_tt(4,:)./EX_tt(2,:);
ex_2_5_tt = EX_2_tt(5,:)./EX_tt(2,:);
ex_2_6_tt = EX_2_tt(6,:)./EX_tt(2,:);
ex_2_7_tt = EX_2_tt(7,:)./EX_tt(2,:);
ex_2_8_tt = EX_2_tt(8,:)./EX_tt(2,:);
ex_2_9_tt = EX_2_tt(9,:)./EX_tt(2,:);

ex_3_1_tt = EX_3_tt(1,:)./EX_tt(3,:);
ex_3_2_tt = EX_3_tt(2,:)./EX_tt(3,:);
ex_3_4_tt = EX_3_tt(4,:)./EX_tt(3,:);
ex_3_5_tt = EX_3_tt(5,:)./EX_tt(3,:);
ex_3_6_tt = EX_3_tt(6,:)./EX_tt(3,:);
ex_3_7_tt = EX_3_tt(7,:)./EX_tt(3,:);
ex_3_8_tt = EX_3_tt(8,:)./EX_tt(3,:);
ex_3_9_tt = EX_3_tt(9,:)./EX_tt(3,:);

ex_4_1_tt = EX_4_tt(1,:)./EX_tt(4,:);
ex_4_2_tt = EX_4_tt(2,:)./EX_tt(4,:);
ex_4_3_tt = EX_4_tt(3,:)./EX_tt(4,:);
ex_4_5_tt = EX_4_tt(5,:)./EX_tt(4,:);
ex_4_6_tt = EX_4_tt(6,:)./EX_tt(4,:);
ex_4_7_tt = EX_4_tt(7,:)./EX_tt(4,:);
ex_4_8_tt = EX_4_tt(8,:)./EX_tt(4,:);
ex_4_9_tt = EX_4_tt(9,:)./EX_tt(4,:);

ex_5_1_tt = EX_5_tt(1,:)./EX_tt(5,:);
ex_5_2_tt = EX_5_tt(2,:)./EX_tt(5,:);
ex_5_3_tt = EX_5_tt(3,:)./EX_tt(5,:);
ex_5_4_tt = EX_5_tt(4,:)./EX_tt(5,:);
ex_5_6_tt = EX_5_tt(6,:)./EX_tt(5,:);
ex_5_7_tt = EX_5_tt(7,:)./EX_tt(5,:);
ex_5_8_tt = EX_5_tt(8,:)./EX_tt(5,:);
ex_5_9_tt = EX_5_tt(9,:)./EX_tt(5,:);

ex_6_1_tt = EX_6_tt(1,:)./EX_tt(6,:);
ex_6_2_tt = EX_6_tt(2,:)./EX_tt(6,:);
ex_6_3_tt = EX_6_tt(3,:)./EX_tt(6,:);
ex_6_4_tt = EX_6_tt(4,:)./EX_tt(6,:);
ex_6_5_tt = EX_6_tt(5,:)./EX_tt(6,:);
ex_6_7_tt = EX_6_tt(7,:)./EX_tt(6,:);
ex_6_8_tt = EX_6_tt(8,:)./EX_tt(6,:);
ex_6_9_tt = EX_6_tt(9,:)./EX_tt(6,:);

ex_7_1_tt = EX_7_tt(1,:)./EX_tt(7,:);
ex_7_2_tt = EX_7_tt(2,:)./EX_tt(7,:);
ex_7_3_tt = EX_7_tt(3,:)./EX_tt(7,:);
ex_7_4_tt = EX_7_tt(4,:)./EX_tt(7,:);
ex_7_5_tt = EX_7_tt(5,:)./EX_tt(7,:);
ex_7_6_tt = EX_7_tt(6,:)./EX_tt(7,:);
ex_7_8_tt = EX_7_tt(8,:)./EX_tt(7,:);
ex_7_9_tt = EX_7_tt(9,:)./EX_tt(7,:);

ex_8_1_tt = EX_8_tt(1,:)./EX_tt(8,:);
ex_8_2_tt = EX_8_tt(2,:)./EX_tt(8,:);
ex_8_3_tt = EX_8_tt(3,:)./EX_tt(8,:);
ex_8_4_tt = EX_8_tt(4,:)./EX_tt(8,:);
ex_8_5_tt = EX_8_tt(5,:)./EX_tt(8,:);
ex_8_6_tt = EX_8_tt(6,:)./EX_tt(8,:);
ex_8_7_tt = EX_8_tt(7,:)./EX_tt(8,:);
ex_8_9_tt = EX_8_tt(9,:)./EX_tt(8,:);

ex_9_1_tt = EX_9_tt(1,:)./EX_tt(9,:);
ex_9_2_tt = EX_9_tt(2,:)./EX_tt(9,:);
ex_9_3_tt = EX_9_tt(3,:)./EX_tt(9,:);
ex_9_4_tt = EX_9_tt(4,:)./EX_tt(9,:);
ex_9_5_tt = EX_9_tt(5,:)./EX_tt(9,:);
ex_9_6_tt = EX_9_tt(6,:)./EX_tt(9,:);
ex_9_7_tt = EX_9_tt(7,:)./EX_tt(9,:);
ex_9_8_tt = EX_9_tt(8,:)./EX_tt(9,:);


im_1_2_tt = IM_1_tt(2,:)./IM_tt(1,:);
im_1_3_tt = IM_1_tt(3,:)./IM_tt(1,:);
im_1_4_tt = IM_1_tt(4,:)./IM_tt(1,:);
im_1_5_tt = IM_1_tt(5,:)./IM_tt(1,:);
im_1_6_tt = IM_1_tt(6,:)./IM_tt(1,:);
im_1_7_tt = IM_1_tt(7,:)./IM_tt(1,:);
im_1_8_tt = IM_1_tt(8,:)./IM_tt(1,:);
im_1_9_tt = IM_1_tt(9,:)./IM_tt(1,:);

im_2_1_tt = IM_2_tt(1,:)./IM_tt(2,:);
im_2_3_tt = IM_2_tt(3,:)./IM_tt(2,:);
im_2_4_tt = IM_2_tt(4,:)./IM_tt(2,:);
im_2_5_tt = IM_2_tt(5,:)./IM_tt(2,:);
im_2_6_tt = IM_2_tt(6,:)./IM_tt(2,:);
im_2_7_tt = IM_2_tt(7,:)./IM_tt(2,:);
im_2_8_tt = IM_2_tt(8,:)./IM_tt(2,:);
im_2_9_tt = IM_2_tt(9,:)./IM_tt(2,:);

im_3_1_tt = IM_3_tt(1,:)./IM_tt(3,:);
im_3_2_tt = IM_3_tt(2,:)./IM_tt(3,:);
im_3_4_tt = IM_3_tt(4,:)./IM_tt(3,:);
im_3_5_tt = IM_3_tt(5,:)./IM_tt(3,:);
im_3_6_tt = IM_3_tt(6,:)./IM_tt(3,:);
im_3_7_tt = IM_3_tt(7,:)./IM_tt(3,:);
im_3_8_tt = IM_3_tt(8,:)./IM_tt(3,:);
im_3_9_tt = IM_3_tt(9,:)./IM_tt(3,:);

im_4_1_tt = IM_4_tt(1,:)./IM_tt(4,:);
im_4_2_tt = IM_4_tt(2,:)./IM_tt(4,:);
im_4_3_tt = IM_4_tt(3,:)./IM_tt(4,:);
im_4_5_tt = IM_4_tt(5,:)./IM_tt(4,:);
im_4_6_tt = IM_4_tt(6,:)./IM_tt(4,:);
im_4_7_tt = IM_4_tt(7,:)./IM_tt(4,:);
im_4_8_tt = IM_4_tt(8,:)./IM_tt(4,:);
im_4_9_tt = IM_4_tt(9,:)./IM_tt(4,:);

im_5_1_tt = IM_5_tt(1,:)./IM_tt(5,:);
im_5_2_tt = IM_5_tt(2,:)./IM_tt(5,:);
im_5_3_tt = IM_5_tt(3,:)./IM_tt(5,:);
im_5_4_tt = IM_5_tt(4,:)./IM_tt(5,:);
im_5_6_tt = IM_5_tt(6,:)./IM_tt(5,:);
im_5_7_tt = IM_5_tt(7,:)./IM_tt(5,:);
im_5_8_tt = IM_5_tt(8,:)./IM_tt(5,:);
im_5_9_tt = IM_5_tt(9,:)./IM_tt(5,:);

im_6_1_tt = IM_6_tt(1,:)./IM_tt(6,:);
im_6_2_tt = IM_6_tt(2,:)./IM_tt(6,:);
im_6_3_tt = IM_6_tt(3,:)./IM_tt(6,:);
im_6_4_tt = IM_6_tt(4,:)./IM_tt(6,:);
im_6_5_tt = IM_6_tt(5,:)./IM_tt(6,:);
im_6_7_tt = IM_6_tt(7,:)./IM_tt(6,:);
im_6_8_tt = IM_6_tt(8,:)./IM_tt(6,:);
im_6_9_tt = IM_6_tt(9,:)./IM_tt(6,:);

im_7_1_tt = IM_7_tt(1,:)./IM_tt(7,:);
im_7_2_tt = IM_7_tt(2,:)./IM_tt(7,:);
im_7_3_tt = IM_7_tt(3,:)./IM_tt(7,:);
im_7_4_tt = IM_7_tt(4,:)./IM_tt(7,:);
im_7_5_tt = IM_7_tt(5,:)./IM_tt(7,:);
im_7_6_tt = IM_7_tt(6,:)./IM_tt(7,:);
im_7_8_tt = IM_7_tt(8,:)./IM_tt(7,:);
im_7_9_tt = IM_7_tt(9,:)./IM_tt(7,:);

im_8_1_tt = IM_8_tt(1,:)./IM_tt(8,:);
im_8_2_tt = IM_8_tt(2,:)./IM_tt(8,:);
im_8_3_tt = IM_8_tt(3,:)./IM_tt(8,:);
im_8_4_tt = IM_8_tt(4,:)./IM_tt(8,:);
im_8_5_tt = IM_8_tt(5,:)./IM_tt(8,:);
im_8_6_tt = IM_8_tt(6,:)./IM_tt(8,:);
im_8_7_tt = IM_8_tt(7,:)./IM_tt(8,:);
im_8_9_tt = IM_8_tt(9,:)./IM_tt(8,:);

im_9_1_tt = IM_9_tt(1,:)./IM_tt(9,:);
im_9_2_tt = IM_9_tt(2,:)./IM_tt(9,:);
im_9_3_tt = IM_9_tt(3,:)./IM_tt(9,:);
im_9_4_tt = IM_9_tt(4,:)./IM_tt(9,:);
im_9_5_tt = IM_9_tt(5,:)./IM_tt(9,:);
im_9_6_tt = IM_9_tt(6,:)./IM_tt(9,:);
im_9_7_tt = IM_9_tt(7,:)./IM_tt(9,:);
im_9_8_tt = IM_9_tt(8,:)./IM_tt(9,:);


%% bilateral trade sectoral of India
ex_5_1_g_tt = squeeze(pi_g_tt(1,5,:))'.*PQ_g_tt(1,:)./EX_5_tt(1,:);
ex_5_1_l_tt = squeeze(pi_l_tt(1,5,:))'.*PQ_l_tt(1,:)./EX_5_tt(1,:);
ex_5_1_h_tt = squeeze(pi_h_tt(1,5,:))'.*PQ_h_tt(1,:)./EX_5_tt(1,:);

ex_5_2_g_tt = squeeze(pi_g_tt(2,5,:))'.*PQ_g_tt(2,:)./EX_5_tt(2,:);
ex_5_2_l_tt = squeeze(pi_l_tt(2,5,:))'.*PQ_l_tt(2,:)./EX_5_tt(2,:);
ex_5_2_h_tt = squeeze(pi_h_tt(2,5,:))'.*PQ_h_tt(2,:)./EX_5_tt(2,:);

ex_5_3_g_tt = squeeze(pi_g_tt(3,5,:))'.*PQ_g_tt(3,:)./EX_5_tt(3,:);
ex_5_3_l_tt = squeeze(pi_l_tt(3,5,:))'.*PQ_l_tt(3,:)./EX_5_tt(3,:);
ex_5_3_h_tt = squeeze(pi_h_tt(3,5,:))'.*PQ_h_tt(3,:)./EX_5_tt(3,:);

ex_5_4_g_tt = squeeze(pi_g_tt(4,5,:))'.*PQ_g_tt(4,:)./EX_5_tt(4,:);
ex_5_4_l_tt = squeeze(pi_l_tt(4,5,:))'.*PQ_l_tt(4,:)./EX_5_tt(4,:);
ex_5_4_h_tt = squeeze(pi_h_tt(4,5,:))'.*PQ_h_tt(4,:)./EX_5_tt(4,:);

ex_5_6_g_tt = squeeze(pi_g_tt(6,5,:))'.*PQ_g_tt(6,:)./EX_5_tt(6,:);
ex_5_6_l_tt = squeeze(pi_l_tt(6,5,:))'.*PQ_l_tt(6,:)./EX_5_tt(6,:);
ex_5_6_h_tt = squeeze(pi_h_tt(6,5,:))'.*PQ_h_tt(6,:)./EX_5_tt(6,:);

ex_5_7_g_tt = squeeze(pi_g_tt(7,5,:))'.*PQ_g_tt(7,:)./EX_5_tt(7,:);
ex_5_7_l_tt = squeeze(pi_l_tt(7,5,:))'.*PQ_l_tt(7,:)./EX_5_tt(7,:);
ex_5_7_h_tt = squeeze(pi_h_tt(7,5,:))'.*PQ_h_tt(7,:)./EX_5_tt(7,:);

ex_5_8_g_tt = squeeze(pi_g_tt(8,5,:))'.*PQ_g_tt(8,:)./EX_5_tt(8,:);
ex_5_8_l_tt = squeeze(pi_l_tt(8,5,:))'.*PQ_l_tt(8,:)./EX_5_tt(8,:);
ex_5_8_h_tt = squeeze(pi_h_tt(8,5,:))'.*PQ_h_tt(8,:)./EX_5_tt(8,:);

ex_5_9_g_tt = squeeze(pi_g_tt(9,5,:))'.*PQ_g_tt(9,:)./EX_5_tt(9,:);
ex_5_9_l_tt = squeeze(pi_l_tt(9,5,:))'.*PQ_l_tt(9,:)./EX_5_tt(9,:);
ex_5_9_h_tt = squeeze(pi_h_tt(9,5,:))'.*PQ_h_tt(9,:)./EX_5_tt(9,:);

im_5_1_g_tt = squeeze(pi_g_tt(5,1,:))'.*PQ_g_tt(5,:)./IM_5_tt(1,:);
im_5_1_l_tt = squeeze(pi_l_tt(5,1,:))'.*PQ_l_tt(5,:)./IM_5_tt(1,:);
im_5_1_h_tt = squeeze(pi_h_tt(5,1,:))'.*PQ_h_tt(5,:)./IM_5_tt(1,:);

im_5_2_g_tt = squeeze(pi_g_tt(5,2,:))'.*PQ_g_tt(5,:)./IM_5_tt(2,:);
im_5_2_l_tt = squeeze(pi_l_tt(5,2,:))'.*PQ_l_tt(5,:)./IM_5_tt(2,:);
im_5_2_h_tt = squeeze(pi_h_tt(5,2,:))'.*PQ_h_tt(5,:)./IM_5_tt(2,:);

im_5_3_g_tt = squeeze(pi_g_tt(5,3,:))'.*PQ_g_tt(5,:)./IM_5_tt(3,:);
im_5_3_l_tt = squeeze(pi_l_tt(5,3,:))'.*PQ_l_tt(5,:)./IM_5_tt(3,:);
im_5_3_h_tt = squeeze(pi_h_tt(5,3,:))'.*PQ_h_tt(5,:)./IM_5_tt(3,:);

im_5_4_g_tt = squeeze(pi_g_tt(5,4,:))'.*PQ_g_tt(5,:)./IM_5_tt(4,:);
im_5_4_l_tt = squeeze(pi_l_tt(5,4,:))'.*PQ_l_tt(5,:)./IM_5_tt(4,:);
im_5_4_h_tt = squeeze(pi_h_tt(5,4,:))'.*PQ_h_tt(5,:)./IM_5_tt(4,:);

im_5_6_g_tt = squeeze(pi_g_tt(5,6,:))'.*PQ_g_tt(5,:)./IM_5_tt(6,:);
im_5_6_l_tt = squeeze(pi_l_tt(5,6,:))'.*PQ_l_tt(5,:)./IM_5_tt(6,:);
im_5_6_h_tt = squeeze(pi_h_tt(5,6,:))'.*PQ_h_tt(5,:)./IM_5_tt(6,:);

im_5_7_g_tt = squeeze(pi_g_tt(5,7,:))'.*PQ_g_tt(5,:)./IM_5_tt(7,:);
im_5_7_l_tt = squeeze(pi_l_tt(5,7,:))'.*PQ_l_tt(5,:)./IM_5_tt(7,:);
im_5_7_h_tt = squeeze(pi_h_tt(5,7,:))'.*PQ_h_tt(5,:)./IM_5_tt(7,:);

im_5_8_g_tt = squeeze(pi_g_tt(5,8,:))'.*PQ_g_tt(5,:)./IM_5_tt(8,:);
im_5_8_l_tt = squeeze(pi_l_tt(5,8,:))'.*PQ_l_tt(5,:)./IM_5_tt(8,:);
im_5_8_h_tt = squeeze(pi_h_tt(5,8,:))'.*PQ_h_tt(5,:)./IM_5_tt(8,:);

im_5_9_g_tt = squeeze(pi_g_tt(5,9,:))'.*PQ_g_tt(5,:)./IM_5_tt(9,:);
im_5_9_l_tt = squeeze(pi_l_tt(5,9,:))'.*PQ_l_tt(5,:)./IM_5_tt(9,:);
im_5_9_h_tt = squeeze(pi_h_tt(5,9,:))'.*PQ_h_tt(5,:)./IM_5_tt(9,:);

ex_5_1_intm_tt = ( squeeze(pi_g_tt(1,5,:))'.*mu_gg_tt(1,:)+squeeze(pi_l_tt(1,5,:))'.*mu_gl_tt(1,:)+squeeze(pi_h_tt(1,5,:))'.*mu_gh_tt(1,:) ).*(1-nu_g_tt(1,:)).*PY_g_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,5,:))'.*mu_lg_tt(1,:)+squeeze(pi_l_tt(1,5,:))'.*mu_ll_tt(1,:)+squeeze(pi_h_tt(1,5,:))'.*mu_lh_tt(1,:) ).*(1-nu_l_tt(1,:)).*PY_l_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,5,:))'.*mu_hg_tt(1,:)+squeeze(pi_l_tt(1,5,:))'.*mu_hl_tt(1,:)+squeeze(pi_h_tt(1,5,:))'.*mu_hh_tt(1,:) ).*(1-nu_h_tt(1,:)).*PY_h_tt(1,:);
ex_5_1_intm_tt = ex_5_1_intm_tt./EX_5_tt(1,:);
ex_5_1_con_tt = ( squeeze(pi_g_tt(1,5,:))'.*omega_cg_tt(1,:)+squeeze(pi_l_tt(1,5,:))'.*omega_cl_tt(1,:)+squeeze(pi_h_tt(1,5,:))'.*omega_ch_tt(1,:) ).*PC_tt(1,:)./EX_5_tt(1,:);
ex_5_1_inv_tt = ( squeeze(pi_g_tt(1,5,:))'.*omega_xg_tt(1,:)+squeeze(pi_l_tt(1,5,:))'.*omega_xl_tt(1,:)+squeeze(pi_h_tt(1,5,:))'.*omega_xh_tt(1,:) ).*PX_tt(1,:)./EX_5_tt(1,:);

ex_5_2_intm_tt = ( squeeze(pi_g_tt(2,5,:))'.*mu_gg_tt(2,:)+squeeze(pi_l_tt(2,5,:))'.*mu_gl_tt(2,:)+squeeze(pi_h_tt(2,5,:))'.*mu_gh_tt(2,:) ).*(1-nu_g_tt(2,:)).*PY_g_tt(2,:)+... 
    ( squeeze(pi_g_tt(2,5,:))'.*mu_lg_tt(2,:)+squeeze(pi_l_tt(2,5,:))'.*mu_ll_tt(2,:)+squeeze(pi_h_tt(2,5,:))'.*mu_lh_tt(2,:) ).*(1-nu_l_tt(2,:)).*PY_l_tt(2,:)+... 
    ( squeeze(pi_g_tt(2,5,:))'.*mu_hg_tt(2,:)+squeeze(pi_l_tt(2,5,:))'.*mu_hl_tt(2,:)+squeeze(pi_h_tt(2,5,:))'.*mu_hh_tt(2,:) ).*(1-nu_h_tt(2,:)).*PY_h_tt(2,:);
ex_5_2_intm_tt = ex_5_2_intm_tt./EX_5_tt(2,:);
ex_5_2_con_tt = ( squeeze(pi_g_tt(2,5,:))'.*omega_cg_tt(2,:)+squeeze(pi_l_tt(2,5,:))'.*omega_cl_tt(2,:)+squeeze(pi_h_tt(2,5,:))'.*omega_ch_tt(2,:) ).*PC_tt(2,:)./EX_5_tt(2,:);
ex_5_2_inv_tt = ( squeeze(pi_g_tt(2,5,:))'.*omega_xg_tt(2,:)+squeeze(pi_l_tt(2,5,:))'.*omega_xl_tt(2,:)+squeeze(pi_h_tt(2,5,:))'.*omega_xh_tt(2,:) ).*PX_tt(2,:)./EX_5_tt(2,:);

ex_5_3_intm_tt = ( squeeze(pi_g_tt(3,5,:))'.*mu_gg_tt(3,:)+squeeze(pi_l_tt(3,5,:))'.*mu_gl_tt(3,:)+squeeze(pi_h_tt(3,5,:))'.*mu_gh_tt(3,:) ).*(1-nu_g_tt(3,:)).*PY_g_tt(3,:)+... 
    ( squeeze(pi_g_tt(3,5,:))'.*mu_lg_tt(3,:)+squeeze(pi_l_tt(3,5,:))'.*mu_ll_tt(3,:)+squeeze(pi_h_tt(3,5,:))'.*mu_lh_tt(3,:) ).*(1-nu_l_tt(3,:)).*PY_l_tt(3,:)+... 
    ( squeeze(pi_g_tt(3,5,:))'.*mu_hg_tt(3,:)+squeeze(pi_l_tt(3,5,:))'.*mu_hl_tt(3,:)+squeeze(pi_h_tt(3,5,:))'.*mu_hh_tt(3,:) ).*(1-nu_h_tt(3,:)).*PY_h_tt(3,:);
ex_5_3_intm_tt = ex_5_3_intm_tt./EX_5_tt(3,:);
ex_5_3_con_tt = ( squeeze(pi_g_tt(3,5,:))'.*omega_cg_tt(3,:)+squeeze(pi_l_tt(3,5,:))'.*omega_cl_tt(3,:)+squeeze(pi_h_tt(3,5,:))'.*omega_ch_tt(3,:) ).*PC_tt(3,:)./EX_5_tt(3,:);
ex_5_3_inv_tt = ( squeeze(pi_g_tt(3,5,:))'.*omega_xg_tt(3,:)+squeeze(pi_l_tt(3,5,:))'.*omega_xl_tt(3,:)+squeeze(pi_h_tt(3,5,:))'.*omega_xh_tt(3,:) ).*PX_tt(3,:)./EX_5_tt(3,:);

ex_5_4_intm_tt = ( squeeze(pi_g_tt(4,5,:))'.*mu_gg_tt(4,:)+squeeze(pi_l_tt(4,5,:))'.*mu_gl_tt(4,:)+squeeze(pi_h_tt(4,5,:))'.*mu_gh_tt(4,:) ).*(1-nu_g_tt(4,:)).*PY_g_tt(4,:)+... 
    ( squeeze(pi_g_tt(4,5,:))'.*mu_lg_tt(4,:)+squeeze(pi_l_tt(4,5,:))'.*mu_ll_tt(4,:)+squeeze(pi_h_tt(4,5,:))'.*mu_lh_tt(4,:) ).*(1-nu_l_tt(4,:)).*PY_l_tt(4,:)+... 
    ( squeeze(pi_g_tt(4,5,:))'.*mu_hg_tt(4,:)+squeeze(pi_l_tt(4,5,:))'.*mu_hl_tt(4,:)+squeeze(pi_h_tt(4,5,:))'.*mu_hh_tt(4,:) ).*(1-nu_h_tt(4,:)).*PY_h_tt(4,:);
ex_5_4_intm_tt = ex_5_4_intm_tt./EX_5_tt(4,:);
ex_5_4_con_tt = ( squeeze(pi_g_tt(4,5,:))'.*omega_cg_tt(4,:)+squeeze(pi_l_tt(4,5,:))'.*omega_cl_tt(4,:)+squeeze(pi_h_tt(4,5,:))'.*omega_ch_tt(4,:) ).*PC_tt(4,:)./EX_5_tt(4,:);
ex_5_4_inv_tt = ( squeeze(pi_g_tt(4,5,:))'.*omega_xg_tt(4,:)+squeeze(pi_l_tt(4,5,:))'.*omega_xl_tt(4,:)+squeeze(pi_h_tt(4,5,:))'.*omega_xh_tt(4,:) ).*PX_tt(4,:)./EX_5_tt(4,:);

ex_5_6_intm_tt = ( squeeze(pi_g_tt(6,5,:))'.*mu_gg_tt(6,:)+squeeze(pi_l_tt(6,5,:))'.*mu_gl_tt(6,:)+squeeze(pi_h_tt(6,5,:))'.*mu_gh_tt(6,:) ).*(1-nu_g_tt(6,:)).*PY_g_tt(6,:)+... 
    ( squeeze(pi_g_tt(6,5,:))'.*mu_lg_tt(6,:)+squeeze(pi_l_tt(6,5,:))'.*mu_ll_tt(6,:)+squeeze(pi_h_tt(6,5,:))'.*mu_lh_tt(6,:) ).*(1-nu_l_tt(6,:)).*PY_l_tt(6,:)+... 
    ( squeeze(pi_g_tt(6,5,:))'.*mu_hg_tt(6,:)+squeeze(pi_l_tt(6,5,:))'.*mu_hl_tt(6,:)+squeeze(pi_h_tt(6,5,:))'.*mu_hh_tt(6,:) ).*(1-nu_h_tt(6,:)).*PY_h_tt(6,:);
ex_5_6_intm_tt = ex_5_6_intm_tt./EX_5_tt(6,:);
ex_5_6_con_tt = ( squeeze(pi_g_tt(6,5,:))'.*omega_cg_tt(6,:)+squeeze(pi_l_tt(6,5,:))'.*omega_cl_tt(6,:)+squeeze(pi_h_tt(6,5,:))'.*omega_ch_tt(6,:) ).*PC_tt(6,:)./EX_5_tt(6,:);
ex_5_6_inv_tt = ( squeeze(pi_g_tt(6,5,:))'.*omega_xg_tt(6,:)+squeeze(pi_l_tt(6,5,:))'.*omega_xl_tt(6,:)+squeeze(pi_h_tt(6,5,:))'.*omega_xh_tt(6,:) ).*PX_tt(6,:)./EX_5_tt(6,:);

ex_5_7_intm_tt = ( squeeze(pi_g_tt(7,5,:))'.*mu_gg_tt(7,:)+squeeze(pi_l_tt(7,5,:))'.*mu_gl_tt(7,:)+squeeze(pi_h_tt(7,5,:))'.*mu_gh_tt(7,:) ).*(1-nu_g_tt(7,:)).*PY_g_tt(7,:)+... 
    ( squeeze(pi_g_tt(7,5,:))'.*mu_lg_tt(7,:)+squeeze(pi_l_tt(7,5,:))'.*mu_ll_tt(7,:)+squeeze(pi_h_tt(7,5,:))'.*mu_lh_tt(7,:) ).*(1-nu_l_tt(7,:)).*PY_l_tt(7,:)+... 
    ( squeeze(pi_g_tt(7,5,:))'.*mu_hg_tt(7,:)+squeeze(pi_l_tt(7,5,:))'.*mu_hl_tt(7,:)+squeeze(pi_h_tt(7,5,:))'.*mu_hh_tt(7,:) ).*(1-nu_h_tt(7,:)).*PY_h_tt(7,:);
ex_5_7_intm_tt = ex_5_7_intm_tt./EX_5_tt(7,:);
ex_5_7_con_tt = ( squeeze(pi_g_tt(7,5,:))'.*omega_cg_tt(7,:)+squeeze(pi_l_tt(7,5,:))'.*omega_cl_tt(7,:)+squeeze(pi_h_tt(7,5,:))'.*omega_ch_tt(7,:) ).*PC_tt(7,:)./EX_5_tt(7,:);
ex_5_7_inv_tt = ( squeeze(pi_g_tt(7,5,:))'.*omega_xg_tt(7,:)+squeeze(pi_l_tt(7,5,:))'.*omega_xl_tt(7,:)+squeeze(pi_h_tt(7,5,:))'.*omega_xh_tt(7,:) ).*PX_tt(7,:)./EX_5_tt(7,:);

ex_5_8_intm_tt = ( squeeze(pi_g_tt(8,5,:))'.*mu_gg_tt(8,:)+squeeze(pi_l_tt(8,5,:))'.*mu_gl_tt(8,:)+squeeze(pi_h_tt(8,5,:))'.*mu_gh_tt(8,:) ).*(1-nu_g_tt(8,:)).*PY_g_tt(8,:)+... 
    ( squeeze(pi_g_tt(8,5,:))'.*mu_lg_tt(8,:)+squeeze(pi_l_tt(8,5,:))'.*mu_ll_tt(8,:)+squeeze(pi_h_tt(8,5,:))'.*mu_lh_tt(8,:) ).*(1-nu_l_tt(8,:)).*PY_l_tt(8,:)+... 
    ( squeeze(pi_g_tt(8,5,:))'.*mu_hg_tt(8,:)+squeeze(pi_l_tt(8,5,:))'.*mu_hl_tt(8,:)+squeeze(pi_h_tt(8,5,:))'.*mu_hh_tt(8,:) ).*(1-nu_h_tt(8,:)).*PY_h_tt(8,:);
ex_5_8_intm_tt = ex_5_8_intm_tt./EX_5_tt(8,:);
ex_5_8_con_tt = ( squeeze(pi_g_tt(8,5,:))'.*omega_cg_tt(8,:)+squeeze(pi_l_tt(8,5,:))'.*omega_cl_tt(8,:)+squeeze(pi_h_tt(8,5,:))'.*omega_ch_tt(8,:) ).*PC_tt(8,:)./EX_5_tt(8,:);
ex_5_8_inv_tt = ( squeeze(pi_g_tt(8,5,:))'.*omega_xg_tt(8,:)+squeeze(pi_l_tt(8,5,:))'.*omega_xl_tt(8,:)+squeeze(pi_h_tt(8,5,:))'.*omega_xh_tt(8,:) ).*PX_tt(8,:)./EX_5_tt(8,:);

ex_5_9_intm_tt = ( squeeze(pi_g_tt(9,5,:))'.*mu_gg_tt(9,:)+squeeze(pi_l_tt(9,5,:))'.*mu_gl_tt(9,:)+squeeze(pi_h_tt(9,5,:))'.*mu_gh_tt(9,:) ).*(1-nu_g_tt(9,:)).*PY_g_tt(9,:)+... 
    ( squeeze(pi_g_tt(9,5,:))'.*mu_lg_tt(9,:)+squeeze(pi_l_tt(9,5,:))'.*mu_ll_tt(9,:)+squeeze(pi_h_tt(9,5,:))'.*mu_lh_tt(9,:) ).*(1-nu_l_tt(9,:)).*PY_l_tt(9,:)+... 
    ( squeeze(pi_g_tt(9,5,:))'.*mu_hg_tt(9,:)+squeeze(pi_l_tt(9,5,:))'.*mu_hl_tt(9,:)+squeeze(pi_h_tt(9,5,:))'.*mu_hh_tt(9,:) ).*(1-nu_h_tt(9,:)).*PY_h_tt(9,:);
ex_5_9_intm_tt = ex_5_9_intm_tt./EX_5_tt(9,:);
ex_5_9_con_tt = ( squeeze(pi_g_tt(9,5,:))'.*omega_cg_tt(9,:)+squeeze(pi_l_tt(9,5,:))'.*omega_cl_tt(9,:)+squeeze(pi_h_tt(9,5,:))'.*omega_ch_tt(9,:) ).*PC_tt(9,:)./EX_5_tt(9,:);
ex_5_9_inv_tt = ( squeeze(pi_g_tt(9,5,:))'.*omega_xg_tt(9,:)+squeeze(pi_l_tt(9,5,:))'.*omega_xl_tt(9,:)+squeeze(pi_h_tt(9,5,:))'.*omega_xh_tt(9,:) ).*PX_tt(9,:)./EX_5_tt(9,:);


im_1_2_intm_tt = ( squeeze(pi_g_tt(1,2,:))'.*mu_gg_tt(1,:)+squeeze(pi_l_tt(1,2,:))'.*mu_gl_tt(1,:)+squeeze(pi_h_tt(1,2,:))'.*mu_gh_tt(1,:) ).*(1-nu_g_tt(1,:)).*PY_g_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,2,:))'.*mu_lg_tt(1,:)+squeeze(pi_l_tt(1,2,:))'.*mu_ll_tt(1,:)+squeeze(pi_h_tt(1,2,:))'.*mu_lh_tt(1,:) ).*(1-nu_l_tt(1,:)).*PY_l_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,2,:))'.*mu_hg_tt(1,:)+squeeze(pi_l_tt(1,2,:))'.*mu_hl_tt(1,:)+squeeze(pi_h_tt(1,2,:))'.*mu_hh_tt(1,:) ).*(1-nu_h_tt(1,:)).*PY_h_tt(1,:);
im_1_2_intm_tt = im_1_2_intm_tt./IM_1_tt(2,:);
im_1_2_con_tt = ( squeeze(pi_g_tt(1,2,:))'.*omega_cg_tt(1,:)+squeeze(pi_l_tt(1,2,:))'.*omega_cl_tt(1,:)+squeeze(pi_h_tt(1,2,:))'.*omega_ch_tt(1,:) ).*PC_tt(1,:)./IM_1_tt(2,:);
im_1_2_inv_tt = ( squeeze(pi_g_tt(1,2,:))'.*omega_xg_tt(1,:)+squeeze(pi_l_tt(1,2,:))'.*omega_xl_tt(1,:)+squeeze(pi_h_tt(1,2,:))'.*omega_xh_tt(1,:) ).*PX_tt(1,:)./IM_1_tt(2,:);

im_1_3_intm_tt = ( squeeze(pi_g_tt(1,3,:))'.*mu_gg_tt(1,:)+squeeze(pi_l_tt(1,3,:))'.*mu_gl_tt(1,:)+squeeze(pi_h_tt(1,3,:))'.*mu_gh_tt(1,:) ).*(1-nu_g_tt(1,:)).*PY_g_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,3,:))'.*mu_lg_tt(1,:)+squeeze(pi_l_tt(1,3,:))'.*mu_ll_tt(1,:)+squeeze(pi_h_tt(1,3,:))'.*mu_lh_tt(1,:) ).*(1-nu_l_tt(1,:)).*PY_l_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,3,:))'.*mu_hg_tt(1,:)+squeeze(pi_l_tt(1,3,:))'.*mu_hl_tt(1,:)+squeeze(pi_h_tt(1,3,:))'.*mu_hh_tt(1,:) ).*(1-nu_h_tt(1,:)).*PY_h_tt(1,:);
im_1_3_intm_tt = im_1_3_intm_tt./IM_1_tt(3,:);
im_1_3_con_tt = ( squeeze(pi_g_tt(1,3,:))'.*omega_cg_tt(1,:)+squeeze(pi_l_tt(1,3,:))'.*omega_cl_tt(1,:)+squeeze(pi_h_tt(1,3,:))'.*omega_ch_tt(1,:) ).*PC_tt(1,:)./IM_1_tt(3,:);
im_1_3_inv_tt = ( squeeze(pi_g_tt(1,3,:))'.*omega_xg_tt(1,:)+squeeze(pi_l_tt(1,3,:))'.*omega_xl_tt(1,:)+squeeze(pi_h_tt(1,3,:))'.*omega_xh_tt(1,:) ).*PX_tt(1,:)./IM_1_tt(3,:);

im_1_4_intm_tt = ( squeeze(pi_g_tt(1,4,:))'.*mu_gg_tt(1,:)+squeeze(pi_l_tt(1,4,:))'.*mu_gl_tt(1,:)+squeeze(pi_h_tt(1,4,:))'.*mu_gh_tt(1,:) ).*(1-nu_g_tt(1,:)).*PY_g_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,4,:))'.*mu_lg_tt(1,:)+squeeze(pi_l_tt(1,4,:))'.*mu_ll_tt(1,:)+squeeze(pi_h_tt(1,4,:))'.*mu_lh_tt(1,:) ).*(1-nu_l_tt(1,:)).*PY_l_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,4,:))'.*mu_hg_tt(1,:)+squeeze(pi_l_tt(1,4,:))'.*mu_hl_tt(1,:)+squeeze(pi_h_tt(1,4,:))'.*mu_hh_tt(1,:) ).*(1-nu_h_tt(1,:)).*PY_h_tt(1,:);
im_1_4_intm_tt = im_1_4_intm_tt./IM_1_tt(4,:);
im_1_4_con_tt = ( squeeze(pi_g_tt(1,4,:))'.*omega_cg_tt(1,:)+squeeze(pi_l_tt(1,4,:))'.*omega_cl_tt(1,:)+squeeze(pi_h_tt(1,4,:))'.*omega_ch_tt(1,:) ).*PC_tt(1,:)./IM_1_tt(4,:);
im_1_4_inv_tt = ( squeeze(pi_g_tt(1,4,:))'.*omega_xg_tt(1,:)+squeeze(pi_l_tt(1,4,:))'.*omega_xl_tt(1,:)+squeeze(pi_h_tt(1,4,:))'.*omega_xh_tt(1,:) ).*PX_tt(1,:)./IM_1_tt(4,:);

im_1_5_intm_tt = ( squeeze(pi_g_tt(1,5,:))'.*mu_gg_tt(1,:)+squeeze(pi_l_tt(1,5,:))'.*mu_gl_tt(1,:)+squeeze(pi_h_tt(1,5,:))'.*mu_gh_tt(1,:) ).*(1-nu_g_tt(1,:)).*PY_g_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,5,:))'.*mu_lg_tt(1,:)+squeeze(pi_l_tt(1,5,:))'.*mu_ll_tt(1,:)+squeeze(pi_h_tt(1,5,:))'.*mu_lh_tt(1,:) ).*(1-nu_l_tt(1,:)).*PY_l_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,5,:))'.*mu_hg_tt(1,:)+squeeze(pi_l_tt(1,5,:))'.*mu_hl_tt(1,:)+squeeze(pi_h_tt(1,5,:))'.*mu_hh_tt(1,:) ).*(1-nu_h_tt(1,:)).*PY_h_tt(1,:);
im_1_5_intm_tt = im_1_5_intm_tt./IM_1_tt(5,:);
im_1_5_con_tt = ( squeeze(pi_g_tt(1,5,:))'.*omega_cg_tt(1,:)+squeeze(pi_l_tt(1,5,:))'.*omega_cl_tt(1,:)+squeeze(pi_h_tt(1,5,:))'.*omega_ch_tt(1,:) ).*PC_tt(1,:)./IM_1_tt(5,:);
im_1_5_inv_tt = ( squeeze(pi_g_tt(1,5,:))'.*omega_xg_tt(1,:)+squeeze(pi_l_tt(1,5,:))'.*omega_xl_tt(1,:)+squeeze(pi_h_tt(1,5,:))'.*omega_xh_tt(1,:) ).*PX_tt(1,:)./IM_1_tt(5,:);

im_1_6_intm_tt = ( squeeze(pi_g_tt(1,6,:))'.*mu_gg_tt(1,:)+squeeze(pi_l_tt(1,6,:))'.*mu_gl_tt(1,:)+squeeze(pi_h_tt(1,6,:))'.*mu_gh_tt(1,:) ).*(1-nu_g_tt(1,:)).*PY_g_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,6,:))'.*mu_lg_tt(1,:)+squeeze(pi_l_tt(1,6,:))'.*mu_ll_tt(1,:)+squeeze(pi_h_tt(1,6,:))'.*mu_lh_tt(1,:) ).*(1-nu_l_tt(1,:)).*PY_l_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,6,:))'.*mu_hg_tt(1,:)+squeeze(pi_l_tt(1,6,:))'.*mu_hl_tt(1,:)+squeeze(pi_h_tt(1,6,:))'.*mu_hh_tt(1,:) ).*(1-nu_h_tt(1,:)).*PY_h_tt(1,:);
im_1_6_intm_tt = im_1_6_intm_tt./IM_1_tt(6,:);
im_1_6_con_tt = ( squeeze(pi_g_tt(1,6,:))'.*omega_cg_tt(1,:)+squeeze(pi_l_tt(1,6,:))'.*omega_cl_tt(1,:)+squeeze(pi_h_tt(1,6,:))'.*omega_ch_tt(1,:) ).*PC_tt(1,:)./IM_1_tt(6,:);
im_1_6_inv_tt = ( squeeze(pi_g_tt(1,6,:))'.*omega_xg_tt(1,:)+squeeze(pi_l_tt(1,6,:))'.*omega_xl_tt(1,:)+squeeze(pi_h_tt(1,6,:))'.*omega_xh_tt(1,:) ).*PX_tt(1,:)./IM_1_tt(6,:);

im_1_7_intm_tt = ( squeeze(pi_g_tt(1,7,:))'.*mu_gg_tt(1,:)+squeeze(pi_l_tt(1,7,:))'.*mu_gl_tt(1,:)+squeeze(pi_h_tt(1,7,:))'.*mu_gh_tt(1,:) ).*(1-nu_g_tt(1,:)).*PY_g_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,7,:))'.*mu_lg_tt(1,:)+squeeze(pi_l_tt(1,7,:))'.*mu_ll_tt(1,:)+squeeze(pi_h_tt(1,7,:))'.*mu_lh_tt(1,:) ).*(1-nu_l_tt(1,:)).*PY_l_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,7,:))'.*mu_hg_tt(1,:)+squeeze(pi_l_tt(1,7,:))'.*mu_hl_tt(1,:)+squeeze(pi_h_tt(1,7,:))'.*mu_hh_tt(1,:) ).*(1-nu_h_tt(1,:)).*PY_h_tt(1,:);
im_1_7_intm_tt = im_1_7_intm_tt./IM_1_tt(7,:);
im_1_7_con_tt = ( squeeze(pi_g_tt(1,7,:))'.*omega_cg_tt(1,:)+squeeze(pi_l_tt(1,7,:))'.*omega_cl_tt(1,:)+squeeze(pi_h_tt(1,7,:))'.*omega_ch_tt(1,:) ).*PC_tt(1,:)./IM_1_tt(7,:);
im_1_7_inv_tt = ( squeeze(pi_g_tt(1,7,:))'.*omega_xg_tt(1,:)+squeeze(pi_l_tt(1,7,:))'.*omega_xl_tt(1,:)+squeeze(pi_h_tt(1,7,:))'.*omega_xh_tt(1,:) ).*PX_tt(1,:)./IM_1_tt(7,:);

im_1_8_intm_tt = ( squeeze(pi_g_tt(1,8,:))'.*mu_gg_tt(1,:)+squeeze(pi_l_tt(1,8,:))'.*mu_gl_tt(1,:)+squeeze(pi_h_tt(1,8,:))'.*mu_gh_tt(1,:) ).*(1-nu_g_tt(1,:)).*PY_g_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,8,:))'.*mu_lg_tt(1,:)+squeeze(pi_l_tt(1,8,:))'.*mu_ll_tt(1,:)+squeeze(pi_h_tt(1,8,:))'.*mu_lh_tt(1,:) ).*(1-nu_l_tt(1,:)).*PY_l_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,8,:))'.*mu_hg_tt(1,:)+squeeze(pi_l_tt(1,8,:))'.*mu_hl_tt(1,:)+squeeze(pi_h_tt(1,8,:))'.*mu_hh_tt(1,:) ).*(1-nu_h_tt(1,:)).*PY_h_tt(1,:);
im_1_8_intm_tt = im_1_8_intm_tt./IM_1_tt(8,:);
im_1_8_con_tt = ( squeeze(pi_g_tt(1,8,:))'.*omega_cg_tt(1,:)+squeeze(pi_l_tt(1,8,:))'.*omega_cl_tt(1,:)+squeeze(pi_h_tt(1,8,:))'.*omega_ch_tt(1,:) ).*PC_tt(1,:)./IM_1_tt(8,:);
im_1_8_inv_tt = ( squeeze(pi_g_tt(1,8,:))'.*omega_xg_tt(1,:)+squeeze(pi_l_tt(1,8,:))'.*omega_xl_tt(1,:)+squeeze(pi_h_tt(1,8,:))'.*omega_xh_tt(1,:) ).*PX_tt(1,:)./IM_1_tt(8,:);

im_1_9_intm_tt = ( squeeze(pi_g_tt(1,9,:))'.*mu_gg_tt(1,:)+squeeze(pi_l_tt(1,9,:))'.*mu_gl_tt(1,:)+squeeze(pi_h_tt(1,9,:))'.*mu_gh_tt(1,:) ).*(1-nu_g_tt(1,:)).*PY_g_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,9,:))'.*mu_lg_tt(1,:)+squeeze(pi_l_tt(1,9,:))'.*mu_ll_tt(1,:)+squeeze(pi_h_tt(1,9,:))'.*mu_lh_tt(1,:) ).*(1-nu_l_tt(1,:)).*PY_l_tt(1,:)+... 
    ( squeeze(pi_g_tt(1,9,:))'.*mu_hg_tt(1,:)+squeeze(pi_l_tt(1,9,:))'.*mu_hl_tt(1,:)+squeeze(pi_h_tt(1,9,:))'.*mu_hh_tt(1,:) ).*(1-nu_h_tt(1,:)).*PY_h_tt(1,:);
im_1_9_intm_tt = im_1_9_intm_tt./IM_1_tt(9,:);
im_1_9_con_tt = ( squeeze(pi_g_tt(1,9,:))'.*omega_cg_tt(1,:)+squeeze(pi_l_tt(1,9,:))'.*omega_cl_tt(1,:)+squeeze(pi_h_tt(1,9,:))'.*omega_ch_tt(1,:) ).*PC_tt(1,:)./IM_1_tt(9,:);
im_1_9_inv_tt = ( squeeze(pi_g_tt(1,9,:))'.*omega_xg_tt(1,:)+squeeze(pi_l_tt(1,9,:))'.*omega_xl_tt(1,:)+squeeze(pi_h_tt(1,9,:))'.*omega_xh_tt(1,:) ).*PX_tt(1,:)./IM_1_tt(9,:);

%% agg bilateral export and import share as in Table 1
figure('Position', get(0, 'Screensize'))
subplot(3,3,1)
plot(1995:1994+ye,ex_1_2_tt(1:ye),'-b',1995:1994+ye,ex_1_3_tt(1:ye),'-r',1995:1994+ye,ex_1_4_tt(1:ye),'-g',1995:1994+ye,ex_1_5_tt(1:ye),'-k',1995:1994+ye,ex_1_6_tt(1:ye),'-c',1995:1994+ye,ex_1_7_tt(1:ye),'-m',1995:1994+ye,ex_1_8_tt(1:ye),'-y',1995:1994+ye,ex_1_9_tt(1:ye),'-o','Linewidth',2.5)
title('AUS export distribution')
legend('BRA','CHN','DEU','IND','JPN','MEX','ROW','USA','location','best')
grid on

subplot(3,3,2)
plot(1995:1994+ye,ex_2_1_tt(1:ye),'-b',1995:1994+ye,ex_2_3_tt(1:ye),'-r',1995:1994+ye,ex_2_4_tt(1:ye),'-g',1995:1994+ye,ex_2_5_tt(1:ye),'-k',1995:1994+ye,ex_2_6_tt(1:ye),'-c',1995:1994+ye,ex_2_7_tt(1:ye),'-m',1995:1994+ye,ex_2_8_tt(1:ye),'-y',1995:1994+ye,ex_2_9_tt(1:ye),'-o','Linewidth',2.5)
title('BRA export distribution')
legend('AUS','CHN','DEU','IND','JPN','MEX','ROW','USA','location','best')
grid on

subplot(3,3,3)
plot(1995:1994+ye,ex_3_1_tt(1:ye),'-b',1995:1994+ye,ex_3_2_tt(1:ye),'-r',1995:1994+ye,ex_3_4_tt(1:ye),'-g',1995:1994+ye,ex_3_5_tt(1:ye),'-k',1995:1994+ye,ex_3_6_tt(1:ye),'-c',1995:1994+ye,ex_3_7_tt(1:ye),'-m',1995:1994+ye,ex_3_8_tt(1:ye),'-y',1995:1994+ye,ex_3_9_tt(1:ye),'-o','Linewidth',2.5)
title('CHN export distribution')
legend('AUS','BRA','DEU','IND','JPN','MEX','ROW','USA','location','best')
grid on

subplot(3,3,4)
plot(1995:1994+ye,ex_4_1_tt(1:ye),'-b',1995:1994+ye,ex_4_2_tt(1:ye),'-r',1995:1994+ye,ex_4_3_tt(1:ye),'-g',1995:1994+ye,ex_4_5_tt(1:ye),'-k',1995:1994+ye,ex_4_6_tt(1:ye),'-c',1995:1994+ye,ex_4_7_tt(1:ye),'-m',1995:1994+ye,ex_4_8_tt(1:ye),'-y',1995:1994+ye,ex_4_9_tt(1:ye),'-o','Linewidth',2.5)
title('DEU export distribution')
legend('AUS','BRA','CHN','IND','JPN','MEX','ROW','USA','location','best')
grid on

subplot(3,3,5)
plot(1995:1994+ye,ex_5_1_tt(1:ye),'-b',1995:1994+ye,ex_5_2_tt(1:ye),'-r',1995:1994+ye,ex_5_3_tt(1:ye),'-g',1995:1994+ye,ex_5_4_tt(1:ye),'-k',1995:1994+ye,ex_5_6_tt(1:ye),'-c',1995:1994+ye,ex_5_7_tt(1:ye),'-m',1995:1994+ye,ex_5_8_tt(1:ye),'-y',1995:1994+ye,ex_5_9_tt(1:ye),'-o','Linewidth',2.5)
title('IND export distribution')
legend('AUS','BRA','CHN','DEU','JPN','MEX','ROW','USA','location','best')
grid on

subplot(3,3,6)
plot(1995:1994+ye,ex_6_1_tt(1:ye),'-b',1995:1994+ye,ex_6_2_tt(1:ye),'-r',1995:1994+ye,ex_6_3_tt(1:ye),'-g',1995:1994+ye,ex_6_4_tt(1:ye),'-k',1995:1994+ye,ex_6_5_tt(1:ye),'-c',1995:1994+ye,ex_6_7_tt(1:ye),'-m',1995:1994+ye,ex_6_8_tt(1:ye),'-y',1995:1994+ye,ex_6_9_tt(1:ye),'-o','Linewidth',2.5)
title('JPN export distribution')
legend('AUS','BRA','CHN','DEU','IND','MEX','ROW','USA','location','best')
grid on

subplot(3,3,7)
plot(1995:1994+ye,ex_7_1_tt(1:ye),'-b',1995:1994+ye,ex_7_2_tt(1:ye),'-r',1995:1994+ye,ex_7_3_tt(1:ye),'-g',1995:1994+ye,ex_7_4_tt(1:ye),'-k',1995:1994+ye,ex_7_5_tt(1:ye),'-c',1995:1994+ye,ex_7_6_tt(1:ye),'-m',1995:1994+ye,ex_7_8_tt(1:ye),'-y',1995:1994+ye,ex_7_9_tt(1:ye),'-o','Linewidth',2.5)
title('MEX export distribution')
legend('AUS','BRA','CHN','DEU','IND','JPN','ROW','USA','location','best')
grid on

subplot(3,3,8)
plot(1995:1994+ye,ex_8_1_tt(1:ye),'-b',1995:1994+ye,ex_8_2_tt(1:ye),'-r',1995:1994+ye,ex_8_3_tt(1:ye),'-g',1995:1994+ye,ex_8_4_tt(1:ye),'-k',1995:1994+ye,ex_8_5_tt(1:ye),'-c',1995:1994+ye,ex_8_6_tt(1:ye),'-m',1995:1994+ye,ex_8_7_tt(1:ye),'-y',1995:1994+ye,ex_8_9_tt(1:ye),'-o','Linewidth',2.5)
title('ROW export distribution')
legend('AUS','BRA','CHN','DEU','IND','JPN','MEX','USA','location','best')
grid on

subplot(3,3,9)
plot(1995:1994+ye,ex_9_1_tt(1:ye),'-b',1995:1994+ye,ex_9_2_tt(1:ye),'-r',1995:1994+ye,ex_9_3_tt(1:ye),'-g',1995:1994+ye,ex_9_4_tt(1:ye),'-k',1995:1994+ye,ex_9_5_tt(1:ye),'-c',1995:1994+ye,ex_9_6_tt(1:ye),'-m',1995:1994+ye,ex_9_7_tt(1:ye),'-y',1995:1994+ye,ex_9_8_tt(1:ye),'-o','Linewidth',2.5)
title('USA export distribution')
legend('AUS','BRA','CHN','DEU','IND','JPN','MEX','USA','location','best')
grid on

sgtitle('Figure. Bilateral export and import share as in Table 1')
saveas(gcf, 'plots_transition/outcome/Exp_imp_share.png');
close(gcf)


% sector export and import share between BRA to paterners as in Table 2 and 3
% figure('Position', get(0, 'Screensize'))
% subplot(2,3,1)
% plot(2000:1999+ye,ex_1_2_a_tt(1:ye),'-b',2000:1999+ye,ex_1_2_m_tt(1:ye),'-r',2000:1999+ye,ex_1_2_h_tt(1:ye),'-g',2000:1999+ye,ex_1_2_s_tt(1:ye),'-m','Linewidth',2.5)
% title('BRA export to CHN by sector')
% legend('P','M','H','S','location','best')
% ylim([0 1])
% grid on
% 
% subplot(2,3,2)
% plot(2000:1999+ye,ex_1_3_a_tt(1:ye),'-b',2000:1999+ye,ex_1_3_m_tt(1:ye),'-r',2000:1999+ye,ex_1_3_h_tt(1:ye),'-g',2000:1999+ye,ex_1_3_s_tt(1:ye),'-m','Linewidth',2.5)
% title('BRA export to USA by sector')
% legend('P','M','H','S','location','best')
% ylim([0 1])
% grid on
% 
% subplot(2,3,3)
% plot(2000:1999+ye,ex_1_4_a_tt(1:ye),'-b',2000:1999+ye,ex_1_4_m_tt(1:ye),'-r',2000:1999+ye,ex_1_4_h_tt(1:ye),'-g',2000:1999+ye,ex_1_4_s_tt(1:ye),'-m','Linewidth',2.5)
% title('BRA export to ROW by sector')
% legend('P','M','H','S','location','best')
% ylim([0 1])
% grid on
% 
% subplot(2,3,4)
% plot(2000:1999+ye,im_1_2_a_tt(1:ye),'-b',2000:1999+ye,im_1_2_m_tt(1:ye),'-r',2000:1999+ye,im_1_2_h_tt(1:ye),'-g',2000:1999+ye,im_1_2_s_tt(1:ye),'-m','Linewidth',2.5)
% title('BRA import from CHN by sector')
% legend('P','M','H','S','location','best')
% ylim([0 1])
% grid on
% 
% subplot(2,3,5)
% plot(2000:1999+ye,im_1_3_a_tt(1:ye),'-b',2000:1999+ye,im_1_3_m_tt(1:ye),'-r',2000:1999+ye,im_1_3_h_tt(1:ye),'-g',2000:1999+ye,im_1_3_s_tt(1:ye),'-m','Linewidth',2.5)
% title('BRA import from USA by sector')
% legend('P','M','H','S','location','best')
% ylim([0 1])
% grid on
% 
% subplot(2,3,6)
% plot(2000:1999+ye,im_1_4_a_tt(1:ye),'-b',2000:1999+ye,im_1_4_m_tt(1:ye),'-r',2000:1999+ye,im_1_4_h_tt(1:ye),'-g',2000:1999+ye,im_1_4_s_tt(1:ye),'-m','Linewidth',2.5)
% title('BRA import from ROW by sector')
% legend('P','M','H','S','location','best')
% ylim([0 1])
% grid on
% 
% sgtitle('Figure. BRA sector export and import share as in Table 2 and 3')
% saveas(gcf, 'plots_transition/BRA_sec_exp_imp_share.png');
% close(gcf)

% export share and import share by use between BRA and paterners
% figure('Position', get(0, 'Screensize'))
% subplot(2,3,1)
% plot(2000:1999+ye,ex_1_2_con_tt(1:ye),'-b',2000:1999+ye,ex_1_2_inv_tt(1:ye),'-r',2000:1999+ye,ex_1_2_intm_tt(1:ye),'-g','Linewidth',2.5)
% title('BRA export to CHN by use')
% legend('Con','Inv','Intm','location','best')
% ylim([0 1])
% grid on
% 
% subplot(2,3,2)
% plot(2000:1999+ye,ex_1_3_con_tt(1:ye),'-b',2000:1999+ye,ex_1_3_inv_tt(1:ye),'-r',2000:1999+ye,ex_1_3_intm_tt(1:ye),'-g','Linewidth',2.5)
% title('BRA export to USA by use')
% legend('Con','Inv','Intm','location','best')
% ylim([0 1])
% grid on
% 
% subplot(2,3,3)
% plot(2000:1999+ye,ex_1_4_con_tt(1:ye),'-b',2000:1999+ye,ex_1_4_inv_tt(1:ye),'-r',2000:1999+ye,ex_1_4_intm_tt(1:ye),'-g','Linewidth',2.5)
% title('BRA export to ROW by use')
% legend('Con','Inv','Intm','location','best')
% ylim([0 1])
% grid on
% 
% subplot(2,3,4)
% plot(2000:1999+ye,im_1_2_con_tt(1:ye),'-b',2000:1999+ye,im_1_2_inv_tt(1:ye),'-r',2000:1999+ye,im_1_2_intm_tt(1:ye),'-g','Linewidth',2.5)
% title('BRA import from CHN by use')
% legend('Con','Inv','Intm','location','best')
% ylim([0 1])
% grid on
% 
% subplot(2,3,5)
% plot(2000:1999+ye,im_1_3_con_tt(1:ye),'-b',2000:1999+ye,im_1_3_inv_tt(1:ye),'-r',2000:1999+ye,im_1_3_intm_tt(1:ye),'-g','Linewidth',2.5)
% title('BRA import from USA by use')
% legend('Con','Inv','Intm','location','best')
% ylim([0 1])
% grid on
% 
% subplot(2,3,6)
% plot(2000:1999+ye,im_1_4_con_tt(1:ye),'-b',2000:1999+ye,im_1_4_inv_tt(1:ye),'-r',2000:1999+ye,im_1_4_intm_tt(1:ye),'-g','Linewidth',2.5)
% title('BRA import from ROW by use')
% legend('Con','Inv','Intm','location','best')
% ylim([0 1])
% grid on
% 
% sgtitle('Figure. BRA export and import share by use')
% saveas(gcf, 'plots_transition/BRA_use_exp_imp_share.png');
% close(gcf)

%% Plots for better illustrations

mkdir plots_transition\final

ye = 2011-1995+1;
% Productivity by sector and country relative to 1995
figure('Position', get(0, 'Screensize'))
subplot(1,3,1)
plot(1995:1994+ye,T_g_tt(1,1:ye)/T_g_tt(1,1),'-b',1995:1994+ye,T_g_tt(2,1:ye)/T_g_tt(2,1),'-r',1995:1994+ye,T_g_tt(3,1:ye)/T_g_tt(3,1),'-g',1995:1994+ye,T_g_tt(4,1:ye)/T_g_tt(4,1),'-c',1995:1994+ye,T_g_tt(5,1:ye)/T_g_tt(5,1),'-k',1995:1994+ye,T_g_tt(7,1:ye)/T_g_tt(7,1),'-m',1995:1994+ye,T_g_tt(9,1:ye)/T_g_tt(9,1),'-y', ...
    'Linewidth',2.5)
title('Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
%ylim([0 3.5])
grid on

subplot(1,3,2)
plot(1995:1994+ye,T_l_tt(1,1:ye)/T_l_tt(1,1),'-b',1995:1994+ye,T_l_tt(2,1:ye)/T_l_tt(2,1),'-r',1995:1994+ye,T_l_tt(3,1:ye)/T_l_tt(3,1),'-g',1995:1994+ye,T_l_tt(4,1:ye)/T_l_tt(4,1),'-c',1995:1994+ye,T_l_tt(5,1:ye)/T_l_tt(5,1),'-k',1995:1994+ye,T_l_tt(7,1:ye)/T_l_tt(7,1),'-m',1995:1994+ye,T_l_tt(9,1:ye)/T_l_tt(9,1),'-y', ...
    'Linewidth',2.5)
title('Low-skill Services')
%ylim([0 3.5])
grid on

subplot(1,3,3)
plot(1995:1994+ye,T_h_tt(1,1:ye)/T_h_tt(1,1),'-b',1995:1994+ye,T_h_tt(2,1:ye)/T_h_tt(2,1),'-r',1995:1994+ye,T_h_tt(3,1:ye)/T_h_tt(3,1),'-g',1995:1994+ye,T_h_tt(4,1:ye)/T_h_tt(4,1),'-c',1995:1994+ye,T_h_tt(5,1:ye)/T_h_tt(5,1),'-k',1995:1994+ye,T_h_tt(7,1:ye)/T_h_tt(7,1),'-m',1995:1994+ye,T_h_tt(9,1:ye)/T_h_tt(9,1),'-y', ...
    'Linewidth',2.5)
title('High-skill Services')
%ylim([0 3.5])
grid on

sgtitle('Figure. Productivity by sector and country relative to 1995')
saveas(gcf, 'plots_transition/final/01_Productivity by sector and country relative to 1995.png');
close(gcf)

% Investment efficiency by country relative to 2000
figure('Position', get(0, 'Screensize'))
subplot(1,2,1)
plot(1995:1994+ye,PX_tt(1,1:ye),'-b',1995:1994+ye,PX_tt(2,1:ye),'-r',1995:1994+ye,PX_tt(3,1:ye),'-g',1995:1994+ye,PX_tt(4,1:ye),'-c',1995:1994+ye,PX_tt(5,1:ye),'-k',1995:1994+ye,PX_tt(7,1:ye),'-m',1995:1994+ye,PX_tt(9,1:ye),'-y', ...
   1995:1995+num_year-1,PX_t,'--black','Linewidth',2.5)
title('Investment spending')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
grid on

subplot(1,2,2)
plot(1995:1994+ye,A_x_tt(1,1:ye)/A_x_tt(1,1),'-b',1995:1994+ye,A_x_tt(2,1:ye)/A_x_tt(2,1),'-r',1995:1994+ye,A_x_tt(3,1:ye)/A_x_tt(3,1),'-g',1995:1994+ye,A_x_tt(4,1:ye)/A_x_tt(4,1),'-c',1995:1994+ye,A_x_tt(5,1:ye)/A_x_tt(5,1),'-k',1995:1994+ye,A_x_tt(7,1:ye)/A_x_tt(7,1),'-m',1995:1994+ye,A_x_tt(9,1:ye)/A_x_tt(9,1),'-y', ...
    'Linewidth',2.5)
title('Investment efficiency')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
grid on

sgtitle('Figure. Investment efficiency by country relative to 1995')
saveas(gcf, 'plots_transition/final/02_Investment efficiency by country relative to 1995.png');
close(gcf)


% % Real exchange rate by country relative to ROW
% figure('Position', get(0, 'Screensize'))
% subplot(1,2,1)
% plot(2000:1999+ye,P_c_tt(1,1:ye)./P_c_tt(4,1:ye),'-b',2000:1999+ye,P_c_tt(2,1:ye)./P_c_tt(4,1:ye),'-r',2000:1999+ye,P_c_tt(3,1:ye)./P_c_tt(4,1:ye),'-g',...
%     'Linewidth',2.5)
% title('Real exchange rate (relative consumption price)')
% legend('BRA','CHN','USA','location','best')
% ylim([0.3 1.3])
% grid on
% 
% subplot(1,2,2)
% plot(2000:1999+ye,P_x_tt(1,1:ye)./P_x_tt(4,1:ye),'-b',2000:1999+ye,P_x_tt(2,1:ye)./P_x_tt(4,1:ye),'-r',2000:1999+ye,P_x_tt(3,1:ye)./P_x_tt(4,1:ye),'-g',...
%     'Linewidth',2.5)
% title('Real exchange rate (relative investment price)')
% ylim([0.3 1.3])
% grid on
% 
% sgtitle('Figure. Real exchange rate by country relative to ROW')
% saveas(gcf, 'plots_transition/final/Real exchange rate by country relative to ROW.png');
% close(gcf)

% % Matching sector ouput 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,PY_a_tt(1,1:ye),'-b',2000:1999+ye,PY_a_tt(2,1:ye),'-r',2000:1999+ye,PY_a_tt(3,1:ye),'-g',2000:1999+ye,PY_a_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,PY_a_t,'--black','Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','ROW','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,PY_m_tt(1,1:ye),'-b',2000:1999+ye,PY_m_tt(2,1:ye),'-r',2000:1999+ye,PY_m_tt(3,1:ye),'-g',2000:1999+ye,PY_m_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,PY_m_t,'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,PY_h_tt(1,1:ye),'-b',2000:1999+ye,PY_h_tt(2,1:ye),'-r',2000:1999+ye,PY_h_tt(3,1:ye),'-g',2000:1999+ye,PY_h_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,PY_h_t,'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,PY_s_tt(1,1:ye),'-b',2000:1999+ye,PY_s_tt(2,1:ye),'-r',2000:1999+ye,PY_s_tt(3,1:ye),'-g',2000:1999+ye,PY_s_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,PY_s_t,'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Matching sector ouput 2000-2014')
% saveas(gcf, 'plots_transition/final/Matching sector ouput 2000-2014.png');
% close(gcf)


% % Matching sector absrption 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,PQ_a_tt(1,1:ye),'-b',2000:1999+ye,PQ_a_tt(2,1:ye),'-r',2000:1999+ye,PQ_a_tt(3,1:ye),'-g',2000:1999+ye,PQ_a_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,PQ_a_t,'--black','Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','ROW','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,PQ_m_tt(1,1:ye),'-b',2000:1999+ye,PQ_m_tt(2,1:ye),'-r',2000:1999+ye,PQ_m_tt(3,1:ye),'-g',2000:1999+ye,PQ_m_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,PQ_m_t,'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,PQ_h_tt(1,1:ye),'-b',2000:1999+ye,PQ_h_tt(2,1:ye),'-r',2000:1999+ye,PQ_h_tt(3,1:ye),'-g',2000:1999+ye,PQ_h_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,PQ_h_t,'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,PQ_s_tt(1,1:ye),'-b',2000:1999+ye,PQ_s_tt(2,1:ye),'-r',2000:1999+ye,PQ_s_tt(3,1:ye),'-g',2000:1999+ye,PQ_s_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,PQ_s_t,'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Matching sector absorption 2000-2014')
% saveas(gcf, 'plots_transition/final/Matching sector absorption 2000-2014.png');
% close(gcf)


% % Matching sector net export 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,NX_a_tt(1,1:ye),'-b',2000:1999+ye,NX_a_tt(2,1:ye),'-r',2000:1999+ye,NX_a_tt(3,1:ye),'-g',2000:1999+ye,NX_a_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,NX_a_t,'--black','Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','ROW','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,NX_m_tt(1,1:ye),'-b',2000:1999+ye,NX_m_tt(2,1:ye),'-r',2000:1999+ye,NX_m_tt(3,1:ye),'-g',2000:1999+ye,NX_m_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,NX_m_t,'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,NX_h_tt(1,1:ye),'-b',2000:1999+ye,NX_h_tt(2,1:ye),'-r',2000:1999+ye,NX_h_tt(3,1:ye),'-g',2000:1999+ye,NX_h_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,NX_h_t,'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,NX_s_tt(1,1:ye),'-b',2000:1999+ye,NX_s_tt(2,1:ye),'-r',2000:1999+ye,NX_s_tt(3,1:ye),'-g',2000:1999+ye,NX_s_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,NX_s_t,'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure.  Matching sector net export 2000-2014')
% saveas(gcf, 'plots_transition/final/Matching sector net export 2000-2014.png');
% close(gcf)


% Fitting bilateral trade in Brazil 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(pi_g_tt(1,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_g_tt(1,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_g_tt(1,4,1:ye)),'-c',2000:1999+ye,squeeze(pi_g_tt(1,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_g_t(1,2:5,:)),'--black','Linewidth',2.5)
% title('Primary')
% legend('CHN','USA','IND','ROW','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(pi_l_tt(1,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_l_tt(1,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_l_tt(1,4,1:ye)),'-c',2000:1999+ye,squeeze(pi_l_tt(1,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_l_t(1,2:5,:)),'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(pi_h_tt(1,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_h_tt(1,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_h_tt(1,4,1:ye)),'-c',2000:1999+ye,squeeze(pi_h_tt(1,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_h_t(1,2:5,:)),'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(pi_s_tt(1,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_s_tt(1,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_s_tt(1,4,1:ye)),'-c',2000:1999+ye,squeeze(pi_s_tt(1,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_s_t(1,2:5,:)),'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Fitting bilateral trade in Brazil 2000-2014')
% saveas(gcf, 'plots_transition/final/Fitting bilateral trade in Brazil 2000-2014.png');
% close(gcf)
% 
% % Fitting bilateral trade in China 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(pi_g_tt(2,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_g_tt(2,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_g_tt(2,4,1:ye)),'-c',2000:1999+ye,squeeze(pi_g_tt(2,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_g_t(2,[1,3,4,5],:)),'--black','Linewidth',2.5)
% title('Primary')
% legend('BRA','USA','IND','ROW','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(pi_l_tt(2,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_l_tt(2,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_l_tt(2,4,1:ye)),'-c',2000:1999+ye,squeeze(pi_l_tt(2,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_l_t(2,[1,3,4,5],:)),'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(pi_h_tt(2,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_h_tt(2,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_h_tt(2,4,1:ye)),'-c',2000:1999+ye,squeeze(pi_h_tt(2,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_h_t(2,[1,3,4,5],:)),'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(pi_s_tt(2,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_s_tt(2,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_s_tt(2,4,1:ye)),'-c',2000:1999+ye,squeeze(pi_s_tt(2,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_s_t(2,[1,3,4,5],:)),'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Fitting bilateral trade in China 2000-2014')
% saveas(gcf, 'plots_transition/final/Fitting bilateral trade in China 2000-2014.png');
% close(gcf)
% 
% 
% % Fitting bilateral trade in USA 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(pi_g_tt(3,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_g_tt(3,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_g_tt(3,4,1:ye)),'-c',2000:1999+ye,squeeze(pi_g_tt(3,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_g_t(3,[1,2,4,5],:)),'--black','Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','IND','ROW','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(pi_l_tt(3,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_l_tt(3,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_l_tt(3,4,1:ye)),'-c',2000:1999+ye,squeeze(pi_l_tt(3,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_l_t(3,[1,2,4,5],:)),'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(pi_h_tt(3,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_h_tt(3,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_h_tt(3,4,1:ye)),'-c',2000:1999+ye,squeeze(pi_h_tt(3,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_h_t(3,[1,2,4,5],:)),'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(pi_s_tt(3,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_s_tt(3,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_s_tt(3,4,1:ye)),'-c',2000:1999+ye,squeeze(pi_s_tt(3,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_s_t(3,[1,2,4,5],:)),'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Fitting bilateral trade in USA 2000-2014')
% saveas(gcf, 'plots_transition/final/Fitting bilateral trade in USA 2000-2014.png');
% close(gcf)
% 
% 
% % Fitting bilateral trade in India 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(pi_g_tt(4,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_g_tt(4,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_g_tt(4,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_g_tt(4,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_g_t(4,[1:3,5],:)),'--black','Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','ROW','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(pi_l_tt(4,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_l_tt(4,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_l_tt(4,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_l_tt(4,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_l_t(4,[1:3,5],:)),'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(pi_h_tt(4,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_h_tt(4,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_h_tt(4,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_h_tt(4,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_h_t(4,[1:3,5],:)),'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(pi_s_tt(4,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_s_tt(4,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_s_tt(4,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_s_tt(4,5,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_s_t(4,[1:3,5],:)),'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Fitting bilateral trade in India 2000-2014')
% saveas(gcf, 'plots_transition/final/Fitting bilateral trade in India 2000-2014.png');
% close(gcf)
% 
% 
% 
% % Fitting bilateral trade in ROW 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(pi_g_tt(5,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_g_tt(5,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_g_tt(5,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_g_tt(5,4,1:ye)),'-c',...
%     2000:2000+num_year-1,squeeze(bt_g_t(5,1:4,:)),'--black','Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','IND','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(pi_l_tt(5,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_l_tt(5,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_l_tt(5,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_l_tt(5,4,1:ye)),'-c',...
%     2000:2000+num_year-1,squeeze(bt_l_t(5,1:4,:)),'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(pi_h_tt(5,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_h_tt(5,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_h_tt(5,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_h_tt(5,4,1:ye)),'-c',...
%     2000:2000+num_year-1,squeeze(bt_h_t(5,1:4,:)),'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(pi_s_tt(5,1,1:ye)),'-b',2000:1999+ye,squeeze(pi_s_tt(5,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_s_tt(5,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_s_tt(5,4,1:ye)),'-c',...
%     2000:2000+num_year-1,squeeze(bt_s_t(5,1:4,:)),'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Fitting bilateral trade in ROW 2000-2014')
% saveas(gcf, 'plots_transition/final/Fitting bilateral trade in ROW 2000-2014.png');
% close(gcf)
% 


% Share of good in consumption 2000-2014
figure('Position', get(0, 'Screensize'))
subplot(2,2,1)
plot(1995:1994+ye,omega_cg_tt(1,1:ye),'-b',1995:1994+ye,omega_cg_tt(2,1:ye),'-r',1995:1994+ye,omega_cg_tt(3,1:ye),'-g',1995:1994+ye,omega_cg_tt(4,1:ye),'-c',1995:1994+ye,omega_cg_tt(5,1:ye),'-k',1995:1994+ye,omega_cg_tt(7,1:ye),'-m',1995:1994+ye,omega_cg_tt(9,1:ye),'-y', ...
     'Linewidth',2.5)
title('Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
grid on

subplot(2,2,2)
plot(1995:1994+ye,omega_cl_tt(1,1:ye),'-b',1995:1994+ye,omega_cl_tt(2,1:ye),'-r',1995:1994+ye,omega_cl_tt(3,1:ye),'-g',1995:1994+ye,omega_cl_tt(4,1:ye),'-c',1995:1994+ye,omega_cl_tt(5,1:ye),'-k',1995:1994+ye,omega_cl_tt(7,1:ye),'-m',1995:1994+ye,omega_cl_tt(9,1:ye),'-y', ...
     'Linewidth',2.5)
title('Low-skill services')
grid on

subplot(2,2,3)
plot(1995:1994+ye,omega_ch_tt(1,1:ye),'-b',1995:1994+ye,omega_ch_tt(2,1:ye),'-r',1995:1994+ye,omega_ch_tt(3,1:ye),'-g',1995:1994+ye,omega_ch_tt(4,1:ye),'-c',1995:1994+ye,omega_ch_tt(5,1:ye),'-k',1995:1994+ye,omega_ch_tt(7,1:ye),'-m',1995:1994+ye,omega_ch_tt(9,1:ye),'-y', ...
     'Linewidth',2.5)
title('High-skill services')
grid on

sgtitle('Figure 2. Shares in consumption 1995-2011')
%sgtitle('Shares of goods in consumption 2000-2014')
saveas(gcf, 'plots_transition/final_2/Figure 2.png');
close(gcf)

% Share of good in investment 2000-2014
figure('Position', get(0, 'Screensize'))
subplot(2,2,1)
plot(1995:1994+ye,omega_xg_tt(1,1:ye),'-b',1995:1994+ye,omega_xg_tt(2,1:ye),'-r',1995:1994+ye,omega_xg_tt(3,1:ye),'-g',1995:1994+ye,omega_xg_tt(4,1:ye),'-c',1995:1994+ye,omega_xg_tt(5,1:ye),'-k',1995:1994+ye,omega_xg_tt(7,1:ye),'-m',1995:1994+ye,omega_xg_tt(9,1:ye),'-y', ...
     'Linewidth',2.5)
title('Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
grid on

subplot(2,2,2)
plot(1995:1994+ye,omega_xl_tt(1,1:ye),'-b',1995:1994+ye,omega_xl_tt(2,1:ye),'-r',1995:1994+ye,omega_xl_tt(3,1:ye),'-g',1995:1994+ye,omega_xl_tt(4,1:ye),'-c',1995:1994+ye,omega_xl_tt(5,1:ye),'-k',1995:1994+ye,omega_xl_tt(7,1:ye),'-m',1995:1994+ye,omega_xl_tt(9,1:ye),'-y', ...
     'Linewidth',2.5)
title('Low-skill services')
grid on

subplot(2,2,3)
plot(1995:1994+ye,omega_xh_tt(1,1:ye),'-b',1995:1994+ye,omega_xh_tt(2,1:ye),'-r',1995:1994+ye,omega_xh_tt(3,1:ye),'-g',1995:1994+ye,omega_xh_tt(4,1:ye),'-c',1995:1994+ye,omega_xh_tt(5,1:ye),'-k',1995:1994+ye,omega_xh_tt(7,1:ye),'-m',1995:1994+ye,omega_xh_tt(9,1:ye),'-y', ...
     'Linewidth',2.5)
title('High-skill services')
grid on

sgtitle('Figure 1. Shares in investment composite good 1995-2011')
%sgtitle('Shares of goods in investment 2000-2014')
saveas(gcf, 'plots_transition/final_2/Figure 1.png');
close(gcf)


% % Bilateral cost of shipping goods to Brazil 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(d_a_tt(1,2,1:ye)),'-r',2000:1999+ye,squeeze(d_a_tt(1,3,1:ye)),'-g',2000:1999+ye,squeeze(d_a_tt(1,4,1:ye)),'-m',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('Primary')
% legend('From CHN','From USA','From ROW','Value of 1','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(d_m_tt(1,2,1:ye)),'-r',2000:1999+ye,squeeze(d_m_tt(1,3,1:ye)),'-g',2000:1999+ye,squeeze(d_m_tt(1,4,1:ye)),'-m',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(d_h_tt(1,2,1:ye)),'-r',2000:1999+ye,squeeze(d_h_tt(1,3,1:ye)),'-g',2000:1999+ye,squeeze(d_h_tt(1,4,1:ye)),'-m',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(d_s_tt(1,2,1:ye)),'-r',2000:1999+ye,squeeze(d_s_tt(1,3,1:ye)),'-g',2000:1999+ye,squeeze(d_s_tt(1,4,1:ye)),'-m',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Bilateral cost of shipping goods to Brazil 2000-2014')
% saveas(gcf, 'plots_transition/final/Bilateral cost of shipping goods to Brazil 2000-2014.png');
% close(gcf)

% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(d_a_tt(1,2,1:ye))./d_a_tt(1,2,1),'-r',2000:1999+ye,squeeze(d_a_tt(1,3,1:ye))./d_a_tt(1,3,1),'-g',2000:1999+ye,squeeze(d_a_tt(1,4,1:ye))./d_a_tt(1,4,1),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('From CHN','From USA','From ROW','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(d_m_tt(1,2,1:ye))./d_m_tt(1,2,1),'-r',2000:1999+ye,squeeze(d_m_tt(1,3,1:ye))./d_m_tt(1,3,1),'-g',2000:1999+ye,squeeze(d_m_tt(1,4,1:ye))./d_m_tt(1,4,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(d_h_tt(1,2,1:ye))./d_h_tt(1,2,1),'-r',2000:1999+ye,squeeze(d_h_tt(1,3,1:ye))./d_h_tt(1,3,1),'-g',2000:1999+ye,squeeze(d_h_tt(1,4,1:ye))./d_h_tt(1,4,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(d_s_tt(1,2,1:ye))./d_s_tt(1,2,1),'-r',2000:1999+ye,squeeze(d_s_tt(1,3,1:ye))./d_s_tt(1,3,1),'-g',2000:1999+ye,squeeze(d_s_tt(1,4,1:ye))./d_s_tt(1,4,1),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Bilateral cost of shipping goods to Brazil relative to 2000')
% saveas(gcf, 'plots_transition/final/Bilateral cost of shipping goods to Brazil relative to 2000.png');
% close(gcf)


% % Bilateral cost of shipping goods to China 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(d_a_tt(2,1,1:ye)),'-b',2000:1999+ye,squeeze(d_a_tt(2,3,1:ye)),'-g',2000:1999+ye,squeeze(d_a_tt(2,4,1:ye)),'-m',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('Primary')
% legend('From BRA','From USA','From ROW','Value of 1','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(d_m_tt(2,1,1:ye)),'-b',2000:1999+ye,squeeze(d_m_tt(2,3,1:ye)),'-g',2000:1999+ye,squeeze(d_m_tt(2,4,1:ye)),'-m',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(d_h_tt(2,1,1:ye)),'-b',2000:1999+ye,squeeze(d_h_tt(2,3,1:ye)),'-g',2000:1999+ye,squeeze(d_h_tt(2,4,1:ye)),'-m',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(d_s_tt(2,1,1:ye)),'-b',2000:1999+ye,squeeze(d_s_tt(2,3,1:ye)),'-g',2000:1999+ye,squeeze(d_s_tt(2,4,1:ye)),'-m',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Bilateral cost of shipping goods to China 2000-2014')
% saveas(gcf, 'plots_transition/final/Bilateral cost of shipping goods to China 2000-2014.png');
% close(gcf)

% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(d_a_tt(2,1,1:ye))./d_a_tt(2,1,1),'-b',2000:1999+ye,squeeze(d_a_tt(2,3,1:ye))./d_a_tt(2,3,1),'-g',2000:1999+ye,squeeze(d_a_tt(2,4,1:ye))./d_a_tt(2,4,1),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('From BRA','From USA','From ROW','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(d_m_tt(2,1,1:ye))./d_m_tt(2,1,1),'-b',2000:1999+ye,squeeze(d_m_tt(2,3,1:ye))./d_m_tt(2,3,1),'-g',2000:1999+ye,squeeze(d_m_tt(2,4,1:ye))./d_m_tt(2,4,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(d_h_tt(2,1,1:ye))./d_h_tt(2,1,1),'-b',2000:1999+ye,squeeze(d_h_tt(2,3,1:ye))./d_h_tt(2,3,1),'-g',2000:1999+ye,squeeze(d_h_tt(2,4,1:ye))./d_h_tt(2,4,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(d_s_tt(2,1,1:ye))./d_s_tt(2,1,1),'-b',2000:1999+ye,squeeze(d_s_tt(2,3,1:ye))./d_s_tt(2,3,1),'-g',2000:1999+ye,squeeze(d_s_tt(2,4,1:ye))./d_s_tt(2,4,1),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Bilateral cost of shipping goods to China relative to 2000')
% saveas(gcf, 'plots_transition/final/Bilateral cost of shipping goods to China relative to 2000.png');
% close(gcf)


% % Bilateral cost of shipping goods to USA 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(d_a_tt(3,1,1:ye)),'-b',2000:1999+ye,squeeze(d_a_tt(3,2,1:ye)),'-r',2000:1999+ye,squeeze(d_a_tt(3,4,1:ye)),'-m',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('Primary')
% legend('From BRA','From CHN','From ROW','Value of 1','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(d_m_tt(3,1,1:ye)),'-b',2000:1999+ye,squeeze(d_m_tt(3,2,1:ye)),'-r',2000:1999+ye,squeeze(d_m_tt(3,4,1:ye)),'-m',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(d_h_tt(3,1,1:ye)),'-b',2000:1999+ye,squeeze(d_h_tt(3,2,1:ye)),'-r',2000:1999+ye,squeeze(d_h_tt(3,4,1:ye)),'-m',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(d_s_tt(3,1,1:ye)),'-b',2000:1999+ye,squeeze(d_s_tt(3,2,1:ye)),'-r',2000:1999+ye,squeeze(d_s_tt(3,4,1:ye)),'-m',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Bilateral cost of shipping goods to USA 2000-2014')
% saveas(gcf, 'plots_transition/final/Bilateral cost of shipping goods to USA 2000-2014.png');
% close(gcf)

% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(d_a_tt(3,1,1:ye))./d_a_tt(3,1,1),'-b',2000:1999+ye,squeeze(d_a_tt(3,2,1:ye))./d_a_tt(3,2,1),'-r',2000:1999+ye,squeeze(d_a_tt(3,4,1:ye))./d_a_tt(3,4,1),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('From BRA','From CHN','From ROW','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(d_m_tt(3,1,1:ye))./d_m_tt(3,1,1),'-b',2000:1999+ye,squeeze(d_m_tt(3,2,1:ye))./d_m_tt(3,2,1),'-r',2000:1999+ye,squeeze(d_m_tt(3,4,1:ye))./d_m_tt(3,4,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(d_h_tt(3,1,1:ye))./d_h_tt(3,1,1),'-b',2000:1999+ye,squeeze(d_h_tt(3,2,1:ye))./d_h_tt(3,2,1),'-r',2000:1999+ye,squeeze(d_h_tt(3,4,1:ye))./d_h_tt(3,4,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(d_s_tt(3,1,1:ye))./d_s_tt(3,1,1),'-b',2000:1999+ye,squeeze(d_s_tt(3,2,1:ye))./d_s_tt(3,2,1),'-r',2000:1999+ye,squeeze(d_s_tt(3,4,1:ye))./d_s_tt(3,4,1),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Bilateral cost of shipping goods to USA relative to 2000')
% saveas(gcf, 'plots_transition/final/Bilateral cost of shipping goods to USA relative to 2000.png');
% close(gcf)

% % Bilateral cost of shipping goods to ROW 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(d_a_tt(4,1,1:ye)),'-b',2000:1999+ye,squeeze(d_a_tt(4,2,1:ye)),'-r',2000:1999+ye,squeeze(d_a_tt(4,3,1:ye)),'-g',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('Primary')
% legend('From BRA','From CHN','From USA','Value of 1','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(d_m_tt(4,1,1:ye)),'-b',2000:1999+ye,squeeze(d_m_tt(4,2,1:ye)),'-r',2000:1999+ye,squeeze(d_m_tt(4,3,1:ye)),'-g',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(d_h_tt(4,1,1:ye)),'-b',2000:1999+ye,squeeze(d_h_tt(4,2,1:ye)),'-r',2000:1999+ye,squeeze(d_h_tt(4,3,1:ye)),'-g',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(d_s_tt(4,1,1:ye)),'-b',2000:1999+ye,squeeze(d_s_tt(4,2,1:ye)),'-r',2000:1999+ye,squeeze(d_s_tt(4,3,1:ye)),'-g',...
%     2000:1999+ye,ones(I,ye),'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Bilateral cost of shipping goods to ROW 2000-2014')
% saveas(gcf, 'plots_transition/final/Bilateral cost of shipping goods to ROW 2000-2014.png');
% close(gcf)

% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(d_a_tt(4,1,1:ye))./d_a_tt(4,1,1),'-b',2000:1999+ye,squeeze(d_a_tt(4,2,1:ye))./d_a_tt(4,2,1),'-r',2000:1999+ye,squeeze(d_a_tt(4,3,1:ye))./d_a_tt(4,3,1),'-g',...
%     'Linewidth',2.5)
% title('Primary')
% legend('From BRA','From CHN','From USA','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(d_m_tt(4,1,1:ye))./d_m_tt(4,1,1),'-b',2000:1999+ye,squeeze(d_m_tt(4,2,1:ye))./d_m_tt(4,2,1),'-r',2000:1999+ye,squeeze(d_m_tt(4,3,1:ye))./d_m_tt(4,3,1),'-g',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(d_h_tt(4,1,1:ye))./d_h_tt(4,1,1),'-b',2000:1999+ye,squeeze(d_h_tt(4,2,1:ye))./d_h_tt(4,2,1),'-r',2000:1999+ye,squeeze(d_h_tt(4,3,1:ye))./d_h_tt(4,3,1),'-g',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(d_s_tt(4,1,1:ye))./d_s_tt(4,1,1),'-b',2000:1999+ye,squeeze(d_s_tt(4,2,1:ye))./d_s_tt(4,2,1),'-r',2000:1999+ye,squeeze(d_s_tt(4,3,1:ye))./d_s_tt(4,3,1),'-g',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Bilateral cost of shipping goods to ROW relative to 2000')
% saveas(gcf, 'plots_transition/final/Bilateral cost of shipping goods to ROW relative to 2000.png');
% close(gcf)

% % Static global portfolio 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(1,2,1)
% plot(2000:1999+num_year,NX_portfolio_t(1,:),'-b',2000:1999+num_year,NX_portfolio_t(2,:),'-r',2000:1999+num_year,NX_portfolio_t(3,:),'-g',2000:1999+num_year,NX_portfolio_t(4,:),'-m',...
%     2000:1999+num_year,NX_t(1,:),'--b',2000:1999+num_year,NX_t(2,:),'--r',2000:1999+num_year,NX_t(3,:),'--g',2000:1999+num_year,NX_t(4,:),'--m','Linewidth',2.5)
% title('Net export')
% legend('BRA', 'CHN', 'USA', 'ROW','Raw data','location','northeast') 
% grid on
% 
% subplot(1,2,2)
% plot(2000:1999+num_year,fi_t(1,:),'-b',2000:1999+num_year,fi_t(2,:),'-r',2000:1999+num_year,fi_t(3,:),'-g',2000:1999+num_year,fi_t(4,:),'-m','Linewidth',2.5)
% title('Static global portfolio')
% legend('BRA', 'CHN', 'USA', 'ROW','Raw data','location','northeast') 
% grid on
% 
% sgtitle('Figure. Static global portfolio')
% saveas(gcf, 'plots_transition/final/Static global portfolio 2000-2014.png');
% close(gcf)


% % Fitting sector price relative to ROW 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,P_a_tt(1,1:ye)./P_a_tt(4,1:ye),'-b',2000:1999+ye,P_a_tt(2,1:ye)./P_a_tt(4,1:ye),'-r',2000:1999+ye,P_a_tt(3,1:ye)./P_a_tt(4,1:ye),'-g',...
%     2000:2000+num_year-1,P_a_t([1,2,3],1:ye)./P_a_t(4,1:ye),'--black','Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,P_m_tt(1,1:ye)./P_m_tt(4,1:ye),'-b',2000:1999+ye,P_m_tt(2,1:ye)./P_m_tt(4,1:ye),'-r',2000:1999+ye,P_m_tt(3,1:ye)./P_m_tt(4,1:ye),'-g',...
%     2000:2000+num_year-1,P_m_t([1,2,3],1:ye)./P_m_t(4,1:ye),'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,P_h_tt(1,1:ye)./P_h_tt(4,1:ye),'-b',2000:1999+ye,P_h_tt(2,1:ye)./P_h_tt(4,1:ye),'-r',2000:1999+ye,P_h_tt(3,1:ye)./P_h_tt(4,1:ye),'-g',...
%     2000:2000+num_year-1,P_h_t([1,2,3],1:ye)./P_h_t(4,1:ye),'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,P_s_tt(1,1:ye)./P_s_tt(4,1:ye),'-b',2000:1999+ye,P_s_tt(2,1:ye)./P_s_tt(4,1:ye),'-r',2000:1999+ye,P_s_tt(3,1:ye)./P_s_tt(4,1:ye),'-g',...
%     2000:2000+num_year-1,P_s_t([1,2,3],1:ye)./P_s_t(4,1:ye),'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Fitting sector price relative to ROW 2000-2014')
% saveas(gcf, 'plots_transition/final/Fitting sector price relative to ROW 2000-2014.png');
% close(gcf)



%% Final 2

ye = 2011-1995+1;
%Productivity by sector and country relative to 1995
figure('Position', get(0, 'Screensize'))
subplot(1,3,1)
plot(1995:1994+ye,T_g_tt(1,1:ye)/T_g_tt(1,1),'-b',1995:1994+ye,T_g_tt(2,1:ye)/T_g_tt(2,1),'-r',1995:1994+ye,T_g_tt(3,1:ye)/T_g_tt(3,1),'-g',1995:1994+ye,T_g_tt(4,1:ye)/T_g_tt(4,1),'-c',1995:1994+ye,T_g_tt(5,1:ye)/T_g_tt(5,1),'-k',1995:1994+ye,T_g_tt(7,1:ye)/T_g_tt(7,1),'-m',1995:1994+ye,T_g_tt(9,1:ye)/T_g_tt(9,1),'-y', ...
     'Linewidth',2.5)
title('Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
%ylim([0 3])
grid on

subplot(1,3,2)
plot(1995:1994+ye,T_l_tt(1,1:ye)/T_l_tt(1,1),'-b',1995:1994+ye,T_l_tt(2,1:ye)/T_l_tt(2,1),'-r',1995:1994+ye,T_l_tt(3,1:ye)/T_l_tt(3,1),'-g',1995:1994+ye,T_l_tt(4,1:ye)/T_l_tt(4,1),'-c',1995:1994+ye,T_l_tt(5,1:ye)/T_l_tt(5,1),'-k',1995:1994+ye,T_l_tt(7,1:ye)/T_l_tt(7,1),'-m',1995:1994+ye,T_l_tt(9,1:ye)/T_l_tt(9,1),'-y', ...
     'Linewidth',2.5)
title('Low-skill services')
%ylim([0 3])
grid on

subplot(1,3,3)
plot(1995:1994+ye,T_h_tt(1,1:ye)/T_h_tt(1,1),'-b',1995:1994+ye,T_h_tt(2,1:ye)/T_h_tt(2,1),'-r',1995:1994+ye,T_h_tt(3,1:ye)/T_h_tt(3,1),'-g',1995:1994+ye,T_h_tt(4,1:ye)/T_h_tt(4,1),'-c',1995:1994+ye,T_h_tt(5,1:ye)/T_h_tt(5,1),'-k',1995:1994+ye,T_h_tt(7,1:ye)/T_h_tt(7,1),'-m',1995:1994+ye,T_h_tt(9,1:ye)/T_h_tt(9,1),'-y', ...
     'Linewidth',2.5)
title('High-skill services')
%ylim([0 3])
grid on

sgtitle('Figure 8. Productivity by sector and country relative to 1995')
%sgtitle('Productivity by sector and country relative to 2000')
saveas(gcf, 'plots_transition/final_2/Figure 8.png');
close(gcf)

% Productivity by sector and country relative to USA
figure('Position', get(0, 'Screensize'))
subplot(1,3,1)
plot(1995:1994+ye,T_g_tt(1,1:ye)./T_g_tt(9,1:ye),'-b',1995:1994+ye,T_g_tt(2,1:ye)./T_g_tt(9,1:ye),'-r',1995:1994+ye,T_g_tt(3,1:ye)./T_g_tt(9,1:ye),'-g',1995:1994+ye,T_g_tt(4,1:ye)./T_g_tt(9,1:ye),'-c',1995:1994+ye,T_g_tt(5,1:ye)./T_g_tt(9,1:ye),'-k',1995:1994+ye,T_g_tt(7,1:ye)./T_g_tt(9,1:ye),'-m', ...
     'Linewidth',2.5)
title('Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','location','best') 
%ylim([0 3])
grid on

subplot(1,3,2)
plot(1995:1994+ye,T_l_tt(1,1:ye)./T_l_tt(9,1:ye),'-b',1995:1994+ye,T_l_tt(2,1:ye)./T_l_tt(9,1:ye),'-r',1995:1994+ye,T_l_tt(3,1:ye)./T_l_tt(9,1:ye),'-g',1995:1994+ye,T_l_tt(4,1:ye)./T_l_tt(9,1:ye),'-c',1995:1994+ye,T_l_tt(5,1:ye)./T_l_tt(9,1:ye),'-k',1995:1994+ye,T_l_tt(7,1:ye)./T_l_tt(9,1:ye),'-m', ...
     'Linewidth',2.5)
title('Low-skill services')
%ylim([0 3])
grid on

subplot(1,3,3)
plot(1995:1994+ye,T_h_tt(1,1:ye)./T_h_tt(9,1:ye),'-b',1995:1994+ye,T_h_tt(2,1:ye)./T_h_tt(9,1:ye),'-r',1995:1994+ye,T_h_tt(3,1:ye)./T_h_tt(9,1:ye),'-g',1995:1994+ye,T_h_tt(4,1:ye)./T_h_tt(9,1:ye),'-c',1995:1994+ye,T_h_tt(5,1:ye)./T_h_tt(9,1:ye),'-k',1995:1994+ye,T_h_tt(7,1:ye)./T_h_tt(9,1:ye),'-m', ...
     'Linewidth',2.5)
title('High-skill services')
%ylim([0 3])
grid on

%sgtitle('Figure 5. Productivity by sector and country relative to 2000')
sgtitle('Productivity relative to USA')
saveas(gcf, 'plots_transition/final_2/Figure 5b.png');
close(gcf)

% Unit cost relative to USA
figure('Position', get(0, 'Screensize'))
subplot(1,3,1)
plot(1995:1994+ye,u_g_tt(1,1:ye)./u_g_tt(9,1:ye),'-b',1995:1994+ye,u_g_tt(2,1:ye)./u_g_tt(9,1:ye),'-r',1995:1994+ye,u_g_tt(3,1:ye)./u_g_tt(9,1:ye),'-g',1995:1994+ye,u_g_tt(4,1:ye)./u_g_tt(9,1:ye),'-c',1995:1994+ye,u_g_tt(5,1:ye)./u_g_tt(9,1:ye),'-k',1995:1994+ye,u_g_tt(7,1:ye)./u_g_tt(9,1:ye),'-m', ...
     'Linewidth',2.5)
title('Goods')
legend('AUS','BRA','CHN','DEU','IND','MEX','location','best') 
%ylim([0 3])
grid on

subplot(1,3,2)
plot(1995:1994+ye,u_l_tt(1,1:ye)./u_l_tt(9,1:ye),'-b',1995:1994+ye,u_l_tt(2,1:ye)./u_l_tt(9,1:ye),'-r',1995:1994+ye,u_l_tt(3,1:ye)./u_l_tt(9,1:ye),'-g',1995:1994+ye,u_l_tt(4,1:ye)./u_l_tt(9,1:ye),'-c',1995:1994+ye,u_l_tt(5,1:ye)./u_l_tt(9,1:ye),'-k',1995:1994+ye,u_l_tt(7,1:ye)./u_l_tt(9,1:ye),'-m', ...
     'Linewidth',2.5)
title('Low-skill services')
%ylim([0 3])
grid on

subplot(1,3,3)
plot(1995:1994+ye,u_h_tt(1,1:ye)./u_h_tt(9,1:ye),'-b',1995:1994+ye,u_h_tt(2,1:ye)./u_h_tt(9,1:ye),'-r',1995:1994+ye,u_h_tt(3,1:ye)./u_h_tt(9,1:ye),'-g',1995:1994+ye,u_h_tt(4,1:ye)./u_h_tt(9,1:ye),'-c',1995:1994+ye,u_h_tt(5,1:ye)./u_h_tt(9,1:ye),'-k',1995:1994+ye,u_h_tt(7,1:ye)./u_h_tt(9,1:ye),'-m', ...
     'Linewidth',2.5)
title('High-skill services')
%ylim([0 3])
grid on

%sgtitle('Figure 5. Productivity by sector and country relative to 2000')
sgtitle('Unit cost relative to USA')
saveas(gcf, 'plots_transition/final_2/Figure 5c.png');
close(gcf)

% Investment efficiency by country relative to 2000
figure('Position', get(0, 'Screensize'))
subplot(1,2,1)
plot(1995:1994+ye,PX_tt(1,1:ye),'-b',1995:1994+ye,PX_tt(2,1:ye),'-r',1995:1994+ye,PX_tt(3,1:ye),'-g',1995:1994+ye,PX_tt(4,1:ye),'-c',1995:1994+ye,PX_tt(5,1:ye),'-k',1995:1994+ye,PX_tt(7,1:ye),'-m',1995:1994+ye,PX_tt(9,1:ye),'-y', ...
    1995:1994+ye,PX_t(1,1:ye),'--b',1995:1994+ye,PX_t(2,1:ye),'--r',1995:1994+ye,PX_t(3,1:ye),'--g',1995:1994+ye,PX_t(4,1:ye),'--c',1995:1994+ye,PX_t(5,1:ye),'--k',1995:1994+ye,PX_t(7,1:ye),'--m',1995:1994+ye,PX_t(9,1:ye),'--y', ...
    'Linewidth',2.5)
title('Investment spending over world GDP')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
grid on

subplot(1,2,2)
plot(1995:1994+ye,A_x_tt(1,1:ye)/A_x_tt(1,1),'-b',1995:1994+ye,A_x_tt(2,1:ye)/A_x_tt(2,1),'-r',1995:1994+ye,A_x_tt(3,1:ye)/A_x_tt(3,1),'-g',1995:1994+ye,A_x_tt(4,1:ye)/A_x_tt(4,1),'-c',1995:1994+ye,A_x_tt(5,1:ye)/A_x_tt(5,1),'-k',1995:1994+ye,A_x_tt(7,1:ye)/A_x_tt(7,1),'-m',1995:1994+ye,A_x_tt(9,1:ye)/A_x_tt(9,1),'-y', ...
     'Linewidth',2.5)
title('Investment efficiency')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
grid on

%sgtitle('Figure 5. Investment efficiency by country relative to 2000')
sgtitle('Investment efficiency by country relative to 2000')
saveas(gcf, 'plots_transition/final_2/Figure 5d.png');
close(gcf)

% Investment efficiency by country relative to 2000
figure('Position', get(0, 'Screensize'))
subplot(1,2,1)
plot(1995:1994+ye,PX_tt(1,1:ye)./GDP_tt(1,1:ye),'-b',1995:1994+ye,PX_tt(2,1:ye)./GDP_tt(2,1:ye),'-r',1995:1994+ye,PX_tt(3,1:ye)./GDP_tt(3,1:ye),'-g',1995:1994+ye,PX_tt(4,1:ye)./GDP_tt(4,1:ye),'-c',1995:1994+ye,PX_tt(5,1:ye)./GDP_tt(5,1:ye),'-k',1995:1994+ye,PX_tt(7,1:ye)./GDP_tt(7,1:ye),'-m',1995:1994+ye,PX_tt(9,1:ye)./GDP_tt(9,1:ye),'-y', ...
    1995:1994+ye,PX_t(1,1:ye)./GDP_t(1,1:ye),'--b',1995:1994+ye,PX_t(2,1:ye)./GDP_t(2,1:ye),'--r',1995:1994+ye,PX_t(3,1:ye)./GDP_t(3,1:ye),'--g',1995:1994+ye,PX_t(4,1:ye)./GDP_t(4,1:ye),'--c',1995:1994+ye,PX_t(5,1:ye)./GDP_t(5,1:ye),'--k',1995:1994+ye,PX_t(7,1:ye)./GDP_t(7,1:ye),'--m',1995:1994+ye,PX_t(9,1:ye)./GDP_t(9,1:ye),'--y', ...
    'Linewidth',2.5)
title('Investment spending over GDP')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','Target','location','best') 
grid on

subplot(1,2,2)
plot(1995:1994+ye,A_x_tt(1,1:ye)/A_x_tt(1,1),'-b',1995:1994+ye,A_x_tt(2,1:ye)/A_x_tt(2,1),'-r',1995:1994+ye,A_x_tt(3,1:ye)/A_x_tt(3,1),'-g',1995:1994+ye,A_x_tt(4,1:ye)/A_x_tt(4,1),'-c',1995:1994+ye,A_x_tt(5,1:ye)/A_x_tt(5,1),'-k',1995:1994+ye,A_x_tt(7,1:ye)/A_x_tt(7,1),'-m',1995:1994+ye,A_x_tt(9,1:ye)/A_x_tt(9,1),'-y', ...
     'Linewidth',2.5)
title('Investment efficiency')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','best') 
grid on

sgtitle('Figure 6. Investment efficiency by country relative to 1995')
%sgtitle('Investment efficiency by country relative to 2000')
saveas(gcf, 'plots_transition/final_2/Figure 6.png');
close(gcf)

%India
figure('Position', get(0, 'Screensize'))
subplot(1,2,1)
plot(1995:1994+ye,PX_tt(3,1:ye),'-b',...
    1995:1995+num_year-1,PX_t(3,1:ye),'--black','Linewidth',2.5)
title('Investment spending over WGDP')
legend('CHN','Target','location','best')
ylim([0 1e-2])
grid on

subplot(1,2,2)
plot(1995:1994+ye,PX_tt(5,1:ye),'-b',...
    1995:1995+num_year-1,PX_t(5,1:ye),'--black','Linewidth',2.5)
title('Investment spending over WGDP')
legend('IND','Target','location','best')
ylim([0 1e-2])
grid on

sgtitle('Figure. Investment efficiency India relative to 1995')
saveas(gcf, 'plots_transition/Investment efficiency India and China relative to 1995.png');
close(gcf)



% Bilateral cost of shipping goods from India

figure('Position', get(0, 'Screensize'))
subplot(2,3,1)
plot(1995:1994+ye,squeeze(d_g_tt(1,5,1:ye))./d_g_tt(1,5,1),'-b',1995:1994+ye,squeeze(d_g_tt(2,5,1:ye))./d_g_tt(2,5,1),'-r',1995:1994+ye,squeeze(d_g_tt(3,5,1:ye))./d_g_tt(3,5,1),'-g',1995:1994+ye,squeeze(d_g_tt(4,5,1:ye))./d_g_tt(4,5,1),'-c',1995:1994+ye,squeeze(d_g_tt(7,5,1:ye))./d_g_tt(7,5,1),'-m',1995:1994+ye,squeeze(d_g_tt(9,5,1:ye))./d_g_tt(9,5,1),'-y', ...
     'Linewidth',2.5)
title('Goods Exports')
legend('To AUS','To BRA','To CHN','To DEU','To MEX','To USA','location','best') 
grid on

subplot(2,3,2)
plot(1995:1994+ye,squeeze(d_l_tt(1,5,1:ye))./d_l_tt(1,5,1),'-b',1995:1994+ye,squeeze(d_l_tt(2,5,1:ye))./d_l_tt(2,5,1),'-r',1995:1994+ye,squeeze(d_l_tt(3,5,1:ye))./d_l_tt(3,5,1),'-g',1995:1994+ye,squeeze(d_l_tt(4,5,1:ye))./d_l_tt(4,5,1),'-c',1995:1994+ye,squeeze(d_l_tt(7,5,1:ye))./d_l_tt(7,5,1),'-m',1995:1994+ye,squeeze(d_l_tt(9,5,1:ye))./d_l_tt(9,5,1),'-y', ...
     'Linewidth',2.5)
title('Low-skill service exports')
grid on

subplot(2,3,3)
plot(1995:1994+ye,squeeze(d_h_tt(1,5,1:ye))./d_h_tt(1,5,1),'-b',1995:1994+ye,squeeze(d_h_tt(2,5,1:ye))./d_h_tt(2,5,1),'-r',1995:1994+ye,squeeze(d_h_tt(3,5,1:ye))./d_h_tt(3,5,1),'-g',1995:1994+ye,squeeze(d_h_tt(4,5,1:ye))./d_h_tt(4,5,1),'-c',1995:1994+ye,squeeze(d_h_tt(7,5,1:ye))./d_h_tt(7,5,1),'-m',1995:1994+ye,squeeze(d_h_tt(9,5,1:ye))./d_h_tt(9,5,1),'-y', ...
     'Linewidth',2.5)
title('High-skill service exports')
grid on

subplot(2,3,4)
plot(1995:1994+ye,squeeze(d_g_tt(5,1,1:ye))./d_g_tt(5,1,1),'-b',1995:1994+ye,squeeze(d_g_tt(5,2,1:ye))./d_g_tt(5,2,1),'-r',1995:1994+ye,squeeze(d_g_tt(5,3,1:ye))./d_g_tt(5,3,1),'-g',1995:1994+ye,squeeze(d_g_tt(5,4,1:ye))./d_g_tt(5,4,1),'-c',1995:1994+ye,squeeze(d_g_tt(5,7,1:ye))./d_g_tt(5,7,1),'-m',1995:1994+ye,squeeze(d_g_tt(5,9,1:ye))./d_g_tt(5,9,1),'-y', ...
     'Linewidth',2.5)
title('Goods Imports')
legend('From AUS','From BRA','From CHN','From DEU','From MEX','From USA','location','best') 
grid on

subplot(2,3,5)
plot(1995:1994+ye,squeeze(d_l_tt(5,1,1:ye))./d_l_tt(5,1,1),'-b',1995:1994+ye,squeeze(d_l_tt(5,2,1:ye))./d_l_tt(5,2,1),'-r',1995:1994+ye,squeeze(d_l_tt(5,3,1:ye))./d_l_tt(5,3,1),'-g',1995:1994+ye,squeeze(d_l_tt(5,4,1:ye))./d_l_tt(5,4,1),'-c',1995:1994+ye,squeeze(d_l_tt(5,7,1:ye))./d_l_tt(5,7,1),'-m',1995:1994+ye,squeeze(d_l_tt(5,9,1:ye))./d_l_tt(5,9,1),'-y', ...
     'Linewidth',2.5)
title('Low-skill service imports')
grid on

subplot(2,3,6)
plot(1995:1994+ye,squeeze(d_h_tt(5,1,1:ye))./d_h_tt(5,1,1),'-b',1995:1994+ye,squeeze(d_h_tt(5,2,1:ye))./d_h_tt(5,2,1),'-r',1995:1994+ye,squeeze(d_h_tt(5,3,1:ye))./d_h_tt(5,3,1),'-g',1995:1994+ye,squeeze(d_h_tt(5,4,1:ye))./d_h_tt(5,4,1),'-c',1995:1994+ye,squeeze(d_h_tt(5,7,1:ye))./d_h_tt(5,7,1),'-m',1995:1994+ye,squeeze(d_h_tt(5,9,1:ye))./d_h_tt(5,9,1),'-y', ...
    'Linewidth',2.5)
title('High-skill service imports')
grid on


sgtitle('Figure 3. Bilateral trade costs for India relative to 1995')
%sgtitle('Bilateral trade costs for Brazil relative to 2000')
saveas(gcf, 'plots_transition/final_2/Figure 3.png');
close(gcf)

figure('Position', get(0, 'Screensize'))
subplot(1,2,1)
plot(1995:1994+ye,squeeze(d_l_tt(2,5,1:ye))./d_l_tt(2,5,1),'-r',1995:1994+ye,squeeze(d_l_tt(8,5,1:ye))./d_l_tt(8,5,1),'-c','Linewidth',2.5)
title('Low-skill service exports')
legend('To BRA','To MEX','location','best') 
grid on

subplot(1,2,2)
plot(1995:1994+ye,squeeze(d_h_tt(2,5,1:ye))./d_h_tt(2,5,1),'-r',1995:1994+ye,squeeze(d_h_tt(8,5,1:ye))./d_h_tt(8,5,1),'-c', ...
     'Linewidth',2.5)
title('High-skill service exports')
grid on

sgtitle('Figure 3. Bilateral trade costs for India relative to 1995')
%sgtitle('Bilateral trade costs for Brazil relative to 2000')
saveas(gcf, 'plots_transition/final_2/Figure 3a.png');
close(gcf)

% % decompose: relative price
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,P_a_t(1,:)./P_a_t(2,:),'-r',2000:1999+ye,P_a_t(1,:)./P_a_t(3,:),'-g',2000:1999+ye,P_a_t(1,:)./P_a_t(4,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('From CHN','From USA','From ROW','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,P_m_t(1,:)./P_m_t(2,:),'-r',2000:1999+ye,P_m_t(1,:)./P_m_t(3,:),'-g',2000:1999+ye,P_m_t(1,:)./P_m_t(4,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,P_h_t(1,:)./P_h_t(2,:),'-r',2000:1999+ye,P_h_t(1,:)./P_h_t(3,:),'-g',2000:1999+ye,P_h_t(1,:)./P_h_t(4,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,P_s_t(1,:)./P_s_t(2,:),'-r',2000:1999+ye,P_s_t(1,:)./P_s_t(3,:),'-g',2000:1999+ye,P_s_t(1,:)./P_s_t(4,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Decompose Bilateral cost to Brazil: relative price')
% saveas(gcf, 'plots_transition/final_2/Bilateral cost to Brazil relative price.png');
% close(gcf)


% % decompose relative trade share
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(pi_a_tt(1,2,1:ye))./squeeze(pi_a_tt(2,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_a_tt(1,3,1:ye))./squeeze(pi_a_tt(3,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_a_tt(1,4,1:ye))./squeeze(pi_a_tt(4,4,1:ye)),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('From CHN','From USA','From ROW','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(pi_m_tt(1,2,1:ye))./squeeze(pi_m_tt(2,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_m_tt(1,3,1:ye))./squeeze(pi_m_tt(3,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_m_tt(1,4,1:ye))./squeeze(pi_m_tt(4,4,1:ye)),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(pi_h_tt(1,2,1:ye))./squeeze(pi_h_tt(2,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_h_tt(1,3,1:ye))./squeeze(pi_h_tt(3,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_h_tt(1,4,1:ye))./squeeze(pi_h_tt(4,4,1:ye)),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(pi_s_tt(1,2,1:ye))./squeeze(pi_s_tt(2,2,1:ye)),'-r',2000:1999+ye,squeeze(pi_s_tt(1,3,1:ye))./squeeze(pi_s_tt(3,3,1:ye)),'-g',2000:1999+ye,squeeze(pi_s_tt(1,4,1:ye))./squeeze(pi_s_tt(4,4,1:ye)),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Decompose Bilateral cost to Brazil: relative trade share')
% saveas(gcf, 'plots_transition/final_2/Bilateral cost to Brazil relative trade share.png');
% close(gcf)
% 
% 
% 
% 
% % Bilateral cost of shipping goods from Brazil 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,4,1)
% plot(2000:1999+ye,squeeze(d_g_tt(2,1,1:ye))./d_g_tt(2,1,1),'-r',2000:1999+ye,squeeze(d_g_tt(3,1,1:ye))./d_g_tt(3,1,1),'-g',2000:1999+ye,squeeze(d_g_tt(4,1,1:ye))./d_g_tt(4,1,1),'-c',2000:1999+ye,squeeze(d_g_tt(5,1,1:ye))./d_g_tt(5,1,1),'-m',...
%     'Linewidth',2.5)
% title('Primary exports')
% legend('To CHN','To USA','To IND','To ROW','location','best')
% grid on
% 
% subplot(2,4,2)
% plot(2000:1999+ye,squeeze(d_l_tt(2,1,1:ye))./d_l_tt(2,1,1),'-r',2000:1999+ye,squeeze(d_l_tt(3,1,1:ye))./d_l_tt(3,1,1),'-g',2000:1999+ye,squeeze(d_l_tt(4,1,1:ye))./d_l_tt(4,1,1),'-c',2000:1999+ye,squeeze(d_l_tt(5,1,1:ye))./d_l_tt(5,1,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing exports')
% grid on
% 
% subplot(2,4,3)
% plot(2000:1999+ye,squeeze(d_h_tt(2,1,1:ye))./d_h_tt(2,1,1),'-r',2000:1999+ye,squeeze(d_h_tt(3,1,1:ye))./d_h_tt(3,1,1),'-g',2000:1999+ye,squeeze(d_h_tt(4,1,1:ye))./d_h_tt(4,1,1),'-c',2000:1999+ye,squeeze(d_h_tt(5,1,1:ye))./d_h_tt(5,1,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing exports')
% grid on
% 
% subplot(2,4,4)
% plot(2000:1999+ye,squeeze(d_s_tt(2,1,1:ye))./d_s_tt(2,1,1),'-r',2000:1999+ye,squeeze(d_s_tt(3,1,1:ye))./d_s_tt(3,1,1),'-g',2000:1999+ye,squeeze(d_s_tt(4,1,1:ye))./d_s_tt(4,1,1),'-c',2000:1999+ye,squeeze(d_s_tt(5,1,1:ye))./d_s_tt(5,1,1),'-m',...
%     'Linewidth',2.5)
% title('Service exports')
% grid on
% 
% subplot(2,4,5)
% plot(2000:1999+ye,squeeze(d_g_tt(1,2,1:ye))./d_g_tt(1,2,1),'-r',2000:1999+ye,squeeze(d_g_tt(1,3,1:ye))./d_g_tt(1,3,1),'-g',2000:1999+ye,squeeze(d_g_tt(1,4,1:ye))./d_g_tt(1,4,1),'-c',2000:1999+ye,squeeze(d_g_tt(1,5,1:ye))./d_g_tt(1,5,1),'-m',...
%     'Linewidth',2.5)
% title('Primary imports')
% legend('From CHN','From USA','From IND','From ROW','location','best')
% grid on
% 
% subplot(2,4,6)
% plot(2000:1999+ye,squeeze(d_l_tt(1,2,1:ye))./d_l_tt(1,2,1),'-r',2000:1999+ye,squeeze(d_l_tt(1,3,1:ye))./d_l_tt(1,3,1),'-g',2000:1999+ye,squeeze(d_l_tt(1,4,1:ye))./d_l_tt(1,4,1),'-c',2000:1999+ye,squeeze(d_l_tt(1,5,1:ye))./d_l_tt(1,5,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing imports')
% grid on
% 
% subplot(2,4,7)
% plot(2000:1999+ye,squeeze(d_h_tt(1,2,1:ye))./d_h_tt(1,2,1),'-r',2000:1999+ye,squeeze(d_h_tt(1,3,1:ye))./d_h_tt(1,3,1),'-g',2000:1999+ye,squeeze(d_h_tt(1,4,1:ye))./d_h_tt(1,4,1),'-c',2000:1999+ye,squeeze(d_h_tt(1,5,1:ye))./d_h_tt(1,5,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing imports')
% grid on
% 
% subplot(2,4,8)
% plot(2000:1999+ye,squeeze(d_s_tt(1,2,1:ye))./d_s_tt(1,2,1),'-r',2000:1999+ye,squeeze(d_s_tt(1,3,1:ye))./d_s_tt(1,3,1),'-g',2000:1999+ye,squeeze(d_s_tt(1,4,1:ye))./d_s_tt(1,4,1),'-c',2000:1999+ye,squeeze(d_s_tt(1,5,1:ye))./d_s_tt(1,5,1),'-m',...
%     'Linewidth',2.5)
% title('Service imports')
% grid on
% 
% sgtitle('Figure 3. Bilateral trade costs for Brazil relative to 2000')
% %sgtitle('Bilateral trade costs for Brazil relative to 2000')
% saveas(gcf, 'plots_transition/final_2/Figure 3.png');
% close(gcf)
% 
% % decompose: relative price
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,P_a_t(2,:)./P_a_t(1,:),'-r',2000:1999+ye,P_a_t(3,:)./P_a_t(1,:),'-g',2000:1999+ye,P_a_t(4,:)./P_a_t(1,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('To CHN','To USA','To ROW','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,P_m_t(2,:)./P_m_t(1,:),'-r',2000:1999+ye,P_m_t(3,:)./P_m_t(1,:),'-g',2000:1999+ye,P_m_t(4,:)./P_m_t(1,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,P_h_t(2,:)./P_h_t(1,:),'-r',2000:1999+ye,P_h_t(3,:)./P_h_t(1,:),'-g',2000:1999+ye,P_h_t(4,:)./P_h_t(1,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,P_s_t(2,:)./P_s_t(1,:),'-r',2000:1999+ye,P_s_t(3,:)./P_s_t(1,:),'-g',2000:1999+ye,P_s_t(4,:)./P_s_t(1,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Decompose Bilateral cost from Brazil: relative price')
% saveas(gcf, 'plots_transition/final_2/Bilateral cost from Brazil relative price.png');
% close(gcf)
% 
% % decompose relative trade share
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(pi_a_tt(2,1,1:ye))./squeeze(pi_a_tt(1,1,1:ye)),'-r',2000:1999+ye,squeeze(pi_a_tt(3,1,1:ye))./squeeze(pi_a_tt(1,1,1:ye)),'-g',2000:1999+ye,squeeze(pi_a_tt(4,1,1:ye))./squeeze(pi_a_tt(1,1,1:ye)),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('To CHN','To USA','To ROW','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(pi_m_tt(2,1,1:ye))./squeeze(pi_m_tt(1,1,1:ye)),'-r',2000:1999+ye,squeeze(pi_m_tt(3,1,1:ye))./squeeze(pi_m_tt(1,1,1:ye)),'-g',2000:1999+ye,squeeze(pi_m_tt(4,1,1:ye))./squeeze(pi_m_tt(1,1,1:ye)),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(pi_h_tt(2,1,1:ye))./squeeze(pi_h_tt(1,1,1:ye)),'-r',2000:1999+ye,squeeze(pi_h_tt(3,1,1:ye))./squeeze(pi_h_tt(1,1,1:ye)),'-g',2000:1999+ye,squeeze(pi_h_tt(4,1,1:ye))./squeeze(pi_h_tt(1,1,1:ye)),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(pi_s_tt(2,1,1:ye))./squeeze(pi_s_tt(1,1,1:ye)),'-r',2000:1999+ye,squeeze(pi_s_tt(3,1,1:ye))./squeeze(pi_s_tt(1,1,1:ye)),'-g',2000:1999+ye,squeeze(pi_s_tt(4,1,1:ye))./squeeze(pi_s_tt(1,1,1:ye)),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Decompose Bilateral cost from Brazil: relative trade share')
% saveas(gcf, 'plots_transition/final_2/Bilateral cost from Brazil relative trade share.png');
% close(gcf)
% 
% 
% % Bilateral cost for CHN 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,4,1)
% plot(2000:1999+ye,squeeze(d_g_tt(1,2,1:ye))./d_g_tt(1,2,1),'-b',2000:1999+ye,squeeze(d_g_tt(3,2,1:ye))./d_g_tt(3,2,1),'-g',2000:1999+ye,squeeze(d_g_tt(4,2,1:ye))./d_g_tt(4,2,1),'-c',2000:1999+ye,squeeze(d_g_tt(5,2,1:ye))./d_g_tt(5,2,1),'-m',...
%      'Linewidth',2.5)
% title('Primary exports')
% legend('To BRA','To USA','To IND','To ROW','location','best')
% grid on
% 
% subplot(2,4,2)
% plot(2000:1999+ye,squeeze(d_l_tt(1,2,1:ye))./d_l_tt(1,2,1),'-b',2000:1999+ye,squeeze(d_l_tt(3,2,1:ye))./d_l_tt(3,2,1),'-g',2000:1999+ye,squeeze(d_l_tt(4,2,1:ye))./d_l_tt(4,2,1),'-c',2000:1999+ye,squeeze(d_l_tt(5,2,1:ye))./d_l_tt(5,2,1),'-m',...
%      'Linewidth',2.5)
% title('Low-tech manufacturing exports')
% grid on
% 
% subplot(2,4,3)
% plot(2000:1999+ye,squeeze(d_h_tt(1,2,1:ye))./d_h_tt(1,2,1),'-b',2000:1999+ye,squeeze(d_h_tt(3,2,1:ye))./d_h_tt(3,2,1),'-g',2000:1999+ye,squeeze(d_h_tt(4,2,1:ye))./d_h_tt(4,2,1),'-c',2000:1999+ye,squeeze(d_h_tt(5,2,1:ye))./d_h_tt(5,2,1),'-m',...
%      'Linewidth',2.5)
% title('High-tech manufacturing exports')
% grid on
% 
% subplot(2,4,4)
% plot(2000:1999+ye,squeeze(d_s_tt(1,2,1:ye))./d_s_tt(1,2,1),'-b',2000:1999+ye,squeeze(d_s_tt(3,2,1:ye))./d_s_tt(3,2,1),'-g',2000:1999+ye,squeeze(d_s_tt(4,2,1:ye))./d_s_tt(4,2,1),'-c',2000:1999+ye,squeeze(d_s_tt(5,2,1:ye))./d_s_tt(5,2,1),'-m',...
%      'Linewidth',2.5)
% title('Service exports')
% grid on
% 
% subplot(2,4,5)
% plot(2000:1999+ye,squeeze(d_g_tt(2,1,1:ye))./d_g_tt(2,1,1),'-b',2000:1999+ye,squeeze(d_g_tt(2,3,1:ye))./d_g_tt(2,3,1),'-g',2000:1999+ye,squeeze(d_g_tt(2,4,1:ye))./d_g_tt(2,4,1),'-c',2000:1999+ye,squeeze(d_g_tt(2,5,1:ye))./d_g_tt(2,5,1),'-m',...
%     'Linewidth',2.5)
% title('Primary imports')
% legend('From BRA','From USA','From IND','From ROW','location','best')
% grid on
% 
% subplot(2,4,6)
% plot(2000:1999+ye,squeeze(d_l_tt(2,1,1:ye))./d_l_tt(2,1,1),'-b',2000:1999+ye,squeeze(d_l_tt(2,3,1:ye))./d_l_tt(2,3,1),'-g',2000:1999+ye,squeeze(d_l_tt(2,4,1:ye))./d_l_tt(2,4,1),'-c',2000:1999+ye,squeeze(d_l_tt(2,5,1:ye))./d_l_tt(2,5,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing imports')
% grid on
% 
% subplot(2,4,7)
% plot(2000:1999+ye,squeeze(d_h_tt(2,1,1:ye))./d_h_tt(2,1,1),'-b',2000:1999+ye,squeeze(d_h_tt(2,3,1:ye))./d_h_tt(2,3,1),'-g',2000:1999+ye,squeeze(d_h_tt(2,4,1:ye))./d_h_tt(2,4,1),'-c',2000:1999+ye,squeeze(d_h_tt(2,5,1:ye))./d_h_tt(2,5,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing imports')
% grid on
% 
% subplot(2,4,8)
% plot(2000:1999+ye,squeeze(d_s_tt(2,1,1:ye))./d_s_tt(2,1,1),'-b',2000:1999+ye,squeeze(d_s_tt(2,3,1:ye))./d_s_tt(2,3,1),'-g',2000:1999+ye,squeeze(d_s_tt(2,4,1:ye))./d_s_tt(2,4,1),'-c',2000:1999+ye,squeeze(d_s_tt(2,5,1:ye))./d_s_tt(2,5,1),'-m',...
%     'Linewidth',2.5)
% title('Service imports')
% grid on
% 
% sgtitle('Figure 4. Bilateral trade costs for China relative to 2000')
% %sgtitle('Bilateral trade costs for China relative to 2000')
% saveas(gcf, 'plots_transition/final_2/Figure 4.png');
% close(gcf)
% 
% % Bilateral cost for IND 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,4,1)
% plot(2000:1999+ye,squeeze(d_g_tt(1,4,1:ye))./d_g_tt(1,4,1),'-b',2000:1999+ye,squeeze(d_g_tt(2,4,1:ye))./d_g_tt(2,4,1),'-r',2000:1999+ye,squeeze(d_g_tt(3,4,1:ye))./d_g_tt(3,4,1),'-g',2000:1999+ye,squeeze(d_g_tt(5,4,1:ye))./d_g_tt(5,4,1),'-m',...
%      'Linewidth',2.5)
% title('Primary exports')
% legend('To BRA','To CHN','To USA','To ROW','location','best')
% grid on
% 
% subplot(2,4,2)
% plot(2000:1999+ye,squeeze(d_l_tt(1,4,1:ye))./d_l_tt(1,4,1),'-b',2000:1999+ye,squeeze(d_l_tt(2,4,1:ye))./d_l_tt(2,4,1),'-r',2000:1999+ye,squeeze(d_l_tt(3,4,1:ye))./d_l_tt(3,4,1),'-g',2000:1999+ye,squeeze(d_l_tt(5,4,1:ye))./d_l_tt(5,4,1),'-m',...
%      'Linewidth',2.5)
% title('Low-tech manufacturing exports')
% grid on
% 
% subplot(2,4,3)
% plot(2000:1999+ye,squeeze(d_h_tt(1,4,1:ye))./d_h_tt(1,4,1),'-b',2000:1999+ye,squeeze(d_h_tt(2,4,1:ye))./d_h_tt(2,4,1),'-r',2000:1999+ye,squeeze(d_h_tt(3,4,1:ye))./d_h_tt(3,4,1),'-g',2000:1999+ye,squeeze(d_h_tt(5,4,1:ye))./d_h_tt(5,4,1),'-m',...
%      'Linewidth',2.5)
% title('High-tech manufacturing exports')
% grid on
% 
% subplot(2,4,4)
% plot(2000:1999+ye,squeeze(d_s_tt(1,4,1:ye))./d_s_tt(1,4,1),'-b',2000:1999+ye,squeeze(d_s_tt(2,4,1:ye))./d_s_tt(2,4,1),'-r',2000:1999+ye,squeeze(d_s_tt(3,4,1:ye))./d_s_tt(3,4,1),'-g',2000:1999+ye,squeeze(d_s_tt(5,4,1:ye))./d_s_tt(5,4,1),'-m',...
%      'Linewidth',2.5)
% title('Service exports')
% grid on
% 
% subplot(2,4,5)
% plot(2000:1999+ye,squeeze(d_g_tt(4,1,1:ye))./d_g_tt(4,1,1),'-b',2000:1999+ye,squeeze(d_g_tt(4,2,1:ye))./d_g_tt(4,2,1),'-r',2000:1999+ye,squeeze(d_g_tt(4,3,1:ye))./d_g_tt(4,3,1),'-g',2000:1999+ye,squeeze(d_g_tt(4,5,1:ye))./d_g_tt(4,5,1),'-m',...
%     'Linewidth',2.5)
% title('Primary imports')
% legend('From BRA','From CHN','From USA','From ROW','location','best')
% grid on
% 
% subplot(2,4,6)
% plot(2000:1999+ye,squeeze(d_l_tt(4,1,1:ye))./d_l_tt(4,1,1),'-b',2000:1999+ye,squeeze(d_l_tt(4,2,1:ye))./d_l_tt(4,2,1),'-r',2000:1999+ye,squeeze(d_l_tt(4,3,1:ye))./d_l_tt(4,3,1),'-g',2000:1999+ye,squeeze(d_l_tt(4,5,1:ye))./d_l_tt(4,5,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing imports')
% grid on
% 
% subplot(2,4,7)
% plot(2000:1999+ye,squeeze(d_h_tt(4,1,1:ye))./d_h_tt(4,1,1),'-b',2000:1999+ye,squeeze(d_h_tt(4,2,1:ye))./d_h_tt(4,2,1),'-r',2000:1999+ye,squeeze(d_h_tt(4,3,1:ye))./d_h_tt(4,3,1),'-g',2000:1999+ye,squeeze(d_h_tt(4,5,1:ye))./d_h_tt(4,5,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing imports')
% grid on
% 
% subplot(2,4,8)
% plot(2000:1999+ye,squeeze(d_s_tt(4,1,1:ye))./d_s_tt(4,1,1),'-b',2000:1999+ye,squeeze(d_s_tt(4,2,1:ye))./d_s_tt(4,2,1),'-r',2000:1999+ye,squeeze(d_s_tt(4,3,1:ye))./d_s_tt(4,3,1),'-g',2000:1999+ye,squeeze(d_s_tt(4,5,1:ye))./d_s_tt(4,5,1),'-m',...
%     'Linewidth',2.5)
% title('Service imports')
% grid on
% 
% %sgtitle('Figure. Bilateral trade costs for India relative to 2000')
% sgtitle('Bilateral trade costs for India relative to 2000')
% saveas(gcf, 'plots_transition/final_2/Trade_cost_India.png');
% close(gcf)
% 
% % Bilateral cost for USA 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,4,1)
% plot(2000:1999+ye,squeeze(d_g_tt(1,3,1:ye))./d_g_tt(1,3,1),'-b',2000:1999+ye,squeeze(d_g_tt(2,3,1:ye))./d_g_tt(2,3,1),'-r',2000:1999+ye,squeeze(d_g_tt(4,3,1:ye))./d_g_tt(4,3,1),'-c',2000:1999+ye,squeeze(d_g_tt(5,3,1:ye))./d_g_tt(5,3,1),'-m',...
%      'Linewidth',2.5)
% title('Primary exports')
% legend('To BRA','To CHN','To IND','To ROW','location','best')
% grid on
% 
% subplot(2,4,2)
% plot(2000:1999+ye,squeeze(d_l_tt(1,3,1:ye))./d_l_tt(1,3,1),'-b',2000:1999+ye,squeeze(d_l_tt(2,3,1:ye))./d_l_tt(2,3,1),'-r',2000:1999+ye,squeeze(d_l_tt(4,3,1:ye))./d_l_tt(4,3,1),'-c',2000:1999+ye,squeeze(d_l_tt(5,3,1:ye))./d_l_tt(5,3,1),'-m',...
%      'Linewidth',2.5)
% title('Low-tech manufacturing exports')
% grid on
% 
% subplot(2,4,3)
% plot(2000:1999+ye,squeeze(d_h_tt(1,3,1:ye))./d_h_tt(1,3,1),'-b',2000:1999+ye,squeeze(d_h_tt(2,3,1:ye))./d_h_tt(2,3,1),'-r',2000:1999+ye,squeeze(d_h_tt(4,3,1:ye))./d_h_tt(4,3,1),'-c',2000:1999+ye,squeeze(d_h_tt(5,3,1:ye))./d_h_tt(5,3,1),'-m',...
%      'Linewidth',2.5)
% title('High-tech manufacturing exports')
% grid on
% 
% subplot(2,4,4)
% plot(2000:1999+ye,squeeze(d_s_tt(1,3,1:ye))./d_s_tt(1,3,1),'-b',2000:1999+ye,squeeze(d_s_tt(2,3,1:ye))./d_s_tt(2,3,1),'-r',2000:1999+ye,squeeze(d_s_tt(4,3,1:ye))./d_s_tt(4,3,1),'-c',2000:1999+ye,squeeze(d_s_tt(5,3,1:ye))./d_s_tt(5,3,1),'-m',...
%      'Linewidth',2.5)
% title('Service exports')
% grid on
% 
% subplot(2,4,5)
% plot(2000:1999+ye,squeeze(d_g_tt(3,1,1:ye))./d_g_tt(3,1,1),'-b',2000:1999+ye,squeeze(d_g_tt(3,2,1:ye))./d_g_tt(3,2,1),'-r',2000:1999+ye,squeeze(d_g_tt(3,4,1:ye))./d_g_tt(3,4,1),'-c',2000:1999+ye,squeeze(d_g_tt(3,5,1:ye))./d_g_tt(3,5,1),'-m',...
%     'Linewidth',2.5)
% title('Primary imports')
% legend('From BRA','From CHN','From IND','From ROW','location','best')
% grid on
% 
% subplot(2,4,6)
% plot(2000:1999+ye,squeeze(d_l_tt(3,1,1:ye))./d_l_tt(3,1,1),'-b',2000:1999+ye,squeeze(d_l_tt(3,2,1:ye))./d_l_tt(3,2,1),'-r',2000:1999+ye,squeeze(d_l_tt(3,4,1:ye))./d_l_tt(3,4,1),'-c',2000:1999+ye,squeeze(d_l_tt(3,5,1:ye))./d_l_tt(3,5,1),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing imports')
% grid on
% 
% subplot(2,4,7)
% plot(2000:1999+ye,squeeze(d_h_tt(3,1,1:ye))./d_h_tt(3,1,1),'-b',2000:1999+ye,squeeze(d_h_tt(3,2,1:ye))./d_h_tt(3,2,1),'-r',2000:1999+ye,squeeze(d_h_tt(3,4,1:ye))./d_h_tt(3,4,1),'-c',2000:1999+ye,squeeze(d_h_tt(3,5,1:ye))./d_h_tt(3,5,1),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing imports')
% grid on
% 
% subplot(2,4,8)
% plot(2000:1999+ye,squeeze(d_s_tt(3,1,1:ye))./d_s_tt(3,1,1),'-b',2000:1999+ye,squeeze(d_s_tt(3,2,1:ye))./d_s_tt(3,2,1),'-r',2000:1999+ye,squeeze(d_s_tt(3,4,1:ye))./d_s_tt(3,4,1),'-c',2000:1999+ye,squeeze(d_s_tt(3,5,1:ye))./d_s_tt(3,5,1),'-m',...
%     'Linewidth',2.5)
% title('Service imports')
% grid on
% 
% sgtitle('Figure. Bilateral trade costs for USA relative to 2000')
% saveas(gcf, 'plots_transition/final_2/Trade_cost_USA.png');
% close(gcf)
% 
% % Bilateral cost for ROW 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,4,1)
% plot(2000:1999+ye,squeeze(d_g_tt(1,5,1:ye))./d_g_tt(1,5,1),'-b',2000:1999+ye,squeeze(d_g_tt(2,5,1:ye))./d_g_tt(2,5,1),'-r',2000:1999+ye,squeeze(d_g_tt(3,5,1:ye))./d_g_tt(3,5,1),'-g',2000:1999+ye,squeeze(d_g_tt(4,5,1:ye))./d_g_tt(4,5,1),'-c',...
%      'Linewidth',2.5)
% title('Primary exports')
% legend('To BRA','To CHN','To USA','To IND','location','best')
% grid on
% 
% subplot(2,4,2)
% plot(2000:1999+ye,squeeze(d_l_tt(1,5,1:ye))./d_l_tt(1,5,1),'-b',2000:1999+ye,squeeze(d_l_tt(2,5,1:ye))./d_l_tt(2,5,1),'-r',2000:1999+ye,squeeze(d_l_tt(3,5,1:ye))./d_l_tt(3,5,1),'-g',2000:1999+ye,squeeze(d_l_tt(4,5,1:ye))./d_l_tt(4,5,1),'-c',...
%      'Linewidth',2.5)
% title('Low-tech manufacturing exports')
% grid on
% 
% subplot(2,4,3)
% plot(2000:1999+ye,squeeze(d_h_tt(1,5,1:ye))./d_h_tt(1,5,1),'-b',2000:1999+ye,squeeze(d_h_tt(2,5,1:ye))./d_h_tt(2,5,1),'-r',2000:1999+ye,squeeze(d_h_tt(3,5,1:ye))./d_h_tt(3,5,1),'-g',2000:1999+ye,squeeze(d_h_tt(4,5,1:ye))./d_h_tt(4,5,1),'-c',...
%      'Linewidth',2.5)
% title('High-tech manufacturing exports')
% grid on
% 
% subplot(2,4,4)
% plot(2000:1999+ye,squeeze(d_s_tt(1,5,1:ye))./d_s_tt(1,5,1),'-b',2000:1999+ye,squeeze(d_s_tt(2,5,1:ye))./d_s_tt(2,5,1),'-r',2000:1999+ye,squeeze(d_s_tt(3,5,1:ye))./d_s_tt(3,5,1),'-g',2000:1999+ye,squeeze(d_s_tt(4,5,1:ye))./d_s_tt(4,5,1),'-c',...
%      'Linewidth',2.5)
% title('Service exports')
% grid on
% 
% subplot(2,4,5)
% plot(2000:1999+ye,squeeze(d_g_tt(5,1,1:ye))./d_g_tt(5,1,1),'-b',2000:1999+ye,squeeze(d_g_tt(5,2,1:ye))./d_g_tt(5,2,1),'-r',2000:1999+ye,squeeze(d_g_tt(5,3,1:ye))./d_g_tt(5,3,1),'-g',2000:1999+ye,squeeze(d_g_tt(5,4,1:ye))./d_g_tt(5,4,1),'-c',...
%     'Linewidth',2.5)
% title('Primary imports')
% legend('From BRA','From CHN','From USA','From IND','location','best')
% grid on
% 
% subplot(2,4,6)
% plot(2000:1999+ye,squeeze(d_l_tt(5,1,1:ye))./d_l_tt(5,1,1),'-b',2000:1999+ye,squeeze(d_l_tt(5,2,1:ye))./d_l_tt(5,2,1),'-r',2000:1999+ye,squeeze(d_l_tt(5,3,1:ye))./d_l_tt(5,3,1),'-g',2000:1999+ye,squeeze(d_l_tt(5,4,1:ye))./d_l_tt(5,4,1),'-c',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing imports')
% grid on
% 
% subplot(2,4,7)
% plot(2000:1999+ye,squeeze(d_h_tt(5,1,1:ye))./d_h_tt(5,1,1),'-b',2000:1999+ye,squeeze(d_h_tt(5,2,1:ye))./d_h_tt(5,2,1),'-r',2000:1999+ye,squeeze(d_h_tt(5,3,1:ye))./d_h_tt(5,3,1),'-g',2000:1999+ye,squeeze(d_h_tt(5,4,1:ye))./d_h_tt(5,4,1),'-c',...
%     'Linewidth',2.5)
% title('High-tech manufacturing imports')
% grid on
% 
% subplot(2,4,8)
% plot(2000:1999+ye,squeeze(d_s_tt(5,1,1:ye))./d_s_tt(5,1,1),'-b',2000:1999+ye,squeeze(d_s_tt(5,2,1:ye))./d_s_tt(5,2,1),'-r',2000:1999+ye,squeeze(d_s_tt(5,3,1:ye))./d_s_tt(5,3,1),'-g',2000:1999+ye,squeeze(d_s_tt(5,4,1:ye))./d_s_tt(5,4,1),'-c',...
%     'Linewidth',2.5)
% title('Service imports')
% grid on
% 
% sgtitle('Figure. Bilateral trade costs for ROW relative to 2000')
% saveas(gcf, 'plots_transition/final_2/Trade_cost_ROW.png');
% close(gcf)
 
% Static global portfolio 2000-2014
figure('Position', get(0, 'Screensize'))
subplot(1,2,1)
plot(1995:1994+ye,NX_portfolio_t(1,:)./GDP_t(1,:),'-b',1995:1994+ye,NX_portfolio_t(2,:)./GDP_t(2,:),'-r',1995:1994+ye,NX_portfolio_t(3,:)./GDP_t(3,:),'-g',1995:1994+ye,NX_portfolio_t(4,:)./GDP_t(4,:),'-c',1995:1994+ye,NX_portfolio_t(5,:)./GDP_t(5,:),'-k',1995:1994+ye,NX_portfolio_t(7,:)./GDP_t(7,:),'-m',1995:1994+ye,NX_portfolio_t(9,:)./GDP_t(9,:),'-y', ...
     1995:1994+ye,NX_t(1,:)./GDP_t(1,:),'--b',1995:1994+ye,NX_t(2,:)./GDP_t(2,:),'--r',1995:1994+ye,NX_t(3,:)./GDP_t(3,:),'--g',1995:1994+ye,NX_t(4,:)./GDP_t(4,:),'--c',1995:1994+ye,NX_t(5,:)./GDP_t(5,:),'--k',1995:1994+ye,NX_t(7,:)./GDP_t(7,:),'--m',1995:1994+ye,NX_t(9,:)./GDP_t(9,:),'--y', ...
     'Linewidth',2.5)
title('Net exports over GDP')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','northeast') 
%legend('BRA', 'CHN', 'USA','IND','ROW','Target','location','northeast') 
%ylim([-0.08 0.12])
grid on

subplot(1,2,2)
plot(1995:1994+ye,fi_t(1,:),'-b',1995:1994+ye,fi_t(2,:),'-r',1995:1994+ye,fi_t(3,:),'-g',1995:1994+ye,fi_t(4,:),'-c',1995:1994+ye,fi_t(5,:),'-k',1995:1994+ye,fi_t(7,:),'-m',1995:1994+ye,fi_t(9,:),'-y', ...
     'Linewidth',2.5)
title('Static global portfolio')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','northeast') 
grid on

sgtitle('Figure 5. Global portfolio share 1995-2011')
%sgtitle('Global portfolio share 2000-2014')
saveas(gcf, 'plots_transition/final_2/Figure 5.png');
close(gcf)

% % Fitting sector price relative to USA 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,P_a_tt(1,1:ye)./P_a_tt(3,1:ye),'-b',2000:1999+ye,P_a_tt(2,1:ye)./P_a_tt(3,1:ye),'-r',2000:1999+ye,P_a_tt(4,1:ye)./P_a_tt(3,1:ye),'-m',...
%     2000:2000+num_year-1,P_a_t([1,2,4],1:ye)./P_a_t(3,1:ye),'--black','Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','ROW','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,P_m_tt(1,1:ye)./P_m_tt(3,1:ye),'-b',2000:1999+ye,P_m_tt(2,1:ye)./P_m_tt(3,1:ye),'-r',2000:1999+ye,P_m_tt(4,1:ye)./P_m_tt(3,1:ye),'-m',...
%     2000:2000+num_year-1,P_m_t([1,2,4],1:ye)./P_m_t(3,1:ye),'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,P_h_tt(1,1:ye)./P_h_tt(3,1:ye),'-b',2000:1999+ye,P_h_tt(2,1:ye)./P_h_tt(3,1:ye),'-r',2000:1999+ye,P_h_tt(4,1:ye)./P_h_tt(3,1:ye),'-m',...
%     2000:2000+num_year-1,P_h_t([1,2,4],1:ye)./P_h_t(3,1:ye),'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,P_s_tt(1,1:ye)./P_s_tt(3,1:ye),'-b',2000:1999+ye,P_s_tt(2,1:ye)./P_s_tt(3,1:ye),'-r',2000:1999+ye,P_s_tt(4,1:ye)./P_s_tt(3,1:ye),'-m',...
%     2000:2000+num_year-1,P_s_t([1,2,4],1:ye)./P_s_t(3,1:ye),'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Fitting sector price relative to USA 2000-2014')
% saveas(gcf, 'plots_transition/final_2/Fitting sector price relative to USA 2000-2014.png');
% close(gcf)


% untargeted capital stock PWT
figure('Position', get(0, 'Screensize'))
subplot(1,2,1)
plot(1995:1994+ye,K_tt(1,1:num_year)./L_tt(1,1:num_year),'-b',1995:1994+ye,K_tt(2,1:num_year)./L_tt(2,1:num_year),'-r',1995:1994+ye,K_tt(3,1:num_year)./L_tt(3,1:num_year),'-g',1995:1994+ye,K_tt(4,1:num_year)./L_tt(4,1:num_year),'-c',1995:1994+ye,K_tt(5,1:num_year)./L_tt(5,1:num_year),'-k',1995:1994+ye,K_tt(7,1:num_year)./L_tt(7,1:num_year),'-m',1995:1994+ye,K_tt(9,1:num_year)./L_tt(9,1:num_year),'-y', ...
     1995:1994+ye,K_t(1,1:num_year)./L_t(1,1:num_year),'--b',1995:1994+ye,K_t(2,1:num_year)./L_t(2,1:num_year),'--r',1995:1994+ye,K_t(3,1:num_year)./L_t(3,1:num_year),'--g',1995:1994+ye,K_t(4,1:num_year)./L_t(4,1:num_year),'--c',1995:1994+ye,K_t(5,1:num_year)./L_t(5,1:num_year),'--k',1995:1994+ye,K_t(7,1:num_year)./L_t(7,1:num_year),'--m',1995:1994+ye,K_t(9,1:num_year)./L_t(9,1:num_year),'--y', ...
     'Linewidth',2.5)
title('Capital stock per worker')
legend('AUS','BRA','CHN','DEU','IND','MEX','USA','location','northeast') 
grid on

subplot(1,2,2)
plot(1995:1994+ye,K_tt(1,1:num_year),'-b',1995:1994+ye,K_tt(2,1:num_year),'-r',1995:1994+ye,K_tt(3,1:num_year),'-g',1995:1994+ye,K_tt(4,1:num_year),'-c',1995:1994+ye,K_tt(5,1:num_year),'-k',1995:1994+ye,K_tt(7,1:num_year),'-m',1995:1994+ye,K_tt(9,1:num_year),'-y', ...
     1995:1994+ye,K_t(1,1:num_year),'--b',1995:1994+ye,K_t(2,1:num_year),'--r',1995:1994+ye,K_t(3,1:num_year),'--g',1995:1994+ye,K_t(4,1:num_year),'--c',1995:1994+ye,K_t(5,1:num_year),'--k',1995:1994+ye,K_t(7,1:num_year),'--m',1995:1994+ye,K_t(9,1:num_year),'--y', ...
     'Linewidth',2.5)
title('Capital stock')
grid on

sgtitle('Figure 7. Capital stock by country 1995-2011')
%sgtitle('Capital stock by country 2000-2014 PWT')
saveas(gcf, 'plots_transition\final_2\Figure 7.png');
close(gcf)

% untargeted capital stock
% figure('Position', get(0, 'Screensize'))
% subplot(1,2,1)
% plot(2000:1999+num_year,K_t(1,1:num_year)./L_t(1,1:num_year),'-b',2000:1999+num_year,K_t(2,1:num_year)./L_t(2,1:num_year),'-r',2000:1999+num_year,K_t(3,1:num_year)./L_t(3,1:num_year),'-g',2000:1999+num_year,K_t(4,1:num_year)./L_t(4,1:num_year),'-c',2000:1999+num_year,K_t(5,1:num_year)./L_t(5,1:num_year),'-m','Linewidth',2.5)
% title('Capital stock per worker')
% grid on
% hold on
% plot(2000:1999+8,k_t(1,1:8)./L_t(1,1:8),'--b',2000:1999+8,k_t(2,1:8)./L_t(2,1:8),'--r',2000:1999+8,k_t(3,1:8)./L_t(3,1:8),'--g',2000:1999+8,k_t(4,1:8)./L_t(4,1:8),'--c',2000:1999+8,k_t(5,1:8)./L_t(5,1:8),'--m','Linewidth',2.5)
% hold off
% legend('BRA', 'CHN', 'USA','IND','ROW','Raw data','location','best')
% 
% subplot(1,2,2)
% plot(2000:1999+num_year,K_t(1,1:num_year),'-b',2000:1999+num_year,K_t(2,1:num_year),'-r',2000:1999+num_year,K_t(3,1:num_year),'-g',2000:1999+num_year,K_t(4,1:num_year),'-c',2000:1999+num_year,K_t(5,1:num_year),'-m','Linewidth',2.5)
% title('Capital stock')
% grid on
% hold on
% plot(2000:1999+8,k_t(1,1:8),'--b',2000:1999+8,k_t(2,1:8),'--r',2000:1999+8,k_t(3,1:8),'--g',2000:1999+8,k_t(4,1:8),'--c',2000:1999+8,k_t(5,1:8),'--m','Linewidth',2.5)
% hold off
% 
% %sgtitle('Figure 6. Capital stock by country 2000-2014')
% sgtitle('Capital stock by country 2000-2014')
% %saveas(gcf, 'plots_transition\final_2\Figure 6.png');
% close(gcf)
% 
% % untargeted capital stock PWT
% figure('Position', get(0, 'Screensize'))
% subplot(1,2,1)
% plot(2000:1999+num_year,K_tt(1,1:num_year)./L_tt(1,1:num_year),'-b',2000:1999+num_year,K_tt(2,1:num_year)./L_tt(2,1:num_year),'-r',2000:1999+num_year,K_tt(3,1:num_year)./L_tt(3,1:num_year),'-g',2000:1999+num_year,K_tt(4,1:num_year)./L_tt(4,1:num_year),'-c',2000:1999+num_year,K_tt(5,1:num_year)./L_tt(5,1:num_year),'-m','Linewidth',2.5)
% title('Capital stock per worker')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(1,1:num_year)./L_t(1,1:num_year),'--b',2000:1999+num_year,real_K_t_simul_PWT(2,1:num_year)./L_t(2,1:num_year),'--r',2000:1999+num_year,real_K_t_simul_PWT(3,1:num_year)./L_t(3,1:num_year),'--g',2000:1999+num_year,real_K_t_simul_PWT(4,1:num_year)./L_t(4,1:num_year),'--c',2000:1999+num_year,real_K_t_simul_PWT(5,1:num_year)./L_t(5,1:num_year),'--m','Linewidth',2.5)
% hold off
% legend('BRA', 'CHN', 'USA','IND','ROW','Raw data','location','best')
% 
% subplot(1,2,2)
% plot(2000:1999+num_year,K_tt(1,1:num_year),'-b',2000:1999+num_year,K_tt(2,1:num_year),'-r',2000:1999+num_year,K_tt(3,1:num_year),'-g',2000:1999+num_year,K_tt(4,1:num_year),'-c',2000:1999+num_year,K_tt(5,1:num_year),'-m','Linewidth',2.5)
% title('Capital stock')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(1,1:num_year),'--b',2000:1999+num_year,real_K_t_simul_PWT(2,1:num_year),'--r',2000:1999+num_year,real_K_t_simul_PWT(3,1:num_year),'--g',2000:1999+num_year,real_K_t_simul_PWT(4,1:num_year),'--c',2000:1999+num_year,real_K_t_simul_PWT(5,1:num_year),'--m','Linewidth',2.5)
% hold off
% 
% sgtitle('Figure 7. Capital stock by country 2000-2014')
% %sgtitle('Capital stock by country 2000-2014 PWT')
% saveas(gcf, 'plots_transition\final_2\Figure 7.png');
% close(gcf)
% 
% % untargeted capital stock PWT
% figure('Position', get(0, 'Screensize'))
% subplot(1,2,1)
% plot(2000:1999+num_year,K_tt(1,1:num_year)./L_tt(1,1:num_year)/(K_tt(1,1)./L_tt(1,1)),'-b',2000:1999+num_year,K_tt(2,1:num_year)./L_tt(2,1:num_year)/(K_tt(2,1)./L_tt(2,1)),'-r',2000:1999+num_year,K_tt(3,1:num_year)./L_tt(3,1:num_year)/(K_tt(3,1)./L_tt(3,1)),'-g',2000:1999+num_year,K_tt(4,1:num_year)./L_tt(4,1:num_year)/(K_tt(4,1)./L_tt(4,1)),'-c',2000:1999+num_year,K_tt(5,1:num_year)./L_tt(5,1:num_year)/(K_tt(5,1)./L_tt(5,1)),'-m','Linewidth',2.5)
% title('Capital stock per worker')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(1,1:num_year)./L_t(1,1:num_year)/(real_K_t_simul_PWT(1,1)./L_t(1,1)),'--b',2000:1999+num_year,real_K_t_simul_PWT(2,1:num_year)./L_t(2,1:num_year)/(real_K_t_simul_PWT(2,1)./L_t(2,1)),'--r',2000:1999+num_year,real_K_t_simul_PWT(3,1:num_year)./L_t(3,1:num_year)/(real_K_t_simul_PWT(3,1)./L_t(3,1)),'--g',2000:1999+num_year,real_K_t_simul_PWT(4,1:num_year)./L_t(4,1:num_year)/(real_K_t_simul_PWT(4,1)./L_t(4,1)),'--c',2000:1999+num_year,real_K_t_simul_PWT(5,1:num_year)./L_t(5,1:num_year)/(real_K_t_simul_PWT(5,1)./L_t(5,1)),'--m','Linewidth',2.5)
% hold off
% legend('BRA', 'CHN', 'USA','IND','ROW','Raw data','location','best')
% 
% subplot(1,2,2)
% plot(2000:1999+num_year,K_tt(1,1:num_year)/K_tt(1,1),'-b',2000:1999+num_year,K_tt(2,1:num_year)/K_tt(2,1),'-r',2000:1999+num_year,K_tt(3,1:num_year)/K_tt(3,1),'-g',2000:1999+num_year,K_tt(4,1:num_year)/K_tt(4,1),'-c',2000:1999+num_year,K_tt(5,1:num_year)/K_tt(5,1),'-m','Linewidth',2.5)
% title('Capital stock')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(1,1:num_year)/real_K_t_simul_PWT(1,1),'--b',2000:1999+num_year,real_K_t_simul_PWT(2,1:num_year)/real_K_t_simul_PWT(2,1),'--r',2000:1999+num_year,real_K_t_simul_PWT(3,1:num_year)/real_K_t_simul_PWT(3,1),'--g',2000:1999+num_year,real_K_t_simul_PWT(4,1:num_year)/real_K_t_simul_PWT(4,1),'--c',2000:1999+num_year,real_K_t_simul_PWT(5,1:num_year)/real_K_t_simul_PWT(5,1),'--m','Linewidth',2.5)
% hold off
% 
% sgtitle('Figure. Capital stock by country relative to 2000')
% %sgtitle('Capital stock by country relative to 2000')
% saveas(gcf, 'plots_transition\final_2\Figure 6b.png');
% close(gcf)
% 
% % untargeted capital stock PWT
% figure('Position', get(0, 'Screensize'))
% subplot(2,3,1)
% plot(2000:1999+num_year,K_tt(1,1:num_year),'-b','Linewidth',2.5)
% title('Brazil')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(1,1:num_year),'--b','Linewidth',2.5)
% hold off
% 
% subplot(2,3,2)
% plot(2000:1999+num_year,K_tt(2,1:num_year),'-r','Linewidth',2.5)
% title('China')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(2,1:num_year),'--r','Linewidth',2.5)
% hold off
% 
% subplot(2,3,3)
% plot(2000:1999+num_year,K_tt(3,1:num_year),'-g','Linewidth',2.5)
% title('USA')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(3,1:num_year),'--g','Linewidth',2.5)
% hold off
% 
% subplot(2,3,4)
% plot(2000:1999+num_year,K_tt(4,1:num_year),'-c','Linewidth',2.5)
% title('IND')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(4,1:num_year),'--c','Linewidth',2.5)
% hold off
% 
% subplot(2,3,5)
% plot(2000:1999+num_year,K_tt(5,1:num_year),'-m','Linewidth',2.5)
% title('ROW')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(5,1:num_year),'--m','Linewidth',2.5)
% hold off
% 
% %sgtitle('Figure 7. Implied capital stock by country 2000-2014')
% sgtitle('Implied capital stock by country 2000-2014 PWT')
% saveas(gcf, 'plots_transition\final_2\Figure 7c.png');
% close(gcf)
% 
% % untargeted capital stock PWT
% figure('Position', get(0, 'Screensize'))
% subplot(2,3,1)
% plot(2000:1999+num_year,K_tt(1,1:num_year)./L_tt(1,1:num_year),'-b','Linewidth',2.5)
% title('Brazil')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(1,1:num_year)./L_t(1,1:num_year),'--b','Linewidth',2.5)
% hold off
% 
% subplot(2,3,2)
% plot(2000:1999+num_year,K_tt(2,1:num_year)./L_tt(2,1:num_year),'-r','Linewidth',2.5)
% title('China')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(2,1:num_year)./L_t(2,1:num_year),'--r','Linewidth',2.5)
% hold off
% 
% subplot(2,3,3)
% plot(2000:1999+num_year,K_tt(3,1:num_year)./L_tt(3,1:num_year),'-g','Linewidth',2.5)
% title('USA')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(3,1:num_year)./L_t(3,1:num_year),'--g','Linewidth',2.5)
% hold off
% 
% subplot(2,3,4)
% plot(2000:1999+num_year,K_tt(4,1:num_year)./L_tt(4,1:num_year),'-c','Linewidth',2.5)
% title('IND')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(4,1:num_year)./L_t(4,1:num_year),'--c','Linewidth',2.5)
% hold off
% 
% subplot(2,3,5)
% plot(2000:1999+num_year,K_tt(5,1:num_year)./L_tt(5,1:num_year),'-m','Linewidth',2.5)
% title('ROW')
% grid on
% hold on
% plot(2000:1999+num_year,real_K_t_simul_PWT(5,1:num_year)./L_t(5,1:num_year),'--m','Linewidth',2.5)
% hold off
% 
% %sgtitle('Figure 7. Implied capital stock by country 2000-2014')
% sgtitle('Implied capital per worker by country 2000-2014 PWT')
% saveas(gcf, 'plots_transition\final_2\Figure 7d.png');
% close(gcf)
% 
% % untargeted sector labor share 
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,L_g_tt(1,1:ye)./L_tt(1,1:ye),'-b',2000:1999+ye,L_g_tt(2,1:ye)./L_tt(2,1:ye),'-r',2000:1999+ye,L_g_tt(3,1:ye)./L_tt(3,1:ye),'-g',2000:1999+ye,L_g_tt(4,1:ye)./L_tt(4,1:ye),'-c',2000:1999+ye,L_g_tt(5,1:ye)./L_tt(5,1:ye),'-m',...
%     2000:2000+num_year-1,lh_g_t,'--black','Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','IND','ROW','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,L_l_tt(1,1:ye)./L_tt(1,1:ye),'-b',2000:1999+ye,L_l_tt(2,1:ye)./L_tt(2,1:ye),'-r',2000:1999+ye,L_l_tt(3,1:ye)./L_tt(3,1:ye),'-g',2000:1999+ye,L_l_tt(4,1:ye)./L_tt(4,1:ye),'-c',2000:1999+ye,L_l_tt(5,1:ye)./L_tt(5,1:ye),'-m',...
%     2000:2000+num_year-1,lh_l_t,'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,L_h_tt(1,1:ye)./L_tt(1,1:ye),'-b',2000:1999+ye,L_h_tt(2,1:ye)./L_tt(2,1:ye),'-r',2000:1999+ye,L_h_tt(3,1:ye)./L_tt(3,1:ye),'-g',2000:1999+ye,L_h_tt(4,1:ye)./L_tt(4,1:ye),'-c',2000:1999+ye,L_h_tt(5,1:ye)./L_tt(5,1:ye),'-m',...
%     2000:2000+num_year-1,lh_h_t,'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,L_s_tt(1,1:ye)./L_tt(1,1:ye),'-b',2000:1999+ye,L_s_tt(2,1:ye)./L_tt(2,1:ye),'-r',2000:1999+ye,L_s_tt(3,1:ye)./L_tt(3,1:ye),'-g',2000:1999+ye,L_s_tt(4,1:ye)./L_tt(4,1:ye),'-c',2000:1999+ye,L_s_tt(5,1:ye)./L_tt(5,1:ye),'-m',...
%     2000:2000+num_year-1,lh_s_t,'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Untargeted sector labor share 2000-2014')
% saveas(gcf, 'plots_transition/Untargeted sector labor share 2000-2014.png');
% close(gcf)
% 
% % untargeted capital stock PWT
% figure('Position', get(0, 'Screensize'))
% subplot(2,3,1)
% plot(2000:1999+num_year,PX_tt(1,1:num_year)./GDP_tt(1,1:num_year),'-b','Linewidth',2.5)
% title('Brazil')
% grid on
% hold on
% plot(2000:1999+num_year,PX_t(1,1:num_year)./GDP_t(1,1:num_year),'--b','Linewidth',2.5)
% hold off
% 
% subplot(2,3,2)
% plot(2000:1999+num_year,PX_tt(2,1:num_year)./GDP_tt(2,1:num_year),'-r','Linewidth',2.5)
% title('China')
% grid on
% hold on
% plot(2000:1999+num_year,PX_t(2,1:num_year)./GDP_t(2,1:num_year),'--r','Linewidth',2.5)
% hold off
% 
% subplot(2,3,3)
% plot(2000:1999+num_year,PX_tt(3,1:num_year)./GDP_tt(3,1:num_year),'-g','Linewidth',2.5)
% title('USA')
% grid on
% hold on
% plot(2000:1999+num_year,PX_t(3,1:num_year)./GDP_t(3,1:num_year),'--g','Linewidth',2.5)
% hold off
% 
% subplot(2,3,4)
% plot(2000:1999+num_year,PX_tt(4,1:num_year)./GDP_tt(4,1:num_year),'-c','Linewidth',2.5)
% title('IND')
% grid on
% hold on
% plot(2000:1999+num_year,PX_t(4,1:num_year)./GDP_t(4,1:num_year),'--c','Linewidth',2.5)
% hold off
% 
% subplot(2,3,5)
% plot(2000:1999+num_year,PX_tt(5,1:num_year)./GDP_tt(5,1:num_year),'-m','Linewidth',2.5)
% title('ROW')
% grid on
% hold on
% plot(2000:1999+num_year,PX_t(5,1:num_year)./GDP_t(5,1:num_year),'--m','Linewidth',2.5)
% hold off
% 
% %sgtitle('Figure 7. Implied capital stock by country 2000-2014')
% sgtitle('PX over GDP by country 2000-2014 PWT')
% saveas(gcf, 'plots_transition\final_2\Figure 3c.png');
% close(gcf)
% 
% %ddd
% 
% % untargeted sector value share (over GDP)
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% % plot(2000:1999+ye,va_a_tt(1,1:ye),'-b',2000:1999+ye,va_a_tt(2,1:ye),'-r',2000:1999+ye,va_a_tt(3,1:ye),'-g',2000:1999+ye,va_a_tt(4,1:ye),'-m',...
% %     2000:2000+num_year-1,va_a_t,'--black','Linewidth',2.5)
% plot(2000:1999+ye,va_g_tt(1,1:ye),'-b',2000:1999+ye,va_g_tt(2,1:ye),'-r',2000:1999+ye,va_g_tt(3,1:ye),'-g',2000:1999+ye,va_g_tt(4,1:ye),'-m',...
%      2000:1999+ye,va_g_t(1,1:ye),'--b',2000:1999+ye,va_g_t(2,1:ye),'--r',2000:1999+ye,va_g_t(3,1:ye),'--g',2000:1999+ye,va_g_t(4,1:ye),'--m','Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','ROW','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,va_l_tt(1,1:ye),'-b',2000:1999+ye,va_l_tt(2,1:ye),'-r',2000:1999+ye,va_l_tt(3,1:ye),'-g',2000:1999+ye,va_l_tt(4,1:ye),'-m',...
%      2000:1999+ye,va_l_t(1,1:ye),'--b',2000:1999+ye,va_l_t(2,1:ye),'--r',2000:1999+ye,va_l_t(3,1:ye),'--g',2000:1999+ye,va_l_t(4,1:ye),'--m','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,va_h_tt(1,1:ye),'-b',2000:1999+ye,va_h_tt(2,1:ye),'-r',2000:1999+ye,va_h_tt(3,1:ye),'-g',2000:1999+ye,va_h_tt(4,1:ye),'-m',...
%      2000:1999+ye,va_h_t(1,1:ye),'--b',2000:1999+ye,va_h_t(2,1:ye),'--r',2000:1999+ye,va_h_t(3,1:ye),'--g',2000:1999+ye,va_h_t(4,1:ye),'--m','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,va_s_tt(1,1:ye),'-b',2000:1999+ye,va_s_tt(2,1:ye),'-r',2000:1999+ye,va_s_tt(3,1:ye),'-g',2000:1999+ye,va_s_tt(4,1:ye),'-m',...
%      2000:1999+ye,va_s_t(1,1:ye),'--b',2000:1999+ye,va_s_t(2,1:ye),'--r',2000:1999+ye,va_s_t(3,1:ye),'--g',2000:1999+ye,va_s_t(4,1:ye),'--m','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure 8. Implied sector value added share by sector and country 2000-2014')
% %saveas(gcf, 'plots_transition\final_2\Figure 8.png');
% close(gcf)
% 
% % untargeted sector capital share 
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,K_g_tt(1,1:ye)./K_tt(1,1:ye),'-b',2000:1999+ye,K_g_tt(2,1:ye)./K_tt(2,1:ye),'-r',2000:1999+ye,K_g_tt(3,1:ye)./K_tt(3,1:ye),'-g',2000:1999+ye,K_g_tt(4,1:ye)./K_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,kh_g_t,'--black','Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','ROW','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,K_l_tt(1,1:ye)./K_tt(1,1:ye),'-b',2000:1999+ye,K_l_tt(2,1:ye)./K_tt(2,1:ye),'-r',2000:1999+ye,K_l_tt(3,1:ye)./K_tt(3,1:ye),'-g',2000:1999+ye,K_l_tt(4,1:ye)./K_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,kh_l_t,'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,K_h_tt(1,1:ye)./K_tt(1,1:ye),'-b',2000:1999+ye,K_h_tt(2,1:ye)./K_tt(2,1:ye),'-r',2000:1999+ye,K_h_tt(3,1:ye)./K_tt(3,1:ye),'-g',2000:1999+ye,K_h_tt(4,1:ye)./K_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,kh_h_t,'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,K_s_tt(1,1:ye)./K_tt(1,1:ye),'-b',2000:1999+ye,K_s_tt(2,1:ye)./K_tt(2,1:ye),'-r',2000:1999+ye,K_s_tt(3,1:ye)./K_tt(3,1:ye),'-g',2000:1999+ye,K_s_tt(4,1:ye)./K_tt(4,1:ye),'-m',...
%     2000:2000+num_year-1,kh_s_t,'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Untargeted sector capital share 2000-2014')
% saveas(gcf, 'plots_transition/Untargeted sector capital share 2000-2014.png');
% close(gcf)
% 
% 
% 
% % Fitting bilateral trade share on Brazil for other countries 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(pi_g_tt(2,1,1:ye)),'-r',2000:1999+ye,squeeze(pi_g_tt(3,1,1:ye)),'-g',2000:1999+ye,squeeze(pi_g_tt(4,1,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_g_t(2:4,1,:)),'--black','Linewidth',2.5)
% title('Primary')
% legend('CHN','USA','ROW','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(pi_l_tt(2,1,1:ye)),'-r',2000:1999+ye,squeeze(pi_l_tt(3,1,1:ye)),'-g',2000:1999+ye,squeeze(pi_l_tt(4,1,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_l_t(2:4,1,:)),'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(pi_h_tt(2,1,1:ye)),'-r',2000:1999+ye,squeeze(pi_h_tt(3,1,1:ye)),'-g',2000:1999+ye,squeeze(pi_h_tt(4,1,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_h_t(2:4,1,:)),'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(pi_s_tt(2,1,1:ye)),'-r',2000:1999+ye,squeeze(pi_s_tt(3,1,1:ye)),'-g',2000:1999+ye,squeeze(pi_s_tt(4,1,1:ye)),'-m',...
%     2000:2000+num_year-1,squeeze(bt_s_t(2:4,1,:)),'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Fitting bilateral trade share on Brazil goods 2000-2014')
% saveas(gcf, 'plots_transition/Fitting bilateral trade share on Brazil goods 2000-2014.png');
% close(gcf)
% 
% % Fitting bilateral import share on Brazil for other countries 2000-2014
% figure('Position', get(0, 'Screensize'))
% subplot(2,2,1)
% plot(2000:1999+ye,squeeze(pi_g_tt(2,1,1:ye))./(1-squeeze(pi_g_tt(2,2,1:ye))),'-r',2000:1999+ye,squeeze(pi_g_tt(3,1,1:ye))./(1-squeeze(pi_g_tt(3,3,1:ye))),'-g',2000:1999+ye,squeeze(pi_g_tt(4,1,1:ye))./(1-squeeze(pi_g_tt(4,4,1:ye))),'-m',...
%      2000:1999+ye,squeeze(bt_g_t(2,1,1:ye))./(1-squeeze(bt_g_t(2,2,1:ye))),'--black',2000:1999+ye,squeeze(bt_g_t(3,1,1:ye))./(1-squeeze(bt_g_t(3,3,1:ye))),'--black',2000:1999+ye,squeeze(bt_g_t(4,1,1:ye))./(1-squeeze(bt_g_t(4,4,1:ye))),'--black','Linewidth',2.5)
% title('Primary')
% legend('CHN','USA','ROW','Target','location','best')
% grid on
% 
% subplot(2,2,2)
% plot(2000:1999+ye,squeeze(pi_l_tt(2,1,1:ye))./(1-squeeze(pi_l_tt(2,2,1:ye))),'-r',2000:1999+ye,squeeze(pi_l_tt(3,1,1:ye))./(1-squeeze(pi_l_tt(3,3,1:ye))),'-g',2000:1999+ye,squeeze(pi_l_tt(4,1,1:ye))./(1-squeeze(pi_l_tt(4,4,1:ye))),'-m',...
%      2000:1999+ye,squeeze(bt_l_t(2,1,1:ye))./(1-squeeze(bt_l_t(2,2,1:ye))),'--black',2000:1999+ye,squeeze(bt_l_t(3,1,1:ye))./(1-squeeze(bt_l_t(3,3,1:ye))),'--black',2000:1999+ye,squeeze(bt_l_t(4,1,1:ye))./(1-squeeze(bt_l_t(4,4,1:ye))),'--black','Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(2,2,3)
% plot(2000:1999+ye,squeeze(pi_h_tt(2,1,1:ye))./(1-squeeze(pi_h_tt(2,2,1:ye))),'-r',2000:1999+ye,squeeze(pi_h_tt(3,1,1:ye))./(1-squeeze(pi_h_tt(3,3,1:ye))),'-g',2000:1999+ye,squeeze(pi_h_tt(4,1,1:ye))./(1-squeeze(pi_h_tt(4,4,1:ye))),'-m',...
%      2000:1999+ye,squeeze(bt_h_t(2,1,1:ye))./(1-squeeze(bt_h_t(2,2,1:ye))),'--black',2000:1999+ye,squeeze(bt_h_t(3,1,1:ye))./(1-squeeze(bt_h_t(3,3,1:ye))),'--black',2000:1999+ye,squeeze(bt_h_t(4,1,1:ye))./(1-squeeze(bt_h_t(4,4,1:ye))),'--black','Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(2,2,4)
% plot(2000:1999+ye,squeeze(pi_s_tt(2,1,1:ye))./(1-squeeze(pi_s_tt(2,2,1:ye))),'-r',2000:1999+ye,squeeze(pi_s_tt(3,1,1:ye))./(1-squeeze(pi_s_tt(3,3,1:ye))),'-g',2000:1999+ye,squeeze(pi_s_tt(4,1,1:ye))./(1-squeeze(pi_s_tt(4,4,1:ye))),'-m',...
%      2000:1999+ye,squeeze(bt_s_t(2,1,1:ye))./(1-squeeze(bt_s_t(2,2,1:ye))),'--black',2000:1999+ye,squeeze(bt_s_t(3,1,1:ye))./(1-squeeze(bt_s_t(3,3,1:ye))),'--black',2000:1999+ye,squeeze(bt_s_t(4,1,1:ye))./(1-squeeze(bt_s_t(4,4,1:ye))),'--black','Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Fitting bilateral import share on Brazil goods 2000-2014')
% saveas(gcf, 'plots_transition/Fitting bilateral import share on Brazil goods.png');
% close(gcf)


%% Decompose trade share
% Decompose trade share on Brazil for other countries 2000-2014
id = 1;
[pi_g_BM,pi_g_T,pi_g_u,pi_g_tao] = tradeshare_decompose(I,id,theta,T_g_tt(:,1:ye),u_g_tt(:,1:ye),d_g_tt(:,:,1:ye));
[pi_l_BM,pi_l_T,pi_l_u,pi_l_tao] = tradeshare_decompose(I,id,theta,T_l_tt(:,1:ye),u_l_tt(:,1:ye),d_l_tt(:,:,1:ye));
[pi_h_BM,pi_h_T,pi_h_u,pi_h_tao] = tradeshare_decompose(I,id,theta,T_h_tt(:,1:ye),u_h_tt(:,1:ye),d_h_tt(:,:,1:ye));

% 
% figure('Position', get(0, 'Screensize'))
% subplot(4,4,1)
% plot(2000:1999+ye,pi_g_BM(2,:),'-r',2000:1999+ye,pi_g_BM(3,:),'-g',2000:1999+ye,pi_g_BM(4,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('CHN','USA','ROW','Target','location','best')
% grid on
% 
% subplot(4,4,2)
% plot(2000:1999+ye,pi_m_BM(2,:),'-r',2000:1999+ye,pi_m_BM(3,:),'-g',2000:1999+ye,pi_m_BM(4,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,3)
% plot(2000:1999+ye,pi_h_BM(2,:),'-r',2000:1999+ye,pi_h_BM(3,:),'-g',2000:1999+ye,pi_h_BM(4,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,4)
% plot(2000:1999+ye,pi_s_BM(2,:),'-r',2000:1999+ye,pi_s_BM(3,:),'-g',2000:1999+ye,pi_s_BM(4,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% subplot(4,4,5)
% plot(2000:1999+ye,pi_a_T(2,:),'-r',2000:1999+ye,pi_a_T(3,:),'-g',2000:1999+ye,pi_a_T(4,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('CHN_T','USA_T','ROW_T','Target','location','best')
% grid on
% 
% subplot(4,4,6)
% plot(2000:1999+ye,pi_m_T(2,:),'-r',2000:1999+ye,pi_m_T(3,:),'-g',2000:1999+ye,pi_m_T(4,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,7)
% plot(2000:1999+ye,pi_h_T(2,:),'-r',2000:1999+ye,pi_h_T(3,:),'-g',2000:1999+ye,pi_h_T(4,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,8)
% plot(2000:1999+ye,pi_s_T(2,:),'-r',2000:1999+ye,pi_s_T(3,:),'-g',2000:1999+ye,pi_s_T(4,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% subplot(4,4,9)
% plot(2000:1999+ye,pi_a_u(2,:),'-r',2000:1999+ye,pi_a_u(3,:),'-g',2000:1999+ye,pi_a_u(4,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('CHN_u','USA_u','ROW_u','Target','location','best')
% grid on
% 
% subplot(4,4,10)
% plot(2000:1999+ye,pi_m_u(2,:),'-r',2000:1999+ye,pi_m_u(3,:),'-g',2000:1999+ye,pi_m_u(4,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,11)
% plot(2000:1999+ye,pi_h_u(2,:),'-r',2000:1999+ye,pi_h_u(3,:),'-g',2000:1999+ye,pi_h_u(4,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,12)
% plot(2000:1999+ye,pi_s_u(2,:),'-r',2000:1999+ye,pi_s_u(3,:),'-g',2000:1999+ye,pi_s_u(4,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% subplot(4,4,13)
% plot(2000:1999+ye,pi_a_tao(2,:),'-r',2000:1999+ye,pi_a_tao(3,:),'-g',2000:1999+ye,pi_a_tao(4,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('CHN_t_a_o','USA_t_a_o','ROW_t_a_o','Target','location','best')
% grid on
% 
% subplot(4,4,14)
% plot(2000:1999+ye,pi_m_tao(2,:),'-r',2000:1999+ye,pi_m_tao(3,:),'-g',2000:1999+ye,pi_m_tao(4,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,15)
% plot(2000:1999+ye,pi_h_tao(2,:),'-r',2000:1999+ye,pi_h_tao(3,:),'-g',2000:1999+ye,pi_h_tao(4,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,16)
% plot(2000:1999+ye,pi_s_tao(2,:),'-r',2000:1999+ye,pi_s_tao(3,:),'-g',2000:1999+ye,pi_s_tao(4,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Decompose trade share on Brazil goods 2000-2014')
% saveas(gcf, 'plots_transition/Decompose trade share on Brazil goods 2000-2014.png');
% close(gcf)
% 
% 
% % Decompose trade share on China for other countries 2000-2014
id = 2;
[pi_g_BM,pi_g_T,pi_g_u,pi_g_tao] = tradeshare_decompose(I,id,theta,T_g_tt(:,1:ye),u_g_tt(:,1:ye),d_g_tt(:,:,1:ye));
[pi_l_BM,pi_l_T,pi_l_u,pi_l_tao] = tradeshare_decompose(I,id,theta,T_l_tt(:,1:ye),u_l_tt(:,1:ye),d_l_tt(:,:,1:ye));
[pi_h_BM,pi_h_T,pi_h_u,pi_h_tao] = tradeshare_decompose(I,id,theta,T_h_tt(:,1:ye),u_h_tt(:,1:ye),d_h_tt(:,:,1:ye));

% 
% figure('Position', get(0, 'Screensize'))
% subplot(4,4,1)
% plot(2000:1999+ye,pi_g_BM(1,:),'-b',2000:1999+ye,pi_g_BM(3,:),'-g',2000:1999+ye,pi_g_BM(4,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA','USA','ROW','Target','location','best')
% grid on
% 
% subplot(4,4,2)
% plot(2000:1999+ye,pi_m_BM(1,:),'-b',2000:1999+ye,pi_m_BM(3,:),'-g',2000:1999+ye,pi_m_BM(4,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,3)
% plot(2000:1999+ye,pi_h_BM(1,:),'-b',2000:1999+ye,pi_h_BM(3,:),'-g',2000:1999+ye,pi_h_BM(4,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,4)
% plot(2000:1999+ye,pi_s_BM(1,:),'-b',2000:1999+ye,pi_s_BM(3,:),'-g',2000:1999+ye,pi_s_BM(4,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% subplot(4,4,5)
% plot(2000:1999+ye,pi_a_T(1,:),'-b',2000:1999+ye,pi_a_T(3,:),'-g',2000:1999+ye,pi_a_T(4,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA_T','USA_T','ROW_T','Target','location','best')
% grid on
% 
% subplot(4,4,6)
% plot(2000:1999+ye,pi_m_T(1,:),'-b',2000:1999+ye,pi_m_T(3,:),'-g',2000:1999+ye,pi_m_T(4,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,7)
% plot(2000:1999+ye,pi_h_T(1,:),'-b',2000:1999+ye,pi_h_T(3,:),'-g',2000:1999+ye,pi_h_T(4,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,8)
% plot(2000:1999+ye,pi_s_T(1,:),'-b',2000:1999+ye,pi_s_T(3,:),'-g',2000:1999+ye,pi_s_T(4,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% subplot(4,4,9)
% plot(2000:1999+ye,pi_a_u(1,:),'-b',2000:1999+ye,pi_a_u(3,:),'-g',2000:1999+ye,pi_a_u(4,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA_u','USA_u','ROW_u','Target','location','best')
% grid on
% 
% subplot(4,4,10)
% plot(2000:1999+ye,pi_m_u(1,:),'-b',2000:1999+ye,pi_m_u(3,:),'-g',2000:1999+ye,pi_m_u(4,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,11)
% plot(2000:1999+ye,pi_h_u(1,:),'-b',2000:1999+ye,pi_h_u(3,:),'-g',2000:1999+ye,pi_h_u(4,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,12)
% plot(2000:1999+ye,pi_s_u(1,:),'-b',2000:1999+ye,pi_s_u(3,:),'-g',2000:1999+ye,pi_s_u(4,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% subplot(4,4,13)
% plot(2000:1999+ye,pi_a_tao(1,:),'-b',2000:1999+ye,pi_a_tao(3,:),'-g',2000:1999+ye,pi_a_tao(4,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA_t_a_o','USA_t_a_o','ROW_t_a_o','Target','location','best')
% grid on
% 
% subplot(4,4,14)
% plot(2000:1999+ye,pi_m_tao(1,:),'-b',2000:1999+ye,pi_m_tao(3,:),'-g',2000:1999+ye,pi_m_tao(4,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,15)
% plot(2000:1999+ye,pi_h_tao(1,:),'-b',2000:1999+ye,pi_h_tao(3,:),'-g',2000:1999+ye,pi_h_tao(4,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,16)
% plot(2000:1999+ye,pi_s_tao(1,:),'-b',2000:1999+ye,pi_s_tao(3,:),'-g',2000:1999+ye,pi_s_tao(4,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Decompose trade share on China goods 2000-2014')
% saveas(gcf, 'plots_transition/Decompose trade share on China goods 2000-2014.png');
% close(gcf)
% 
% 
% 
% % Decompose trade share on USA for other countries 2000-2014
id = 9;
[pi_g_BM,pi_g_T,pi_g_u,pi_g_tao] = tradeshare_decompose(I,id,theta,T_g_tt(:,1:ye),u_g_tt(:,1:ye),d_g_tt(:,:,1:ye));
[pi_l_BM,pi_l_T,pi_l_u,pi_l_tao] = tradeshare_decompose(I,id,theta,T_l_tt(:,1:ye),u_l_tt(:,1:ye),d_l_tt(:,:,1:ye));
[pi_h_BM,pi_h_T,pi_h_u,pi_h_tao] = tradeshare_decompose(I,id,theta,T_h_tt(:,1:ye),u_h_tt(:,1:ye),d_h_tt(:,:,1:ye));

% 
% figure('Position', get(0, 'Screensize'))
% subplot(4,4,1)
% plot(2000:1999+ye,pi_g_BM(1,:),'-b',2000:1999+ye,pi_g_BM(2,:),'-r',2000:1999+ye,pi_g_BM(4,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','ROW','Target','location','best')
% grid on
% 
% subplot(4,4,2)
% plot(2000:1999+ye,pi_m_BM(1,:),'-b',2000:1999+ye,pi_m_BM(2,:),'-r',2000:1999+ye,pi_m_BM(4,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,3)
% plot(2000:1999+ye,pi_h_BM(1,:),'-b',2000:1999+ye,pi_h_BM(2,:),'-r',2000:1999+ye,pi_h_BM(4,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,4)
% plot(2000:1999+ye,pi_s_BM(1,:),'-b',2000:1999+ye,pi_s_BM(2,:),'-r',2000:1999+ye,pi_s_BM(4,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% subplot(4,4,5)
% plot(2000:1999+ye,pi_a_T(1,:),'-b',2000:1999+ye,pi_a_T(2,:),'-r',2000:1999+ye,pi_a_T(4,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA_T','CHN_T','ROW_T','Target','location','best')
% grid on
% 
% subplot(4,4,6)
% plot(2000:1999+ye,pi_m_T(1,:),'-b',2000:1999+ye,pi_m_T(2,:),'-r',2000:1999+ye,pi_m_T(4,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,7)
% plot(2000:1999+ye,pi_h_T(1,:),'-b',2000:1999+ye,pi_h_T(2,:),'-r',2000:1999+ye,pi_h_T(4,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,8)
% plot(2000:1999+ye,pi_s_T(1,:),'-b',2000:1999+ye,pi_s_T(2,:),'-r',2000:1999+ye,pi_s_T(4,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% subplot(4,4,9)
% plot(2000:1999+ye,pi_a_u(1,:),'-b',2000:1999+ye,pi_a_u(2,:),'-r',2000:1999+ye,pi_a_u(4,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA_u','CHN_u','ROW_u','Target','location','best')
% grid on
% 
% subplot(4,4,10)
% plot(2000:1999+ye,pi_m_u(1,:),'-b',2000:1999+ye,pi_m_u(2,:),'-r',2000:1999+ye,pi_m_u(4,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,11)
% plot(2000:1999+ye,pi_h_u(1,:),'-b',2000:1999+ye,pi_h_u(2,:),'-r',2000:1999+ye,pi_h_u(4,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,12)
% plot(2000:1999+ye,pi_s_u(1,:),'-b',2000:1999+ye,pi_s_u(2,:),'-r',2000:1999+ye,pi_s_u(4,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% subplot(4,4,13)
% plot(2000:1999+ye,pi_a_tao(1,:),'-b',2000:1999+ye,pi_a_tao(2,:),'-r',2000:1999+ye,pi_a_tao(4,:),'-m',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA_t_a_o','CHN_t_a_o','ROW_t_a_o','Target','location','best')
% grid on
% 
% subplot(4,4,14)
% plot(2000:1999+ye,pi_m_tao(1,:),'-b',2000:1999+ye,pi_m_tao(2,:),'-r',2000:1999+ye,pi_m_tao(4,:),'-m',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,15)
% plot(2000:1999+ye,pi_h_tao(1,:),'-b',2000:1999+ye,pi_h_tao(2,:),'-r',2000:1999+ye,pi_h_tao(4,:),'-m',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,16)
% plot(2000:1999+ye,pi_s_tao(1,:),'-b',2000:1999+ye,pi_s_tao(2,:),'-r',2000:1999+ye,pi_s_tao(4,:),'-m',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Decompose trade share on USA goods 2000-2014')
% saveas(gcf, 'plots_transition/Decompose trade share on USA goods 2000-2014.png');
% close(gcf)



% Decompose trade share on ROW for other countries 2000-2014
id = 8;
[pi_g_BM,pi_g_T,pi_g_u,pi_g_tao] = tradeshare_decompose(I,id,theta,T_g_tt(:,1:ye),u_g_tt(:,1:ye),d_g_tt(:,:,1:ye));
[pi_l_BM,pi_l_T,pi_l_u,pi_l_tao] = tradeshare_decompose(I,id,theta,T_l_tt(:,1:ye),u_l_tt(:,1:ye),d_l_tt(:,:,1:ye));
[pi_h_BM,pi_h_T,pi_h_u,pi_h_tao] = tradeshare_decompose(I,id,theta,T_h_tt(:,1:ye),u_h_tt(:,1:ye),d_h_tt(:,:,1:ye));

% 
% figure('Position', get(0, 'Screensize'))
% subplot(4,4,1)
% plot(2000:1999+ye,pi_a_BM(1,:),'-b',2000:1999+ye,pi_a_BM(2,:),'-r',2000:1999+ye,pi_a_BM(3,:),'-g',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA','CHN','USA','Target','location','best')
% grid on
% 
% subplot(4,4,2)
% plot(2000:1999+ye,pi_m_BM(1,:),'-b',2000:1999+ye,pi_m_BM(2,:),'-r',2000:1999+ye,pi_m_BM(3,:),'-g',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,3)
% plot(2000:1999+ye,pi_h_BM(1,:),'-b',2000:1999+ye,pi_h_BM(2,:),'-r',2000:1999+ye,pi_h_BM(3,:),'-g',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,4)
% plot(2000:1999+ye,pi_s_BM(1,:),'-b',2000:1999+ye,pi_s_BM(2,:),'-r',2000:1999+ye,pi_s_BM(3,:),'-g',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% subplot(4,4,5)
% plot(2000:1999+ye,pi_a_T(1,:),'-b',2000:1999+ye,pi_a_T(2,:),'-r',2000:1999+ye,pi_a_T(3,:),'-g',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA_T','CHN_T','USA_T','Target','location','best')
% grid on
% 
% subplot(4,4,6)
% plot(2000:1999+ye,pi_m_T(1,:),'-b',2000:1999+ye,pi_m_T(2,:),'-r',2000:1999+ye,pi_m_T(3,:),'-g',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,7)
% plot(2000:1999+ye,pi_h_T(1,:),'-b',2000:1999+ye,pi_h_T(2,:),'-r',2000:1999+ye,pi_h_T(3,:),'-g',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,8)
% plot(2000:1999+ye,pi_s_T(1,:),'-b',2000:1999+ye,pi_s_T(2,:),'-r',2000:1999+ye,pi_s_T(3,:),'-g',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% subplot(4,4,9)
% plot(2000:1999+ye,pi_a_u(1,:),'-b',2000:1999+ye,pi_a_u(2,:),'-r',2000:1999+ye,pi_a_u(3,:),'-g',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA_u','CHN_u','USA_u','Target','location','best')
% grid on
% 
% subplot(4,4,10)
% plot(2000:1999+ye,pi_m_u(1,:),'-b',2000:1999+ye,pi_m_u(2,:),'-r',2000:1999+ye,pi_m_u(3,:),'-g',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,11)
% plot(2000:1999+ye,pi_h_u(1,:),'-b',2000:1999+ye,pi_h_u(2,:),'-r',2000:1999+ye,pi_h_u(3,:),'-g',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,12)
% plot(2000:1999+ye,pi_s_u(1,:),'-b',2000:1999+ye,pi_s_u(2,:),'-r',2000:1999+ye,pi_s_u(3,:),'-g',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% subplot(4,4,13)
% plot(2000:1999+ye,pi_a_tao(1,:),'-b',2000:1999+ye,pi_a_tao(2,:),'-r',2000:1999+ye,pi_a_tao(3,:),'-g',...
%     'Linewidth',2.5)
% title('Primary')
% legend('BRA_t_a_o','CHN_t_a_o','USA_t_a_o','Target','location','best')
% grid on
% 
% subplot(4,4,14)
% plot(2000:1999+ye,pi_m_tao(1,:),'-b',2000:1999+ye,pi_m_tao(2,:),'-r',2000:1999+ye,pi_m_tao(3,:),'-g',...
%     'Linewidth',2.5)
% title('Low-tech manufacturing')
% grid on
% 
% subplot(4,4,15)
% plot(2000:1999+ye,pi_h_tao(1,:),'-b',2000:1999+ye,pi_h_tao(2,:),'-r',2000:1999+ye,pi_h_tao(3,:),'-g',...
%     'Linewidth',2.5)
% title('High-tech manufacturing')
% grid on
% 
% subplot(4,4,16)
% plot(2000:1999+ye,pi_s_tao(1,:),'-b',2000:1999+ye,pi_s_tao(2,:),'-r',2000:1999+ye,pi_s_tao(3,:),'-g',...
%     'Linewidth',2.5)
% title('Service')
% grid on
% 
% sgtitle('Figure. Decompose trade share on ROW goods 2000-2014')
% saveas(gcf, 'plots_transition/Decompose trade share on ROW goods 2000-2014.png');
% close(gcf)




%% Save calibrated shocks

mkdir plots_transition\parameters

% parameters
diary(['plots_transition\parameters\parameters.txt']);
disp('I');
[I]
disp('theta');
[theta]
disp('eta');
[eta]
disp('del');
[del]
disp('lambda');
[lambda]
disp('beta');
beta
disp('transition length');
TT

diary off;

K_t_raw = k_t;

%% shocks
save('plots_transition\parameters\K_1.txt','K_1','-ASCII');
save('plots_transition\parameters\L_t.txt','l_t','-ASCII');
save('plots_transition\parameters\K_t_raw.txt','K_t_raw','-ASCII');

save('plots_transition\parameters\omega_bar.txt','omega_bar','-ASCII');
save('plots_transition\parameters\fi_t.txt','fi_t','-ASCII');
save('plots_transition\parameters\omega_cg_t.txt','omega_cg_t','-ASCII');
save('plots_transition\parameters\omega_cl_t.txt','omega_cl_t','-ASCII');
save('plots_transition\parameters\omega_ch_t.txt','omega_ch_t','-ASCII');
save('plots_transition\parameters\omega_xg_t.txt','omega_xg_t','-ASCII');
save('plots_transition\parameters\omega_xl_t.txt','omega_xl_t','-ASCII');
save('plots_transition\parameters\omega_xh_t.txt','omega_xh_t','-ASCII');

save('plots_transition\parameters\Ax_t.txt','Ax_t','-ASCII');
save('plots_transition\parameters\T_g_t.txt','T_g_t','-ASCII');
save('plots_transition\parameters\T_l_t.txt','T_l_t','-ASCII');
save('plots_transition\parameters\T_h_t.txt','T_h_t','-ASCII');

save('plots_transition\parameters\alpha_g_t.txt','alpha_g_t','-ASCII');
save('plots_transition\parameters\alpha_m_t.txt','alpha_l_t','-ASCII');
save('plots_transition\parameters\alpha_h_t.txt','alpha_h_t','-ASCII');
save('plots_transition\parameters\nu_g_t.txt','nu_g_t','-ASCII');
save('plots_transition\parameters\nu_m_t.txt','nu_l_t','-ASCII');
save('plots_transition\parameters\nu_h_t.txt','nu_h_t','-ASCII');

save('plots_transition\parameters\mu_gg_t.txt','mu_gg_t','-ASCII');
save('plots_transition\parameters\mu_gl_t.txt','mu_gl_t','-ASCII');
save('plots_transition\parameters\mu_gh_t.txt','mu_gh_t','-ASCII');
save('plots_transition\parameters\mu_lg_t.txt','mu_lg_t','-ASCII');
save('plots_transition\parameters\mu_ll_t.txt','mu_ll_t','-ASCII');
save('plots_transition\parameters\mu_lh_t.txt','mu_lh_t','-ASCII');
save('plots_transition\parameters\mu_hg_t.txt','mu_hg_t','-ASCII');
save('plots_transition\parameters\mu_hl_t.txt','mu_hl_t','-ASCII');
save('plots_transition\parameters\mu_hh_t.txt','mu_hh_t','-ASCII');

trade_cost = {tao_g_t,tao_l_t,tao_h_t};
save('plots_transition\parameters\trade_cost.mat','trade_cost');

% save initial guess based from data
save('plots_transition\parameters\w_t.txt','w_t','-ASCII');
save('plots_transition\parameters\lambda.txt','lambda','-ASCII');
save('plots_transition\parameters\del.txt','del','-ASCII');
save('plots_transition\parameters\TT.txt','TT','-ASCII');


%% fig 1
years = 1995:2011;
T = length(years);

p25 = zeros(T,3);
p50 = zeros(T,3);
p75 = zeros(T,3);

for t = 1:T
    % Goods
    p25(t,1) = prctile(tao_g_t(:,t)./tao_g_t(:,1),25);
    p50(t,1) = prctile(tao_g_t(:,t)./tao_g_t(:,1),50);
    p75(t,1) = prctile(tao_g_t(:,t)./tao_g_t(:,1),75);

    % Low-skill services
    p25(t,2) = prctile(tao_l_t(:,t)./tao_l_t(:,1),25);
    p50(t,2) = prctile(tao_l_t(:,t)./tao_l_t(:,1),50);
    p75(t,2) = prctile(tao_l_t(:,t)./tao_l_t(:,1),75);

    % High-skill services
    p25(t,3) = prctile(tao_h_t(:,t)./tao_h_t(:,1),25);
    p50(t,3) = prctile(tao_h_t(:,t)./tao_h_t(:,1),50);
    p75(t,3) = prctile(tao_h_t(:,t)./tao_h_t(:,1),75);
end

sector_names = {'Goods','Low-skill services','High-skill services'};

figure('Position',[100 100 1000 350]);

for s = 1:3
    subplot(1,3,s); hold on

    x = years;

    % ---- Shaded 25–75 band ----
    fill([x fliplr(x)], ...
         [p25(:,s)' fliplr(p75(:,s)')], ...
         [0.8 0.85 1], ...
         'EdgeColor','none','FaceAlpha',0.8);

    % ---- Percentile lines ----
    plot(x, p25(:,s), 'b--', 'LineWidth',1.2);
    plot(x, p75(:,s), 'b--', 'LineWidth',1.2);

    % ---- Median ----
    plot(x, p50(:,s), 'b-', 'LineWidth',2);

    title(sector_names{s});
    xlabel('Year');
    if s==1
        ylabel('Trade cost');
    end
    xlim([min(x) max(x)]);
    grid on
end

sgtitle({'Distribution of Bilateral Trade Costs by Sector and Year', ...
         '(25th percentile, median and 75th percentile)'});

saveas(gcf, 'trade_cost_distribution_3sectors.png');


%% fig 1a
years = 1995:2011;
T = length(years);

% squeeze(d_h_tt(2,5,1:ye))./d_h_tt(2,5,1)

p25 = zeros(T,6);
p50 = zeros(T,6);
p75 = zeros(T,6);

for t = 1:T
    % Goods Export
    p25(t,1) = prctile(squeeze(d_g_tt(:,5,t))./d_g_tt(:,5,1),25);
    p50(t,1) = prctile(squeeze(d_g_tt(:,5,t))./d_g_tt(:,5,1),50);
    p75(t,1) = prctile(squeeze(d_g_tt(:,5,t))./d_g_tt(:,5,1),75);

    % Low-skill services
    p25(t,2) = prctile(squeeze(d_l_tt(:,5,t))./d_l_tt(:,5,1),25);
    p50(t,2) = prctile(squeeze(d_l_tt(:,5,t))./d_l_tt(:,5,1),50);
    p75(t,2) = prctile(squeeze(d_l_tt(:,5,t))./d_l_tt(:,5,1),75);

    % High-skill services
    p25(t,3) = prctile(squeeze(d_h_tt(:,5,t))./d_h_tt(:,5,1),25);
    p50(t,3) = prctile(squeeze(d_h_tt(:,5,t))./d_h_tt(:,5,1),50);
    p75(t,3) = prctile(squeeze(d_h_tt(:,5,t))./d_h_tt(:,5,1),75);

    % Goods Import
    p25(t,4) = prctile(squeeze(d_g_tt(5,:,t))./d_g_tt(5,:,1),25);
    p50(t,4) = prctile(squeeze(d_g_tt(5,:,t))./d_g_tt(5,:,1),50);
    p75(t,4) = prctile(squeeze(d_g_tt(5,:,t))./d_g_tt(5,:,1),75);

    p25(t,5) = prctile(squeeze(d_l_tt(5,:,t))./d_l_tt(5,:,1),25);
    p50(t,5) = prctile(squeeze(d_l_tt(5,:,t))./d_l_tt(5,:,1),50);
    p75(t,5) = prctile(squeeze(d_l_tt(5,:,t))./d_l_tt(5,:,1),75);

    p25(t,6) = prctile(squeeze(d_h_tt(5,:,t))./d_h_tt(5,:,1),25);
    p50(t,6) = prctile(squeeze(d_h_tt(5,:,t))./d_h_tt(5,:,1),50);
    p75(t,6) = prctile(squeeze(d_h_tt(5,:,t))./d_h_tt(5,:,1),75);
end

sector_names = {'Goods export','Low-skill services export','High-skill services export','Goods import','Low-skill services import','High-skill services import'};

figure('Position',[100 100 1000 350]);

for s = 1:6
    subplot(2,3,s); hold on

    x = years;

    % ---- Shaded 25–75 band ----
    fill([x fliplr(x)], ...
         [p25(:,s)' fliplr(p75(:,s)')], ...
         [0.8 0.85 1], ...
         'EdgeColor','none','FaceAlpha',0.8);

    % ---- Percentile lines ----
    plot(x, p25(:,s), 'b--', 'LineWidth',1.2);
    plot(x, p75(:,s), 'b--', 'LineWidth',1.2);

    % ---- Median ----
    plot(x, p50(:,s), 'b-', 'LineWidth',2);

    title(sector_names{s});
    xlabel('Year');
    if s==1
        ylabel('Trade cost');
    end
    xlim([min(x) max(x)]);
    grid on
end

sgtitle({'Distribution of Bilateral Trade Costs for India by Sector and Year', ...
         '(25th percentile, median and 75th percentile)'});

saveas(gcf, 'trade_cost_distribution_india.png');



%% fig 2

years = 1995:2011;                 % adjust if needed
T = length(years);

p25 = zeros(T,1);
p50 = zeros(T,1);
p75 = zeros(T,1);

for t = 1:T
    x = A_x_tt(:,t)./A_x_tt(:,1);
    x = x(~isnan(x));              % drop NaNs if any

    p25(t) = prctile(x,25);
    p50(t) = prctile(x,50);
    p75(t) = prctile(x,75);
end

figure('Position',[200 200 650 450]);
hold on

x = years;

% ---- Shaded 25–75 band ----
fill([x fliplr(x)], ...
     [p25' fliplr(p75')], ...
     [0.8 0.85 1], ...
     'EdgeColor','none','FaceAlpha',0.8);

% ---- Percentile bounds (dashed) ----
plot(x, p25, 'b--', 'LineWidth',1.5);
plot(x, p75, 'b--', 'LineWidth',1.5);

% ---- Median (solid) ----
plot(x, p50, 'b-', 'LineWidth',2);

xlabel('Year');
ylabel('Investment efficiency');
title({'Distribution of Investment Efficiency by Year relative to 1995', ...
       '(25th percentile, median and 75th percentile)'});

xlim([min(x) max(x)]);
grid on
box on

saveas(gcf, 'investment_efficiency_distribution.png');



%% Fig 3

years = 1995:2011;
T = length(years);

p25 = zeros(T,3);
p50 = zeros(T,3);
p75 = zeros(T,3);



for t = 1:T
    % Goods
    p25(t,1) = prctile(T_g_tt(:,t)./T_g_tt(:,1),25);
    p50(t,1) = prctile(T_g_tt(:,t)./T_g_tt(:,1),50);
    p75(t,1) = prctile(T_g_tt(:,t)./T_g_tt(:,1),75);

    % Low-skill services
    p25(t,2) = prctile(T_l_tt(:,t)./T_l_tt(:,1),25);
    p50(t,2) = prctile(T_l_tt(:,t)./T_l_tt(:,1),50);
    p75(t,2) = prctile(T_l_tt(:,t)./T_l_tt(:,1),75);

    % High-skill services
    p25(t,3) = prctile(T_h_tt(:,t)./T_h_tt(:,1),25);
    p50(t,3) = prctile(T_h_tt(:,t)./T_h_tt(:,1),50);
    p75(t,3) = prctile(T_h_tt(:,t)./T_h_tt(:,1),75);
end

sector_names = {'Goods','Low-skill services','High-skill services'};

figure('Position',[100 100 1000 350]);

for s = 1:3
    subplot(1,3,s); hold on

    x = years;

    % ---- Shaded 25–75 band ----
    fill([x fliplr(x)], ...
         [p25(:,s)' fliplr(p75(:,s)')], ...
         [0.8 0.85 1], ...
         'EdgeColor','none','FaceAlpha',0.8);

    % ---- Percentile lines ----
    plot(x, p25(:,s), 'b--', 'LineWidth',1.2);
    plot(x, p75(:,s), 'b--', 'LineWidth',1.2);

    % ---- Median ----
    plot(x, p50(:,s), 'b-', 'LineWidth',2);

    title(sector_names{s});
    xlabel('Year');
    if s==1
        ylabel('Trade cost');
    end
    xlim([min(x) max(x)]);
    grid on
end

sgtitle({'Distribution of Productivity by Sector and Year relative to 1995', ...
         '(25th percentile, median and 75th percentile)'});

saveas(gcf, 'productivity_distribution_3sectors.png');

%% sectoral export and import of India relative to 1995
figure('Position', get(0, 'Screensize'))
subplot(2,3,1)
plot(1995:1994+ye,ex_5_1_g_tt(1:ye)./ex_5_1_g_tt(1),'-b',1995:1994+ye,ex_5_2_g_tt(1:ye)./ex_5_2_g_tt(1),'-r',1995:1994+ye,ex_5_3_g_tt(1:ye)./ex_5_3_g_tt(1),'-g',1995:1994+ye,ex_5_4_g_tt(1:ye)./ex_5_4_g_tt(1),'-k',1995:1994+ye,ex_5_6_g_tt(1:ye)./ex_5_6_g_tt(1),'-c',1995:1994+ye,ex_5_7_g_tt(1:ye)./ex_5_7_g_tt(1),'-m',1995:1994+ye,ex_5_8_g_tt(1:ye)./ex_5_8_g_tt(1),'--k',1995:1994+ye,ex_5_9_g_tt(1:ye)./ex_5_9_g_tt(1),'-y','Linewidth',2.5)
title('Goods Export')
legend('To AUS','To BRA','To CHN','To DEU','To JPN','To MEX','To ROW','To USA','location','best')
grid on

subplot(2,3,2)
plot(1995:1994+ye,ex_5_1_l_tt(1:ye)./ex_5_1_l_tt(1),'-b',1995:1994+ye,ex_5_2_l_tt(1:ye)./ex_5_2_l_tt(1),'-r',1995:1994+ye,ex_5_3_l_tt(1:ye)./ex_5_3_l_tt(1),'-g',1995:1994+ye,ex_5_4_l_tt(1:ye)./ex_5_4_l_tt(1),'-k',1995:1994+ye,ex_5_6_l_tt(1:ye)./ex_5_6_l_tt(1),'-c',1995:1994+ye,ex_5_7_l_tt(1:ye)./ex_5_7_l_tt(1),'-m',1995:1994+ye,ex_5_8_l_tt(1:ye)./ex_5_8_l_tt(1),'--k',1995:1994+ye,ex_5_9_l_tt(1:ye)./ex_5_9_l_tt(1),'-y','Linewidth',2.5)
title('Low-skill service Export')
grid on

subplot(2,3,3)
plot(1995:1994+ye,ex_5_1_h_tt(1:ye)./ex_5_1_h_tt(1),'-b',1995:1994+ye,ex_5_2_h_tt(1:ye)./ex_5_2_h_tt(1),'-r',1995:1994+ye,ex_5_3_h_tt(1:ye)./ex_5_3_h_tt(1),'-g',1995:1994+ye,ex_5_4_h_tt(1:ye)./ex_5_4_h_tt(1),'-k',1995:1994+ye,ex_5_6_h_tt(1:ye)./ex_5_6_h_tt(1),'-c',1995:1994+ye,ex_5_7_h_tt(1:ye)./ex_5_7_h_tt(1),'-m',1995:1994+ye,ex_5_8_h_tt(1:ye)./ex_5_8_h_tt(1),'--k',1995:1994+ye,ex_5_9_h_tt(1:ye)./ex_5_9_h_tt(1),'-y','Linewidth',2.5)
title('High-skill service Export')
grid on

subplot(2,3,4)
plot(1995:1994+ye,im_5_1_g_tt(1:ye)./im_5_1_g_tt(1),'-b',1995:1994+ye,im_5_2_g_tt(1:ye)./im_5_2_g_tt(1),'-r',1995:1994+ye,im_5_3_g_tt(1:ye)./im_5_3_g_tt(1),'-g',1995:1994+ye,im_5_4_g_tt(1:ye)./im_5_4_g_tt(1),'-k',1995:1994+ye,im_5_6_g_tt(1:ye)./im_5_6_g_tt(1),'-c',1995:1994+ye,im_5_7_g_tt(1:ye)./im_5_7_g_tt(1),'-m',1995:1994+ye,im_5_8_g_tt(1:ye)./im_5_8_g_tt(1),'--k',1995:1994+ye,im_5_9_g_tt(1:ye)./im_5_9_g_tt(1),'-y','Linewidth',2.5)
title('Goods Import')
legend('From AUS','From BRA','From CHN','From DEU','From JPN','From MEX','From ROW','From USA','location','best')
grid on

subplot(2,3,5)
plot(1995:1994+ye,im_5_1_l_tt(1:ye)./im_5_1_l_tt(1),'-b',1995:1994+ye,im_5_2_l_tt(1:ye)./im_5_2_l_tt(1),'-r',1995:1994+ye,im_5_3_l_tt(1:ye)./im_5_3_l_tt(1),'-g',1995:1994+ye,im_5_4_l_tt(1:ye)./im_5_4_l_tt(1),'-k',1995:1994+ye,im_5_6_l_tt(1:ye)./im_5_6_l_tt(1),'-c',1995:1994+ye,im_5_7_l_tt(1:ye)./im_5_7_l_tt(1),'-m',1995:1994+ye,im_5_8_l_tt(1:ye)./im_5_8_l_tt(1),'--k',1995:1994+ye,im_5_9_l_tt(1:ye)./im_5_9_l_tt(1),'-y','Linewidth',2.5)
title('Low-skill service Import')
grid on

subplot(2,3,6)
plot(1995:1994+ye,im_5_1_h_tt(1:ye)./im_5_1_h_tt(1),'-b',1995:1994+ye,im_5_2_h_tt(1:ye)./im_5_2_h_tt(1),'-r',1995:1994+ye,im_5_3_h_tt(1:ye)./im_5_3_h_tt(1),'-g',1995:1994+ye,im_5_4_h_tt(1:ye)./im_5_4_h_tt(1),'-k',1995:1994+ye,im_5_6_h_tt(1:ye)./im_5_6_h_tt(1),'-c',1995:1994+ye,im_5_7_h_tt(1:ye)./im_5_7_h_tt(1),'-m',1995:1994+ye,im_5_8_h_tt(1:ye)./im_5_8_h_tt(1),'--k',1995:1994+ye,im_5_9_h_tt(1:ye)./im_5_9_h_tt(1),'-y','Linewidth',2.5)
title('High-skill service Import')
grid on

sgtitle('Figure. Bilateral export and import share as in Table 1')
saveas(gcf, 'plots_transition/outcome/01_export_import_ind.png');
close(gcf)
