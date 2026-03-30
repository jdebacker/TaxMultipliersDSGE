// tax_dsge.mod
// =====================================================================
// Zubairy (2014, IER) extended with full tax structure — Bayesian
// estimation version (standalone, no external .m file).
//
// Base: "On Fiscal Multipliers: Estimates from a Medium Scale DSGE Model"
//       International Economic Review, 55(1), 169-195.
//
// Extension adds tau_c (consumption tax), tau_ic (investment tax credit),
//   tau_i (interest income tax), tau_d (dividend tax), e_tau (expensing
//   rate, stochastic), plus tax basis tracking (k_tau) and firm dividend
//   accounting (div).
//
// Estimation:
//   - Steady state computed in steady_state_model + model-local (#)
//     variables so derived quantities update at every MCMC draw.
//   - Data file: US_Data_Matlab.mat (200 quarterly obs, 14 variables, 12 used as observables)
//   - prefilter = 1 demeaning handles SS ≠ sample mean mismatches
//
// Design:
//   - Tax channels enter HH equations following tax_dsge_fiscal_shocks_b
//   - Tax rate processes follow Zubairy-style rules (AR + debt/output feedback)
//   - All non-tax structure (deep habits, Rotemberg pricing, wage rigidity,
//     utilization, investment adjustment costs, monetary rule) is pure Zubairy
//   - Capital is predetermined: production uses k(-1)
//   - Present value multipliers computed outside Dynare (Matlab post-processing)
//
// /***************************************************************/
// /* Preamble                                                    */
// /***************************************************************/
//
// /* Endogenous Variables (46 structural + 13 measurement = 59)  */
//
// --- Macro aggregates ---
// c        = private consumption
// g        = government purchases
// i        = investment
// k        = physical capital stock (predetermined: k(-1) in production)
// u        = capital utilization rate
// h        = hours worked
// y        = output (GDP)
//
// --- Deep habits ---
// xc       = habit-adjusted private consumption surplus
// xg       = habit-adjusted government consumption surplus
// sC       = private consumption habit stock
// sG       = government consumption habit stock
//
// --- Prices and shadow values ---
// lambda   = marginal utility of consumption (shadow value of income)
// q        = Tobin's q (shadow price of installed capital)
// rk       = rental rate of capital services
// w        = real wage
// mc       = real marginal cost
//
// --- Nominal variables ---
// pi       = gross price inflation
// piw      = gross wage inflation
// R        = gross nominal interest rate (policy rate)
//
// --- Government ---
// b        = real government debt
// tr       = lump-sum transfers to households
// taxrev   = total government tax revenue
//
// --- Tax rates (Zubairy original) ---
// tau_k    = capital income tax rate
// tau_w    = wage/labor income tax rate
//
// --- Tax rates (extension) ---
// tau_c    = consumption tax rate
// tau_ic   = investment tax credit rate
// tau_i    = interest income tax rate
// tau_d    = dividend income tax rate
//
// --- Tax accounting ---
// k_tau    = tax basis of capital stock
// v_tau    = shadow value of tax basis (PV of future depreciation deductions)
// div      = firm dividends (markup profits, positive in SS since mc < 1)
//
// --- Tax depreciation (stochastic) ---
// e_tau    = immediate expensing rate (promoted from parameter to AR(1) process)
//
// --- Deep habits pricing ---
// nu_c     = shadow value of private consumption customer base
// nu_g     = shadow value of government consumption customer base
// nu_i     = shadow value of investment goods pricing
//
// --- Exogenous processes ---
// z        = technology shock (level, AR(1))
// d        = preference shock (level, AR(1))
// mu       = investment efficiency shock (level, AR(1))
//
// --- Wage rigidity ---
// wtilde   = wage markup
//
// --- Fiscal revenue accounting ---
// rev_c    = consumption tax revenue
// rev_w    = wage/labor tax revenue
// rev_k    = capital income tax revenue
// rev_i    = interest income tax revenue
// rev_d    = dividend tax revenue
// spend_ic = investment tax credit outlays
// rev_net  = net tax revenue (sum of revenues minus ITC outlays)
//
// --- Measurement variables (13) ---
// chat      = 100 × (c - c_ss)/c_ss          [% deviation]
// ihat      = 100 × (i - i_ss)/i_ss          [% deviation]
// ghat      = 100 × (g - g_ss)/g_ss          [% deviation]
// pi_obs    = pi - pi_ss                      [level deviation]
// r_obs     = R - R_ss                        [level deviation]
// tau_k_hat = (tau_k - tauk_ss)/tauk_ss       [fractional deviation]
// tau_l_hat = (tau_w - tauw_ss)/tauw_ss       [fractional deviation]
// bhat      = 100 × (b - b_ss)/b_ss          [% deviation]
// trhat     = 100 × (tr - tr_ss)/tr_ss       [% deviation]
// tau_c_hat = (tau_c - tauc_ss)/tauc_ss       [fractional deviation]
// tau_d_hat = (tau_d - taud_ss)/taud_ss       [fractional deviation]
// tau_itc_hat = tau_ic - tauic_ss             [level deviation, SS may be 0]
// e_tau_hat   = e_tau - etau_ss               [level deviation, SS may be 0]
//
// /***************************************************************/
//
// /* Exogenous Shocks (13)                                       */
//
// eps_g    = government spending shock
// eps_tr   = transfer shock
// eps_k    = capital tax rate shock (structural innovation)
// eps_w    = labor tax rate shock (structural innovation)
// eps_z    = technology shock
// eps_d    = preference shock
// eps_mu   = investment efficiency shock
// eps_m    = monetary policy shock
// eps_c    = consumption tax rate shock
// eps_ic   = investment tax credit shock
// eps_ti   = interest income tax rate shock
// eps_td   = dividend tax rate shock
// eps_etau = expensing rate shock
//
// /***************************************************************/
//
// /* Parameters                                                  */
//
// --- Preferences ---
// beta     = discount factor
// gamma    = risk aversion / inverse EIS
//
// --- Technology ---
// theta    = capital share in Cobb-Douglas production
// delta    = economic depreciation rate of capital
// kappa    = investment adjustment cost parameter
// sig_u    = curvature of capital utilization cost
//
// --- Goods market ---
// eta      = elasticity of substitution across goods varieties
// alphaP   = Rotemberg price adjustment cost
//
// --- Labor market ---
// etaw     = elasticity of substitution across labor varieties
// alphaW   = Rotemberg wage adjustment cost
//
// --- Deep habits ---
// bc       = deep habit parameter, private consumption
// thetac   = habit stock adjustment speed, private consumption
// bg       = deep habit parameter, government consumption
// thetag   = habit stock adjustment speed, government consumption
//
// --- Steady-state tax rates (calibrated) ---
// tauk_ss  = SS capital income tax rate
// tauw_ss  = SS wage/labor income tax rate
// tauc_ss  = SS consumption tax rate
// tauic_ss = SS investment tax credit rate
// taui_ss  = SS interest income tax rate
// taud_ss  = SS dividend income tax rate
// etau_ss  = SS immediate expensing rate
//
// --- Tax depreciation (structural) ---
// delta_tau = tax depreciation rate
//
// --- Calibration targets ---
// pi_ss    = SS gross inflation
// h_ss     = SS hours worked
// g_y      = SS government spending share of output
// b_y      = SS debt-to-output ratio
//
// --- Taylor rule ---
// alphaR   = interest rate smoothing
// alphapi  = Taylor rule response to inflation
// alphaY   = Taylor rule response to output
//
// --- Exogenous process persistence ---
// rho_z    = AR(1) coefficient, technology
// rho_d    = AR(1) coefficient, preference
// rho_mu   = AR(1) coefficient, investment efficiency
//
// --- Fiscal rule: government spending ---
// rho_g    = AR(1) persistence
// rho_g_b  = response to lagged debt (negative sign in rule)
// rho_g_y  = response to lagged output
//
// --- Fiscal rule: transfers ---
// rho_tr   = AR(1) persistence
// rho_tr_b = response to lagged debt (negative sign in rule)
// rho_tr_y = response to lagged output
//
// --- Fiscal rule: capital tax ---
// rho_k    = AR(1) persistence
// rho_k_b  = response to lagged debt
// rho_k_y  = response to lagged output
//
// --- Fiscal rule: wage tax ---
// rho_w    = AR(1) persistence
// rho_w_b  = response to lagged debt
// rho_w_y  = response to lagged output
// phi      = cross-correlation: eps_k enters tau_w rule
// phi_wd   = cross-correlation: eps_w enters tau_d rule (tau_l/tau_d corr = 0.90)
//
// --- Fiscal rule: consumption tax ---
// rho_c    = AR(1) persistence
// rho_c_b  = response to lagged debt
// rho_c_y  = response to lagged output
//
// --- Fiscal rule: investment tax credit ---
// rho_ic   = AR(1) persistence (level-deviation form)
// rho_ic_b = response to lagged debt
// rho_ic_y = response to lagged output
//
// --- Fiscal rule: interest income tax ---
// rho_ti   = AR(1) persistence
// rho_ti_b = response to lagged debt
// rho_ti_y = response to lagged output
//
// --- Fiscal rule: dividend tax ---
// rho_td   = AR(1) persistence
// rho_td_b = response to lagged debt
// rho_td_y = response to lagged output
//
// --- Fiscal rule: expensing rate ---
// rho_etau   = AR(1) persistence (level-deviation form)
// rho_etau_b = response to lagged debt
// rho_etau_y = response to lagged output
//
// /***************************************************************/

var
    c g i k u h y
    xc xg sC sG
    lambda q rk w mc
    pi piw R
    b tr taxrev
    tau_k tau_w
    tau_c tau_ic tau_i tau_d
    k_tau v_tau div
    e_tau
    nu_c nu_g nu_i
    z d mu
    wtilde
    rev_c rev_w rev_k rev_i rev_d spend_ic rev_net
    // Measurement variables (match data column names)
    chat ihat ghat pi_obs r_obs tau_k_hat tau_l_hat bhat
    trhat tau_c_hat tau_d_hat tau_itc_hat e_tau_hat
;

varexo
    eps_g eps_tr
    eps_k eps_w
    eps_z eps_d eps_mu
    eps_m
    eps_c eps_ic eps_ti eps_td
    eps_etau
;

parameters
    beta delta theta eta etaw gamma
    kappa sig_u
    alphaP alphaW
    bc thetac bg thetag
    tauk_ss tauw_ss pi_ss h_ss g_y b_y
    alphaR alphapi alphaY
    rho_z rho_d rho_mu
    rho_g rho_g_b rho_g_y
    rho_tr rho_tr_b rho_tr_y
    rho_k rho_k_b rho_k_y
    rho_w rho_w_b rho_w_y
    phi phi_wd
    tauc_ss  rho_c  rho_c_b  rho_c_y
    tauic_ss rho_ic rho_ic_b rho_ic_y
    taui_ss  rho_ti rho_ti_b rho_ti_y
    taud_ss  rho_td rho_td_b rho_td_y
    etau_ss  rho_etau rho_etau_b rho_etau_y
    delta_tau divshare
;


// =====================================================================
// Calibration: fixed parameters (not estimated)
// =====================================================================

beta    = 1/(1.03^(1/4));
delta   = 0.025;
theta   = 0.30;
eta     = 5.3;
etaw    = 21.0;
pi_ss   = 1.0;
h_ss    = 0.50;
g_y     = 0.18;
b_y     = 0.33;
tauk_ss = 0.41;
tauw_ss = 0.23;

// Extension: SS tax rates (calibrated)
tauc_ss  = 0.10;
tauic_ss = 0.0;
taui_ss  = 0.20;
taud_ss  = 0.15;
etau_ss  = 0.0;
delta_tau = 0.025;

// Dividend share target (fraction of output going to firm profits)
// Old model (Calvo explicit profits): ~0.11. Set slightly higher to
// account for Rotemberg pricing absorbing some markup.
divshare  = 0.11;


// =====================================================================
// Estimated parameters: initial values
// =====================================================================

gamma   = 2.49;
sig_u   = 2.52;
kappa   = 1.99;
alphaP  = 48.90;
alphaW  = 102.0;
bc      = 0.96;
thetac  = 0.53;
bg      = 0.85;
thetag  = 0.98;
alphaR  = 0.60;
alphapi = 1.67;
alphaY  = 0.09;
rho_z   = 0.79;
rho_d   = 0.72;
rho_mu  = 0.67;
rho_g   = 0.92;
rho_tr  = 0.64;
rho_k   = 0.91;
rho_w   = 0.90;
phi     = 0.22;
phi_wd  = 0.25;
rho_k_b = 0.017;
rho_w_b = 0.020;
rho_g_b = 0.009;
rho_tr_b = 0.439;
rho_k_y = 0.148;
rho_w_y = 0.132;
rho_g_y = -0.039;
rho_tr_y = -0.079;

// Extension fiscal rule initial values
rho_c    = 0.90;  rho_c_b  = 0.0;    rho_c_y  = 0.0;
rho_td   = 0.90;  rho_td_b = 0.0;    rho_td_y = 0.0;
rho_ic   = 0.90;  rho_ic_b = 0.0;    rho_ic_y = 0.0;
rho_ti   = 0.90;  rho_ti_b = 0.0;    rho_ti_y = 0.0;
rho_etau = 0.90;  rho_etau_b = 0.0;  rho_etau_y = 0.0;


// =====================================================================
// MODEL (46 structural + 13 measurement = 59 equations)
// =====================================================================
// Model-local (#) variables recompute derived SS quantities at every
// evaluation, so they update when estimated parameters change during MCMC.
// This replaces the role of the external _steadystate.m file.
// =====================================================================
model;

    // ---- Derived SS quantities as model-local variables ----
    // Dividends = divshare * Y_ss > 0. A partial fixed cost psi absorbs
    // the remainder of markup profits. With partial psi, the actual
    // demand shares are: I/Y = sharei*(1-divshare), G/Y = g_y,
    // C/Y = 1 - G/Y - I/Y. These feed the deep habits pricing condition.
    #v_tau_ss = beta*tauk_ss*delta_tau / (1 - beta*(1 - delta_tau));
    #q_ss     = 1 - tauic_ss - tauk_ss*etau_ss - v_tau_ss*(1 - etau_ss);
    #rk_ss    = q_ss*(1/beta - 1 + delta) / (1 - tauk_ss);
    #sharei   = delta*theta / rk_ss;
    #sharec   = 1 - g_y - sharei*(1 - divshare);
    #aac      = 1 - bc;
    #aag      = 1 - bg;
    #bbc      = (beta*bc*(thetac - 1)) / (beta*thetac - 1) - 1;
    #bbg      = (beta*bg*(thetag - 1)) / (beta*thetag - 1) - 1;
    #mc_ss    = 1/(eta*(sharec*aac/bbc + g_y*aag/bbg - sharei*(1 - divshare))) + 1;
    #K_ss     = (rk_ss/(mc_ss*theta))^(1/(theta - 1)) * h_ss;
    #w_ss     = mc_ss*(1 - theta)*(K_ss/h_ss)^theta;
    #Y_ss     = (w_ss*h_ss + rk_ss*K_ss) / (1 - divshare);
    #psi_val  = K_ss^theta*h_ss^(1 - theta) - Y_ss;
    #I_ss     = delta*K_ss;
    #G_ss     = g_y*Y_ss;
    #C_ss     = Y_ss - I_ss - G_ss;
    #k_tau_ss = I_ss*(1 - etau_ss)/delta_tau;
    #R_ss     = 1 + (pi_ss/beta - 1) / (1 - taui_ss);
    #b_ss     = b_y*Y_ss;
    #div_ss   = divshare*Y_ss;
    #taxrev_ss = tauc_ss*C_ss
               + tauw_ss*w_ss*h_ss
               + tauk_ss*(rk_ss*K_ss - delta_tau*k_tau_ss - etau_ss*I_ss)
               + taui_ss*(R_ss - 1)*b_ss/pi_ss
               + taud_ss*div_ss
               - tauic_ss*I_ss;
    #tr_ss    = b_ss*(1 - R_ss/pi_ss) - G_ss + taxrev_ss;
    #xc_ss    = (1 - bc)*C_ss;
    #xg_ss    = (1 - bg)*G_ss;
    #wtilde_ss = etaw/((etaw - 1)*(1 - tauw_ss));
    #a_share  = xc_ss*(1 + tauc_ss) / (xc_ss*(1 + tauc_ss) + (w_ss/wtilde_ss)*(1 - h_ss));
    #lambda_ss = a_share * ( (xc_ss^a_share)*((1 - h_ss)^(1 - a_share)) )^(1 - gamma) / (xc_ss*(1 + tauc_ss));
    #gamma1   = (1 - tauk_ss)*rk_ss;


    // (1-4) Deep habits
    xc = c - bc*sC(-1);
    xg = g - bg*sG(-1);
    sC = thetac*sC(-1) + (1 - thetac)*c;
    sG = thetag*sG(-1) + (1 - thetag)*g;

    // (5) Marginal utility of consumption [MODIFIED: /(1+tau_c)]
    lambda = d * a_share * ( (xc^a_share)*((1 - h)^(1 - a_share)) )^(1 - gamma) / (xc*(1 + tau_c));

    // (6) Labor-leisure intratemporal condition
    d * (1 - a_share) * ( (xc^a_share)*((1 - h)^(1 - a_share)) )^(1 - gamma) / (1 - h) = lambda*w / wtilde;

    // (7-9) Production and factor prices
    y = z*(u*k(-1))^theta*h^(1 - theta) - psi_val - (alphaP/2)*(pi - pi_ss)^2;
    w  = mc*z*(1 - theta)*(u*k(-1))^theta*h^(-theta);
    rk = mc*z*theta*(u*k(-1))^(theta - 1)*h^(1 - theta);

    // (10) Resource constraint
    z*(u*k(-1))^theta*h^(1 - theta) - psi_val
        = c + g + i
        + (gamma1/(1 + sig_u))*(u^(1 + sig_u) - 1)*k(-1)
        + (alphaP/2)*(pi - pi_ss)^2
        + (alphaW/2)*(piw - pi_ss)^2*w;

    // (11) Capital accumulation
    k = (1 - delta)*k(-1) + i*(1 - (kappa/2)*(mu*i/i(-1) - 1)^2);

    // (12) Tax basis of capital [EXTENSION]
    k_tau = (1 - delta_tau)*k_tau(-1) + i*(1 - e_tau);

    // (12b) Shadow value of tax basis [EXTENSION]
    //   v_tau is the PV (in consumption units) of an extra unit of tax basis,
    //   which yields future depreciation deductions at rate delta_tau taxed at tau_k.
    lambda*v_tau = beta*lambda(+1) * (
        (1 - delta_tau)*v_tau(+1) + tau_k(+1)*delta_tau
    );

    // (13) Utilization FOC
    (1 - tau_k)*rk = gamma1*u^sig_u;

    // (14) Capital Euler
    lambda*q = beta*lambda(+1) * (
        (1 - tau_k(+1))*rk(+1)*u(+1)
        - (gamma1/(1 + sig_u))*(u(+1)^(1 + sig_u) - 1)
        + q(+1)*(1 - delta)
    );

    // (15) Investment FOC [MODIFIED: (1-tau_ic-tau_k*e_tau) on LHS]
    //   LHS: marginal cost of investing $1, net of:
    //     - tau_ic: investment tax credit (direct subsidy)
    //     - tau_k*e_tau: immediate expensing deduction (fraction e_tau is
    //       expensed at capital tax rate tau_k, reducing cost today)
    //   The FUTURE depreciation benefit on the non-expensed portion
    //   (1-e_tau) is carried by v_tau (shadow value of tax basis).
    //   Correspondingly, the Zubairy-style depreciation terms are removed
    //   from the utilization FOC and the capital Euler to avoid double counting.
    lambda*(1 - tau_ic - tau_k*e_tau) = lambda*q * (
        1 - (kappa/2)*(mu*i/i(-1) - 1)^2
          - (mu*i/i(-1)) * kappa*(mu*i/i(-1) - 1)
    )
    + beta*lambda(+1)*q(+1) * (
        mu(+1)*(i(+1)/i)^2 * kappa*(mu(+1)*i(+1)/i - 1)
    )
    + lambda*v_tau*(1 - e_tau);

    // (16) Bond Euler [MODIFIED: after-tax return with tau_i(+1)]
    lambda = beta * (1 + (R - 1)*(1 - tau_i(+1))) * lambda(+1) / pi(+1);

    // (17) Wage inflation identity
    piw = (w/w(-1))*pi;

    // (18) Wage Phillips curve
    (etaw - 1)*(1 - tau_w)*h + alphaW*piw*(piw - pi_ss) - etaw*h/wtilde
        = alphaW*beta*(lambda(+1)/lambda)*piw(+1)*(piw(+1) - pi_ss);

    // (19-21) Deep habits pricing
    1 - mc - nu_i = 0;

    (1 - mc - nu_c)/(thetac - 1)
        = beta*(lambda(+1)/lambda) * (
            bc*nu_c(+1) + thetac/(thetac - 1)*(1 - mc(+1) - nu_c(+1))
        );

    (1 - mc - nu_g)/(thetag - 1)
        = beta*(lambda(+1)/lambda) * (
            bg*nu_g(+1) + thetag/(thetag - 1)*(1 - mc(+1) - nu_g(+1))
        );

    // (22) Price Phillips curve
    eta*(nu_c*xc + nu_g*xg + nu_i*(y - c - g))
        + alphaP*pi*(pi - pi_ss) - y
        = alphaP*beta*(lambda(+1)/lambda)*pi(+1)*(pi(+1) - pi_ss);

    // (23) Firm dividends [EXTENSION]
    div = y - w*h - rk*u*k(-1);

    // (24-25) Government budget and tax revenue [MODIFIED]
    b = R(-1)*b(-1)/pi + g + tr - taxrev;

    taxrev = tau_c*c
           + tau_w*w*h
           + tau_k*(rk*u*k(-1) - delta_tau*k_tau(-1) - e_tau*i)
           + tau_i*(R(-1) - 1)*b(-1)/pi
           + tau_d*div
           - tau_ic*i;

    // (26-29) Zubairy fiscal rules (unchanged)
    log(tau_k/tauk_ss) = rho_k*log(tau_k(-1)/tauk_ss)
        + rho_k_b*log(b(-1)/b_ss) + rho_k_y*log(y(-1)/Y_ss) + eps_k;

    log(tau_w/tauw_ss) = rho_w*log(tau_w(-1)/tauw_ss)
        + rho_w_b*log(b(-1)/b_ss) + rho_w_y*log(y(-1)/Y_ss) + phi*eps_k + eps_w;

    log(g/G_ss) = rho_g*log(g(-1)/G_ss)
        - rho_g_b*log(b(-1)/b_ss) + rho_g_y*log(y(-1)/Y_ss) + eps_g;

    log(tr/tr_ss) = rho_tr*log(tr(-1)/tr_ss)
        - rho_tr_b*log(b(-1)/b_ss) + rho_tr_y*log(y(-1)/Y_ss) + eps_tr;

    // (30) Taylor rule (unchanged)
    log(R/R_ss) = alphaR*log(R(-1)/R_ss)
        + (1 - alphaR)*(alphapi*log(pi/pi_ss) + alphaY*log(y/Y_ss)) + eps_m;

    // (31-34) Extension fiscal rules [EXTENSION]
    log(tau_c/tauc_ss) = rho_c*log(tau_c(-1)/tauc_ss)
        + rho_c_b*log(b(-1)/b_ss) + rho_c_y*log(y(-1)/Y_ss) + eps_c;

    tau_ic - tauic_ss = rho_ic*(tau_ic(-1) - tauic_ss)
        + rho_ic_b*log(b(-1)/b_ss) + rho_ic_y*log(y(-1)/Y_ss) + eps_ic;

    log(tau_i/taui_ss) = rho_ti*log(tau_i(-1)/taui_ss)
        + rho_ti_b*log(b(-1)/b_ss) + rho_ti_y*log(y(-1)/Y_ss) + eps_ti;

    log(tau_d/taud_ss) = rho_td*log(tau_d(-1)/taud_ss)
        + rho_td_b*log(b(-1)/b_ss) + rho_td_y*log(y(-1)/Y_ss) + phi_wd*eps_w + eps_td;

    // (35) Expensing rate fiscal rule [NEW: e_tau promoted to stochastic]
    //   Level-deviation form (like tau_ic) since SS may be zero.
    e_tau - etau_ss = rho_etau*(e_tau(-1) - etau_ss)
        + rho_etau_b*log(b(-1)/b_ss) + rho_etau_y*log(y(-1)/Y_ss) + eps_etau;

    // (36-38) Exogenous processes (unchanged)
    log(z)  = rho_z*log(z(-1))   + eps_z;
    log(d)  = rho_d*log(d(-1))   + eps_d;
    log(mu) = rho_mu*log(mu(-1)) + eps_mu;

    // (39-45) Fiscal revenue accounting [EXTENSION]
    rev_c    = tau_c*c;
    rev_w    = tau_w*w*h;
    rev_k    = tau_k*(rk*u*k(-1) - delta_tau*k_tau(-1) - e_tau*i);
    rev_i    = tau_i*(R(-1) - 1)*b(-1)/pi;
    rev_d    = tau_d*div;
    spend_ic = tau_ic*i;
    rev_net  = rev_c + rev_w + rev_k + rev_i + rev_d - spend_ic;

    // (46-58) Measurement equations
    // Quantities: 100 × percent deviation from SS
    // Tax rates:  fractional deviation from SS (or level deviation when SS ≈ 0)
    // Rates:      level deviations from SS
    chat        = 100 * (c - C_ss) / C_ss;
    ihat        = 100 * (i - I_ss) / I_ss;
    ghat        = 100 * (g - G_ss) / G_ss;
    pi_obs      = pi - pi_ss;
    r_obs       = R - R_ss;
    tau_k_hat   = (tau_k - tauk_ss) / tauk_ss;
    tau_l_hat   = (tau_w - tauw_ss) / tauw_ss;
    bhat        = 100 * (b - b_ss) / b_ss;
    trhat       = 100 * (tr - tr_ss) / tr_ss;
    tau_c_hat   = (tau_c - tauc_ss) / tauc_ss;
    tau_d_hat   = (tau_d - taud_ss) / taud_ss;
    tau_itc_hat = tau_ic - tauic_ss;
    e_tau_hat   = e_tau - etau_ss;

end;


// =====================================================================
// STEADY STATE MODEL
// =====================================================================
// Recomputes SS analytically from current parameters at each MCMC draw.
// Uses the same formulas as the # model-local variables above.
// =====================================================================
steady_state_model;

    // Derived SS quantities (divshare/psi approach, positive dividends)
    v_tau_ss_v = beta*tauk_ss*delta_tau / (1 - beta*(1 - delta_tau));
    q_ss_v     = 1 - tauic_ss - tauk_ss*etau_ss - v_tau_ss_v*(1 - etau_ss);
    rk_ss_v    = q_ss_v*(1/beta - 1 + delta) / (1 - tauk_ss);
    sharei_v   = delta*theta / rk_ss_v;
    sharec_v   = 1 - g_y - sharei_v*(1 - divshare);
    aac_v      = 1 - bc;
    aag_v      = 1 - bg;
    bbc_v      = (beta*bc*(thetac - 1)) / (beta*thetac - 1) - 1;
    bbg_v      = (beta*bg*(thetag - 1)) / (beta*thetag - 1) - 1;
    mc_ss_v    = 1/(eta*(sharec_v*aac_v/bbc_v + g_y*aag_v/bbg_v - sharei_v*(1 - divshare))) + 1;
    K_ss_v     = (rk_ss_v/(mc_ss_v*theta))^(1/(theta - 1)) * h_ss;
    w_ss_v     = mc_ss_v*(1 - theta)*(K_ss_v/h_ss)^theta;
    Y_ss_v     = (w_ss_v*h_ss + rk_ss_v*K_ss_v) / (1 - divshare);
    I_ss_v     = delta*K_ss_v;
    G_ss_v     = g_y*Y_ss_v;
    C_ss_v     = Y_ss_v - I_ss_v - G_ss_v;
    k_tau_ss_v = I_ss_v*(1 - etau_ss)/delta_tau;
    R_ss_v     = 1 + (pi_ss/beta - 1) / (1 - taui_ss);
    b_ss_v     = b_y*Y_ss_v;
    div_ss_v   = divshare*Y_ss_v;
    taxrev_ss_v = tauc_ss*C_ss_v
                + tauw_ss*w_ss_v*h_ss
                + tauk_ss*(rk_ss_v*K_ss_v - delta_tau*k_tau_ss_v - etau_ss*I_ss_v)
                + taui_ss*(R_ss_v - 1)*b_ss_v/pi_ss
                + taud_ss*div_ss_v
                - tauic_ss*I_ss_v;
    tr_ss_v    = b_ss_v*(1 - R_ss_v/pi_ss) - G_ss_v + taxrev_ss_v;
    xc_ss_v    = (1 - bc)*C_ss_v;
    xg_ss_v    = (1 - bg)*G_ss_v;
    wtilde_ss_v = etaw/((etaw - 1)*(1 - tauw_ss));
    a_share_v  = xc_ss_v*(1 + tauc_ss) / (xc_ss_v*(1 + tauc_ss) + (w_ss_v/wtilde_ss_v)*(1 - h_ss));
    lambda_ss_v = a_share_v * ( (xc_ss_v^a_share_v)*((1 - h_ss)^(1 - a_share_v)) )^(1 - gamma) / (xc_ss_v*(1 + tauc_ss));

    // Assign endogenous variable SS values
    z  = 1;
    d  = 1;
    mu = 1;
    c  = C_ss_v;
    g  = G_ss_v;
    i  = I_ss_v;
    k  = K_ss_v;
    u  = 1;
    h  = h_ss;
    y  = Y_ss_v;
    xc = xc_ss_v;
    xg = xg_ss_v;
    sC = C_ss_v;
    sG = G_ss_v;
    lambda = lambda_ss_v;
    q  = q_ss_v;
    rk = rk_ss_v;
    w  = w_ss_v;
    mc = mc_ss_v;
    pi = pi_ss;
    piw = pi_ss;
    R  = R_ss_v;
    b  = b_ss_v;
    tr = tr_ss_v;
    taxrev = taxrev_ss_v;
    tau_k = tauk_ss;
    tau_w = tauw_ss;
    tau_c = tauc_ss;
    tau_ic = tauic_ss;
    tau_i = taui_ss;
    tau_d = taud_ss;
    e_tau = etau_ss;
    nu_c = (mc_ss_v - 1)/bbc_v;
    nu_g = (mc_ss_v - 1)/bbg_v;
    nu_i = 1 - mc_ss_v;
    wtilde = wtilde_ss_v;
    k_tau = k_tau_ss_v;
    v_tau = v_tau_ss_v;
    div = div_ss_v;
    rev_c = tauc_ss*C_ss_v;
    rev_w = tauw_ss*w_ss_v*h_ss;
    rev_k = tauk_ss*(rk_ss_v*K_ss_v - delta_tau*k_tau_ss_v - etau_ss*I_ss_v);
    rev_i = taui_ss*(R_ss_v - 1)*b_ss_v/pi_ss;
    rev_d = taud_ss*div_ss_v;
    spend_ic = tauic_ss*I_ss_v;
    rev_net = tauc_ss*C_ss_v + tauw_ss*w_ss_v*h_ss
            + tauk_ss*(rk_ss_v*K_ss_v - delta_tau*k_tau_ss_v - etau_ss*I_ss_v)
            + taui_ss*(R_ss_v - 1)*b_ss_v/pi_ss + taud_ss*div_ss_v - tauic_ss*I_ss_v;

    // Measurement variables (all zero in SS)
    chat        = 0;
    ihat        = 0;
    ghat        = 0;
    pi_obs      = 0;
    r_obs       = 0;
    tau_k_hat   = 0;
    tau_l_hat   = 0;
    bhat        = 0;
    trhat       = 0;
    tau_c_hat   = 0;
    tau_d_hat   = 0;
    tau_itc_hat = 0;
    e_tau_hat   = 0;

end;

steady;
check;


// =====================================================================
// OBSERVABLES — 12 observables, 13 shocks (over-identified by 1)
// =====================================================================
// Zubairy 8 + tau_c_hat + tau_d_hat + tau_itc_hat + e_tau_hat.
// trhat excluded: transfer data is noisy/unreliable. eps_tr is still
// in the model and identified indirectly through debt dynamics.
// The unmatched shock is eps_ti (no tau_i observable in the data).
// =====================================================================
varobs chat ihat ghat pi_obs r_obs tau_k_hat tau_l_hat bhat
       tau_c_hat tau_d_hat tau_itc_hat e_tau_hat;


// =====================================================================
// PRIORS
// =====================================================================
// Zubairy original parameters (28): estimated as in Table 1.
// All 5 extension fiscal rules are also estimated:
//   tau_c  (rho_c, rho_c_b, rho_c_y)
//   tau_ic (rho_ic, rho_ic_b, rho_ic_y)
//   tau_i  (rho_ti, rho_ti_b, rho_ti_y)
//   tau_d  (rho_td, rho_td_b, rho_td_y)
//   e_tau  (rho_etau, rho_etau_b, rho_etau_y)
// Plus phi_wd (tau_l/tau_d cross-correlation) and all 13 shock stds.
// Total: 57 estimated parameters.
// =====================================================================
estimated_params;

    // --- Deep habits ---
    bc,       beta_pdf,       0.7,   0.1;
    thetac,   beta_pdf,       0.8,   0.1;
    bg,       beta_pdf,       0.7,   0.1;
    thetag,   beta_pdf,       0.8,   0.1;

    // --- Structural ---
    gamma,    normal_pdf,     2.0,   1.0;
    sig_u,    normal_pdf,     2.5,   0.5;
    kappa,    normal_pdf,     2.0,   0.5;

    // --- Nominal rigidities ---
    alphaP,   normal_pdf,     17.0,  5.0;
    alphaW,   normal_pdf,     100.0, 30.0;

    // --- Monetary policy ---
    alphaR,   beta_pdf,       0.8,   0.2;  // 0.8,   0.15; //paper: 0.8,   0.2;
    alphapi,  normal_pdf,     1.6,   0.2;
    alphaY,   normal_pdf,     0.1,   0.05;

    // --- AR coefficients (Zubairy original) ---
    rho_k,    beta_pdf,       0.7,   0.2;
    rho_w,    beta_pdf,       0.7,   0.2;
    rho_g,    beta_pdf,       0.8,   0.2;  // 0.8,   0.15; //paper: 0.8,   0.2;
    rho_d,    beta_pdf,       0.7,   0.2;
    rho_tr,   beta_pdf,       0.7,   0.2;
    rho_z,    beta_pdf,       0.7,   0.2;
    rho_mu,   beta_pdf,       0.7,   0.2;

    // --- Fiscal rule: debt response (Zubairy original) ---
    rho_k_b,  gamma_pdf,      0.5,   0.25;
    rho_w_b,  gamma_pdf,      0.5,   0.25;
    rho_g_b,  gamma_pdf,      0.5,   0.25;
    rho_tr_b, gamma_pdf,      0.5,   0.25;

    // --- Fiscal rule: output response (Zubairy original) ---
    rho_k_y,  gamma_pdf,      0.15,  0.1;
    rho_w_y,  gamma_pdf,      0.15,  0.1;
    rho_g_y,  normal_pdf,    -0.05,  0.05;
    rho_tr_y, normal_pdf,    -0.1,   0.05;

    // --- Tax correlations ---
    phi,      normal_pdf,     0.25,  0.1;
    phi_wd,   normal_pdf,     0.25,  0.1;   // eps_w enters tau_d (data corr = 0.90)

    // --- AR coefficients (extension) ---
    rho_c,    beta_pdf,       0.7,   0.2;
    rho_ic,   beta_pdf,       0.7,   0.2;
    rho_ti,   beta_pdf,       0.7,   0.2;
    rho_td,   beta_pdf,       0.7,   0.2;
    rho_etau, beta_pdf,       0.7,   0.2;

    // --- Fiscal rule: debt response (extension) ---
    rho_c_b,    gamma_pdf,    0.5,   0.25;
    rho_ti_b,   gamma_pdf,    0.5,   0.25;
    rho_td_b,   gamma_pdf,    0.5,   0.25;
    rho_ic_b,   normal_pdf,   0.0,   0.25;
    rho_etau_b, normal_pdf,   0.0,   0.25;

    // --- Fiscal rule: output response (extension) ---
    rho_c_y,    normal_pdf,   0.0,   0.1;
    rho_ic_y,   normal_pdf,   0.0,   0.1;
    rho_ti_y,   normal_pdf,   0.0,   0.1;
    rho_td_y,   normal_pdf,   0.0,   0.1;
    rho_etau_y, normal_pdf,   0.0,   0.1;

    // --- Shock standard deviations (Zubairy original) ---
    stderr eps_k,   inv_gamma_pdf,  0.05,  0.1;
    stderr eps_w,   inv_gamma_pdf,  0.05,  0.1;
    stderr eps_g,   inv_gamma_pdf,  0.05,  0.1;
    stderr eps_d,   inv_gamma_pdf,  0.05,  0.1;
    stderr eps_tr,  inv_gamma_pdf,  0.05,  0.1;
    stderr eps_z,   inv_gamma_pdf,  0.05,  0.1;
    stderr eps_m,   inv_gamma_pdf,  0.05,  0.1;
    stderr eps_mu,  inv_gamma_pdf,  0.05,  0.1;

    // --- Shock standard deviations (extension) ---
    stderr eps_c,    inv_gamma_pdf, 0.05,  0.1;
    stderr eps_ic,   inv_gamma_pdf, 0.05,  0.1;
    stderr eps_ti,   inv_gamma_pdf, 0.05,  0.1;
    stderr eps_td,   inv_gamma_pdf, 0.05,  0.1;
    stderr eps_etau, inv_gamma_pdf, 0.05,  0.1;

end;


// =====================================================================
// ESTIMATION
// =====================================================================

estimation(
    datafile     = US_Data_Matlab,
    prefilter    = 1,
    lik_init     = 1,
    mh_replic    = 1500000, // 500000, 20000, 1500000
    mh_nblocks   = 2,
    mh_drop      = 0.333,
    mh_jscale    = 0.3,
    mode_compute = 6, // 4 is csminwel, 6 is a Monte Carlo based mode-finding routine
    // posterior_sampling_method = 6,
    // bayesian_irf,
    irf          = 40,
    nodisplay
);
