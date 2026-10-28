// =====================================================================
// tax_dsge_kgoods2_tifix: COPY of tax_dsge_kgoods2, 2026-10-22. Only change: the interest income
// tax is exogenous and constant (tau_i = taui_ss). Its rule, its shock
// eps_ti and its estimated parameters are commented out. Reason: tau_i
// is not observed (tau_i_hat dropped from varobs) and its rule was the
// slowest-mixing block in all four models, failing C1 in the feedback
// models. The interest income tax multiplier is no longer computed.
// =====================================================================
// =====================================================================
// tax_dsge_kgoods2_tifix: COPY made 2026-10-21 for re-estimation with the shared
// calibration of thetac, sig_u, alphaW, bg (see CALIBRATED 2026-10-21).
// Old estimation blocks disabled; RUN BLOCKS at the end of the file.
// =====================================================================
// tax_dsge_est_full.mod

// /***************************************************************/
// /* Preamble                                                    */
// /***************************************************************/
//
// /* Endogenous Variables (49 structural + 14 measurement = 63)  */
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
// --- Measurement variables (14) ---
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
// tau_i_hat = (tau_i - taui_ss)/taui_ss       [fractional deviation]
// tau_itc_hat = tau_ic - tauic_ss             [level deviation, SS may be 0]
// e_tau_hat   = e_tau - etau_ss               [level deviation, SS may be 0]
//
// /***************************************************************/
//
// /* Exogenous Shocks (14)                                       */
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
// eps_b    = government budget / debt residual shock
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
//
// --- Fiscal rule: transfers ---
// rho_tr   = AR(1) persistence
// rho_tr_b = response to lagged debt (negative sign in rule)
//
// --- Fiscal rule: capital tax ---
// rho_k    = AR(1) persistence
//
// --- Fiscal rule: wage tax ---
// rho_w    = AR(1) persistence
// phi      = cross-correlation: eps_k enters tau_w rule
//
// --- Fiscal rule: consumption tax ---
// rho_c    = AR(1) persistence
//
// --- Fiscal rule: investment tax credit ---
// rho_ic   = AR(1) persistence (level-deviation form)
//
// --- Fiscal rule: interest income tax ---
// rho_ti   = AR(1) persistence
//
// --- Fiscal rule: dividend tax ---
// rho_td   = AR(1) persistence
//
// --- Fiscal rule: expensing rate ---
// rho_etau   = AR(1) persistence (level-deviation form)
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
    PI XI profK
    rev_c rev_w rev_k rev_i rev_d spend_ic rev_net
    // Measurement variables (match data column names)
    chat ihat ghat pi_obs r_obs tau_k_hat tau_l_hat bhat
    trhat tau_c_hat tau_d_hat tau_itc_hat e_tau_hat tau_i_hat
;

varexo
    eps_g eps_tr
    eps_k eps_w
    eps_z eps_d eps_mu
    eps_m
    eps_c eps_ic eps_td   // eps_ti removed 2026-10-22 (interest tax exogenous)
    eps_etau eps_b
;

parameters
    beta delta theta eta etaw gamma
    kappa sig_u
    alphaP alphaW
    bc thetac bg thetag
    tauk_ss tauw_ss pi_ss h_ss g_y b_y
    alphaR alphapi alphaY
    rho_z rho_d rho_mu
    rho_g
    rho_tr rho_tr_b
    rho_k
    rho_w
    phi
    tauc_ss  rho_c
    tauic_ss rho_ic
    taui_ss  rho_ti
    taud_ss  rho_td
    etau_ss  rho_etau
    delta_tau divshare
    epsK divKG
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
g_y      = 0.2073;     // 0.2053 (1958Q1-2019Q4);  0.2073 (1961Q1-2013Q4)
b_y      = 0.3220;     // 0.3549 (1958Q1-2019Q4);  0.3220 (1961Q1-2013Q4)
tauw_ss  = 0.217283;     // 0.216196 (1958Q1-2019Q4);  0.217283 (1961Q1-2013Q4)
tauk_ss  = 0.419945;     // 0.416325 (1958Q1-2019Q4);  0.419945 (1961Q1-2013Q4)

// Extension: SS tax rates (calibrated)
delta_tau = 0.025;     
tauc_ss  = 0.058625;     // 0.058396 (1958Q1-2019Q4); 0.058625 (1961Q1-2013Q4)
tauic_ss = 0.017568;     // 0.015018 (1958Q1-2019Q4); 0.017568 (1961Q1-2013Q4)
etau_ss  = 0.055903;     // 0.088618 (1958Q1-2019Q4); 0.055903 (1961Q1-2013Q4)
taui_ss  = 0.249417;         // 0.20 calibration;         0.249417 (TAXSIM observed years / 1961Q1-2013Q4)
taud_ss  = 0.299868;     // 0.138509 (1958Q1-2019Q4); 0.299868 (taxsim 1961Q1-2013Q4)

// Steady-state dividend share. Targets the postwar average of NIPA
// net corporate dividends / GDP, which centers around 0.05. This is
// the natural empirical counterpart to the model's div = y - w*h - rk*u*k(-1),
// the residual flow accruing to firm owners.
divshare = 0.05;

// Capital-goods sector (ported from the markup paper, its eq. 9).
// epsK is the supply elasticity of capital goods: House and Shapiro (2008).
// FK is a fixed RESOURCE cost, set below so that the capital-goods sector
// earns zero quasi-rent in steady state AND the total goods absorbed by
// investment, XI + FK, equals I exactly. The steady state is therefore
// IDENTICAL to the no-capital-goods version; the sector only bites out of
// steady state, where PI moves with investment.
divKG = 1;         // 1: capital-goods rents are taxed as dividends.
                   // 0: they reach households untaxed. profK_ss = 0 so
                   // the steady state is identical either way; this only
                   // moves the tau_d tax base out of steady state.
epsK = 4.0;        // Capital-goods supply elasticity.
                   // DETERMINACY FLOOR, measured 2026-10-18: the model is
                   // indeterminate for epsK <= 3 and determinate for
                   // epsK >= 4 across kappa in [0.10, 3.00]. epsK = 2 fails
                   // at EVERY kappa, and taxing vs not taxing the rents
                   // (divKG) makes no difference, so the bound comes from
                   // the friction itself. 4.0 is the smallest safe value
                   // and implies a 2.4% rise in the relative price of
                   // investment goods for a 10% rise in investment.
                   // See sweep_epsK.m / sweep_grid.m / sweep_div.m.


// ---------------------------------------------------------------------
// CALIBRATED 2026-10-21: thetac, sig_u, alphaW, bg. IDENTICAL values in
// every model of the 2x2 (new2_cal, new_cal, kgoods2, kgoods-feedback).
// Why calibrate: in tax_dsge_new2 (750k run) posterior sd / prior sd is
//   thetac 1.16, alphaW 1.03 (data silent), sig_u 0.88, bg 0.81 (weak).
// What depends on them (tax_dsge_new2_ridge, 2026-10-21, H20 pv-own):
//   alphaW: multipliers flat over 3 to 120.
//   thetac: ITC and consumption tax move ~10% over its posterior range.
//   sig_u : drives the capital tax multiplier (4.64 at 1.5, 3.09 at 2.5,
//           2.39 at 3.5). Report capital tax at sig_u = 2.0 / 2.5 / 3.0.
//   bg    : drives the spending multiplier (0.54 at 0.55, 0.77 at 0.85).
//           Report spending at bg = 0.55 / 0.70 / 0.85.
// The estimated_params lines are kept below, commented out, so this is
// reversible.
// ---------------------------------------------------------------------
thetac  = 0.53;     // Zubairy (2014, Table 1) posterior mean; data silent here.
                    // new2 750k median 0.479 [0.274, 0.649].
sig_u   = 2.5;      // Conventional curvature; Zubairy prior mean 2.5.
                    // new2 750k median 2.306 [1.630, 3.105].
alphaW  = 100;      // ~4-quarter wage duration (Rotemberg-Calvo mapping, as
                    // in Zubairy's prior). new2 750k median 41.1 [2.99, 97.1].
bg      = 0.70;     // Literature prior mean (Zubairy 2014, Table 1).
                    // new2 750k median 0.715 [0.565, 0.823].


// =====================================================================
// Estimated parameters: initial values
// =====================================================================

gamma   = 2.49;
kappa   = 1.99;
alphaP  = 48.90;
bc      = 0.96;
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
rho_tr_b = 0.439;

// Extension fiscal rule initial values
rho_c    = 0.90;
rho_td   = 0.90;
rho_ic   = 0.90;
rho_ti   = 0.90;
rho_etau = 0.90;


// =====================================================================
// MODEL (49 structural + 14 measurement = 63 equations)
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
    // Capital-goods block. PI_ss is normalised to 1, so the marginal
    // resource cost of investment at the steady state is also 1:
    //   XI_ss = I_ss/(1 + 1/epsK),  FK_val = I_ss/(1 + epsK),
    //   XI_ss + FK_val = I_ss  and  profK_ss = 0.
    #XI_ss   = I_ss/(1 + 1/epsK);
    #FK_val  = I_ss/(1 + epsK);


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
        = c + g + XI + FK_val
        + (gamma1/(1 + sig_u))*(u^(1 + sig_u) - 1)*k(-1)
        + (alphaP/2)*(pi - pi_ss)^2
        + (alphaW/2)*(piw - pi_ss)^2*w;

    // (11) Capital accumulation
    k = (1 - delta)*k(-1) + i*(1 - (kappa/2)*(mu*i/i(-1) - 1)^2);

    // (11b-11d) Capital-goods sector [NEW]
    //   Convex transformation cost gives a finite supply elasticity epsK and
    //   a relative price of investment that rises with the quantity produced.
    //   The producer takes PI as given and maximises profK, so PI = XI'(i).
    PI    = (i/I_ss)^(1/epsK);
    XI    = XI_ss*(i/I_ss)^(1 + 1/epsK);
    profK = PI*i - XI - FK_val;

    // (12) Tax basis of capital [MODIFIED: accumulates the PURCHASE PRICE]
    k_tau = (1 - delta_tau)*k_tau(-1) + PI*i*(1 - e_tau);

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
    lambda*PI*(1 - tau_ic - tau_k*e_tau) = lambda*q * (
        1 - (kappa/2)*(mu*i/i(-1) - 1)^2
          - (mu*i/i(-1)) * kappa*(mu*i/i(-1) - 1)
    )
    + beta*lambda(+1)*q(+1) * (
        mu(+1)*(i(+1)/i)^2 * kappa*(mu(+1)*i(+1)/i - 1)
    )
    + lambda*v_tau*(1 - e_tau)*PI;

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
    div = y - w*h - rk*u*k(-1) + divKG*profK;

    // (24-25) Government budget and tax revenue [MODIFIED]
    b = R(-1)*b(-1)/pi + g + tr - taxrev + eps_b;

    taxrev = tau_c*c
           + tau_w*w*h
           + tau_k*(rk*u*k(-1) - delta_tau*k_tau(-1) - e_tau*PI*i)
           + tau_i*(R(-1) - 1)*b(-1)/pi
           + tau_d*div
           - tau_ic*PI*i;

    // (26-29) Fiscal rules: pure AR(1), no debt or output feedback
    log(tau_k/tauk_ss) = rho_k*log(tau_k(-1)/tauk_ss) + eps_k;

    log(tau_w/tauw_ss) = rho_w*log(tau_w(-1)/tauw_ss) + phi*eps_k + eps_w;

    log(g/G_ss) = rho_g*log(g(-1)/G_ss) + eps_g;

    // Transfers are the ONLY instrument retaining debt feedback. The
    // negative sign is what stabilizes: transfers fall when debt is above
    // steady state. rho_tr_b > 0 (gamma prior).
    log(tr/tr_ss) = rho_tr*log(tr(-1)/tr_ss)
        - rho_tr_b*log(b(-1)/b_ss) + eps_tr;

    // (30) Taylor rule (unchanged)
    log(R/R_ss) = alphaR*log(R(-1)/R_ss)
        + (1 - alphaR)*(alphapi*log(pi/pi_ss) + alphaY*log(y/Y_ss)) + eps_m;

    // (31-34) Extension fiscal rules [EXTENSION]
    log(tau_c/tauc_ss) = rho_c*log(tau_c(-1)/tauc_ss) + eps_c;

    tau_ic - tauic_ss = rho_ic*(tau_ic(-1) - tauic_ss) + eps_ic;

    // INTEREST INCOME TAX EXOGENOUS (constant), 2026-10-22. Was:

    //   

    //   log(tau_i/taui_ss) = rho_ti*log(tau_i(-1)/taui_ss) + eps_ti;

    tau_i = taui_ss;

    log(tau_d/taud_ss) = rho_td*log(tau_d(-1)/taud_ss) + eps_td;

    // (35) Expensing rate fiscal rule [NEW: e_tau promoted to stochastic]
    //   Level-deviation form (like tau_ic) since SS may be zero.
    e_tau - etau_ss = rho_etau*(e_tau(-1) - etau_ss) + eps_etau;

    // (36-38) Exogenous processes (unchanged)
    log(z)  = rho_z*log(z(-1))   + eps_z;
    log(d)  = rho_d*log(d(-1))   + eps_d;
    log(mu) = rho_mu*log(mu(-1)) + eps_mu;

    // (39-45) Fiscal revenue accounting [EXTENSION]
    rev_c    = tau_c*c;
    rev_w    = tau_w*w*h;
    rev_k    = tau_k*(rk*u*k(-1) - delta_tau*k_tau(-1) - e_tau*PI*i);
    rev_i    = tau_i*(R(-1) - 1)*b(-1)/pi;
    rev_d    = tau_d*div;
    spend_ic = tau_ic*PI*i;
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
    tau_i_hat   = (tau_i - taui_ss) / taui_ss;
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
    XI_ss_v    = I_ss_v/(1 + 1/epsK);
    FK_val_v   = I_ss_v/(1 + epsK);
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
    PI    = 1;
    XI    = XI_ss_v;
    profK = 0;
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
    tau_i_hat   = 0;
    tau_itc_hat = 0;
    e_tau_hat   = 0;

end;

steady;
check;


// =====================================================================
// OBSERVABLES — 12 observables, 14 shocks (Run 1: data trim)
// =====================================================================
// Dropped: trhat (2008-2010 ARRA outliers)
//          tau_i_hat (0.77 collinearity with tau_d_hat)
// Kept all 14 shocks (eps_tr and eps_ti become unobserved structural shocks).
// eps_b is included as a government budget / debt residual shock.
// =====================================================================
varobs chat ihat ghat bhat        // trhat
       pi_obs r_obs 
       tau_k_hat tau_l_hat
       tau_c_hat tau_itc_hat e_tau_hat tau_d_hat;   // tau_i_hat  

// =====================================================================
// PRIORS
// =====================================================================
// Fiscal rules are pure AR(1) processes. The debt and output feedback
// terms of the Zubairy specification have been dropped from every rule
// except transfers, which retains debt feedback as the sole anchor on
// government debt. The Zubairy-style specification is archived in
// ../tax_dsge_new.
//
// Persistences estimated for all nine instruments:
//   rho_k, rho_w, rho_g, rho_tr, rho_c, rho_ic, rho_ti, rho_td, rho_etau
// Plus rho_tr_b (transfer debt feedback) and the cross-shock term phi.
// All 14 shock stds, including eps_b.
// Total: 40 estimated parameters (was 57).
// =====================================================================
estimated_params;

    // --- Deep habits ---
    bc,       beta_pdf,       0.7,   0.1;
    // thetac,   beta_pdf,       0.8,   0.1;   // CALIBRATED 2026-10-21
    // bg,       beta_pdf,       0.7,   0.1;   // CALIBRATED 2026-10-21
    thetag,   beta_pdf,       0.8,   0.1;

    // --- Structural ---
    gamma,    normal_pdf,     2.0,   1.0;
    // sig_u,    normal_pdf,     2.5,   0.5;   // CALIBRATED 2026-10-21
    kappa,    normal_pdf,     2.0,   0.5;

    // --- Nominal rigidities ---
    alphaP,   normal_pdf,     17.0,  5.0;
    // alphaW,   normal_pdf,     100.0, 30.0;  // CALIBRATED 2026-10-21
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

    // --- Fiscal rule: debt response ---
    // Transfers carry the only debt feedback; every other rule is a pure
    // AR(1). NOTE: with transfers as the sole anchor the (b, tr) block is
    // stable only for rho_tr_b in roughly [0.005, 0.55] at rho_tr = 0.96,
    // and the upper bound tightens as rho_tr rises. This gamma prior puts
    // about a third of its mass above that bound, so expect some
    // Blanchard-Kahn rejections. Tighten toward 0.15 if it degrades
    // acceptance.
    rho_tr_b, gamma_pdf,      0.5,   0.25;

    // --- Tax correlations ---
    phi,      normal_pdf,     0.25,  0.1;

    // --- AR coefficients (extension) ---
    rho_c,    beta_pdf,       0.7,   0.2;
    rho_ic,   beta_pdf,       0.7,   0.2;
    // rho_ti,   beta_pdf,       0.7,   0.2;   // interest tax exogenous, 2026-10-22
    rho_td,   beta_pdf,       0.7,   0.2;
    rho_etau, beta_pdf,       0.7,   0.2;

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
    // stderr eps_ti,   inv_gamma_pdf, 0.05,  0.1;   // interest tax exogenous, 2026-10-22
    stderr eps_td,   inv_gamma_pdf, 0.05,  0.1;
    stderr eps_etau, inv_gamma_pdf, 0.05,  0.1;
    stderr eps_b, inv_gamma_pdf, 0.05, 0.1;

end;


// =====================================================================
// ESTIMATION
// =====================================================================

// ---------------------------------------------------------------------
// BLOCK 1: full run. This is the live block.
// Reuses the mode and proposal covariance from mode_compute = 6, so it
// skips the slow mode-finding step.
// ---------------------------------------------------------------------

// mh_replic = 3M. This is what the live run of 2026-10-14 19:35 UTC used,
// so the value is kept here to match the chains on disk. At the pilot's
// median tau of 3928 it gives median ESS ~1500 and p90 ~840.
// 1.5M (Zubairy's count) would have given median ~760 / p90 ~420, also
// adequate; 3M was chosen as insurance and costs only wall clock.
// If diagnostics ever come back thin, chains can be EXTENDED rather than
// rerun: add load_mh_file and raise mh_replic, and do NOT delete
// tax_dsge_est/metropolis in that case.

// DISABLED 2026-10-18 for the determinacy test. Re-enable after check passes.
/*
estimation(
    datafile     = US_Data_Matlab_1961_2013,
    prefilter    = 1,
    lik_init     = 1,
    mh_replic    = 750000,
    mh_nblocks   = 4,
    mh_drop      = 0.5,
    mh_jscale    = 0.20,
    mode_compute = 0,
    mode_file    = tax_dsge_est_mode_from_chains,
    bayesian_irf,
    irf          = 40,
    nodisplay
) y c i h pi R g tr div taxrev
  rev_c rev_w rev_k rev_d rev_i spend_ic;
*/


 
/*
// ---------------------------------------------------------------------
// BLOCK 2: reload existing chains and mode. Samples nothing.
// Use this before compute_multipliers, compute_multipliers_reveq,
// posterior_table, mcmc_diagnostics, check_chains, sensitivity_bc.
// mh_nblocks and mh_drop must match the run that produced the chains.
// ---------------------------------------------------------------------
estimation(
    datafile     = US_Data_Matlab_1961_2013,
    prefilter    = 1,
    lik_init     = 1,
    mh_replic    = 0,                 // skip MCMC
    load_mh_file,                     // load existing chains
    mode_compute = 0,                 // skip mode-finding
    mode_file    = tax_dsge_est_mode,
    mh_nblocks   = 4,
    mh_drop      = 0.5,
    mh_jscale    = 0.05,
    bayesian_irf,
    irf          = 40,
    nodisplay
) y c i h pi R g tr div taxrev
  rev_c rev_w rev_k rev_d rev_i spend_ic;
*/
 
 
// ---------------------------------------------------------------------
// BLOCK 2k: RELOAD the converged chains. Samples nothing.  <-- LIVE
//
// Run this before compute_multipliers, compute_multipliers_reveq,
// posterior_table, mcmc_diagnostics, check_chains, sensitivity_bc.
//
// compute_multipliers reads the metropolis folder itself and builds each
// draw's IRF with resol + simult_, so bayesian_irf is NOT needed and is
// deliberately omitted: it cost ~7 hours of post-processing on new2.
//
// CONVERGED, 2026-10-20, 1M x 4 draws, 7h27m, acceptance 0.54:
//   Rhat        max 1.0094 (rho_ti), every parameter under 1.05
//   drift       +0.6 log points over the retained half
//   block spread 0.7
//   ESS         median 2061, min 902
//   tau         median 333, p90 516
//   block maxima 3335.81 / 3336.23 / 3337.14 / 3335.64
//   block means  3323.13 / 3322.86 / 3322.80 / 3322.80
// mh_nblocks and mh_drop must match the run that produced the chains.
// ---------------------------------------------------------------------
/* DISABLED in tax_dsge_kgoods2_tifix until its estimation run is configured.
estimation(
    datafile     = US_Data_Matlab_1961_2013,
    prefilter    = 1,
    lik_init     = 1,
    mh_replic    = 0,
    load_mh_file,
    mh_nblocks   = 4,
    mh_drop      = 0.5,
    mh_jscale    = 0.15455,
    mode_compute = 0,
    mode_file    = tax_dsge_est_mode36,
    nodisplay
) y c i h pi R g tr div taxrev
  rev_c rev_w rev_k rev_d rev_i spend_ic;
*/


// ---------------------------------------------------------------------
// BLOCK 1k: PRODUCTION RUN, capital-goods model.  <-- LIVE
//
// Decision 2026-10-19: do NOT change the specification. The kgoods posterior
// has a genuine ridge that tax_dsge_new2 does not (within-block correlations
// on drift-free windows: median 0.20 and ~10 pairs above 0.7, against 0.07
// and 2 pairs for new2). Hub is sig_u <-> alphaY at -0.75, with an
// eps_ti / bg / thetag cluster. Economic reading: with capital-goods supply
// inelastic, investment cannot adjust freely, so capital UTILIZATION becomes
// the marginal margin and sig_u does work it did not do before.
//
// mode_compute = 5 FAILED here (min Hessian eigenvalue -9.4e7, not positive
// definite) because a ridge has no well-defined peak. So the mode and the
// proposal both come from the chains, as on new2.
//
// By this point the chains had essentially arrived in LOCATION: drift +8.9
// log points over the retained half, block spread 4.7. What remains is slow
// mixing ALONG the ridge, and the answer to that is draws.
//
// 2M x 4 at ~385 draws/sec ~= 6 hours. NO bayesian_irf: it added ~7 hours of
// post-processing on new2 and IRFs regenerate from saved chains afterwards.
//
// Judge this run by whether the MULTIPLIERS are stable, not by whether every
// nuisance parameter mixes. bg, thetag and sig_u are already reported as
// prior-driven in new2.
// ---------------------------------------------------------------------
/*
estimation(
    datafile     = US_Data_Matlab_1961_2013,
    prefilter    = 1,
    lik_init     = 1,
    mh_replic      = 1000000,
    mh_nblocks     = 4,
    mh_drop        = 0.5,
    mh_jscale      = 0.15455,   // tuned value from the 300k pilot,
                                // acceptance 0.29 across all four blocks
    mode_compute   = 0,
    mode_file      = tax_dsge_est_mode36,
    nodisplay
) y c i h pi R g tr div taxrev
  rev_c rev_w rev_k rev_d rev_i spend_ic;
*/


// ---------------------------------------------------------------------
// BLOCK 3c: mode REFINEMENT from the chain argmax.  <-- LIVE
//
// WHY. make_mode_from_chains does two jobs and they are not equally valid
// while a chain is still travelling:
//   - the ARGMAX is fine: a good point is a good point.
//   - the COVARIANCE is not: a drifting chain's covariance is dominated by
//     the direction of travel, hugely inflated along it. Feeding that back
//     as a proposal made mixing WORSE, tau 2231 -> 6965, and forced the
//     tuner down to jscale 0.0316.
//
// So: take the mode from the chains, take the covariance from a Hessian
// computed AT that mode. The original Hessian was not bad because Hessians
// are bad; it was bad because it sat 400 log points below the truth.
//
// mode_file supplies the STARTING VALUES (chain argmax, log posterior
// 3101.96); mode_compute = 5 (newrat) optimises from there and writes a
// fresh Hessian. mh_replic = 0, so no sampling and no stuck-chain illusion.
// ---------------------------------------------------------------------
/*
estimation(
    datafile     = US_Data_Matlab_1961_2013,
    prefilter    = 1,
    lik_init     = 1,
    mh_replic    = 0,
    mh_nblocks   = 4,
    mh_drop      = 0.5,
    mode_compute = 5,
    mode_file    = tax_dsge_est_mode_from_chains,
    nodisplay
) y c i h pi R g tr div taxrev
  rev_c rev_w rev_k rev_d rev_i spend_ic;
*/


// ---------------------------------------------------------------------
// BLOCK 3b: jscale pilot on the CAPITAL-GOODS model.  <-- LIVE
// Mode from BLOCK 3 (mode_compute = 6, 2026-10-18, 21m45s,
// log posterior 2700.87, Hessian positive definite).
//
// Purpose is twofold:
//   (a) tune mh_jscale, and
//   (b) give make_mode_from_chains something to work with. On the
//       tax_dsge_new2 model, mode_compute = 6 missed the mode by 441 log
//       points, so the Hessian proposal described curvature the chains
//       never visited. Expect the same here and plan to rebuild.
//
// mh_drop = 0 keeps every draw for the autocorrelation estimate.
// Start at mh_jscale = 0.06: that is what gave acceptance 0.23 on the
// Hessian-derived proposal in tax_dsge_new2. LOWER jscale RAISES acceptance.
//
// After it finishes:  check_chains  then  mcmc_diagnostics
// Delete tax_dsge_est/metropolis before re-running this at another scale.
// ---------------------------------------------------------------------
/*
estimation(
    datafile     = US_Data_Matlab_1961_2013,
    prefilter    = 1,
    lik_init     = 1,
    mh_replic       = 100000,
    mh_nblocks      = 4,
    mh_drop         = 0,
    mh_tune_jscale  = 0.25,   // auto-tune the scale.
                              // ROUND 1, Hessian proposal from mode_compute=6:
                              //   jscale 0.06  -> accept 0.80, median ESS 32
                              //   jscale 0.387 -> accept 0.24, median ESS 38
                              // Scale was never the real problem: the Hessian
                              // had the SHAPE wrong. Chain covariance needs
                              // steps 12x wider in some directions (eps_ti,
                              // rho_ti, alphaY, alphaR) and 5x narrower in
                              // others (eps_c, eps_z), which no scalar jscale
                              // can deliver.
                              // ROUND 2: proposal rebuilt from chains, mode
                              // 2720.78. Chains then took off: best 3016.90,
                              // +296, still climbing at the end. The earlier
                              // "flat" runs were stuck chains, not converged
                              // ones.
                              // ROUND 3 (this one): rebuilt again from the top
                              // 20% of those chains. 100k draws for more travel
                              // distance per round.
    mode_compute    = 0,
    mode_file       = tax_dsge_est_mode_from_chains,
    nodisplay
) y c i h pi R g tr div taxrev
  rev_c rev_w rev_k rev_d rev_i spend_ic;
*/


// ---------------------------------------------------------------------
// BLOCK 3: mode-finding from scratch. Done 2026-10-14, 19m53s.
// Log posterior 2862.54, Hessian positive definite, 40 parameters.
// Required because the fiscal rules changed: the old 57-parameter mode
// file does not apply to this 40-parameter model.
// mode_compute = 6 is slow, on the order of a day. mh_replic = 0 so it
// finds the mode and stops without sampling; the proposal scale is tuned
// separately afterwards with BLOCK 4, because mode_compute = 6 overrides
// mh_jscale with its own value.
// ---------------------------------------------------------------------
/*
estimation(
    datafile     = US_Data_Matlab_1961_2013,
    prefilter    = 1,
    lik_init     = 1,
    mh_replic    = 0,
    mh_nblocks   = 4,
    mh_drop      = 0.5,
    mode_compute = 6,   // 4 is csminwel, 6 is Monte Carlo based
    nodisplay
) y c i h pi R g tr div taxrev
  rev_c rev_w rev_k rev_d rev_i spend_ic;
*/




// ---------------------------------------------------------------------
// BLOCK 4: jscale pilot. Done. Results on the 40-parameter model:
//   jscale 0.04 -> acceptance 0.33, median ESS 39, median tau 5183
//   jscale 0.06 -> acceptance 0.23, median ESS 51, median tau 3928  <- adopted
// Acceptance is now at the 0.234 optimum; scale tuning is exhausted. Run after BLOCK 3 writes a new mode file.
// 50k x 4 blocks, mh_drop = 0 to keep every draw for the autocorrelation
// estimate. Read realized acceptance from check_chains and the implied
// mh_replic from mcmc_diagnostics, then set BLOCK 1 and run it once.
// Target acceptance ~0.234. LOWER mh_jscale to RAISE acceptance.
// Delete tax_dsge_est/metropolis before this and before BLOCK 1.
// ---------------------------------------------------------------------
/*
estimation(
    datafile     = US_Data_Matlab_1961_2013,
    prefilter    = 1,
    lik_init     = 1,
    mh_replic    = 50000,
    mh_nblocks   = 4,
    mh_drop      = 0,
    mh_jscale    = 0.20,   // Scale sweep on the chain-covariance proposal.
                           // 0.38 -> acceptance 0.038, median ESS 203,
                           //         median tau 991, zero block gaps.
                           // Acceptance is NOT the target: the posterior is
                           // far from Gaussian so 0.234 does not apply here.
                           // Compare median ESS across scales instead.
    mode_compute = 0,
    mode_file    = tax_dsge_est_mode_from_chains,
    nodisplay
) y c i h pi R g tr div taxrev
  rev_c rev_w rev_k rev_d rev_i spend_ic;
*/

 

// =====================================================================
// RUN BLOCKS for tax_dsge_kgoods2_tifix, 2026-10-21. The stage is set on the command line
// by overnight_2026_10_21.m; with no -DSTAGE nothing is estimated.
//   -DSTAGE=1  mode search (newrat) starting from start_mode.mat
//              (built from the old chains: best draw, chain covariance).
//   -DSTAGE=2  MCMC 4 x 1M from tax_dsge_est_mode_cal.mat = stage-1 mode
//              and fval (so C4 is valid) with the chain-covariance
//              proposal that gave tau 333 on kgoods.
// =====================================================================
@#ifndef STAGE
  @#define STAGE = 0
@#endif

@#if STAGE == 1
estimation(
    datafile     = US_Data_Matlab_1961_2013,
    prefilter    = 1,
    lik_init     = 1,
    mh_nblocks   = 4,
    mh_drop      = 0.5,
    mh_replic    = 0,
    mode_compute = 5,
    optim        = ('MaxIter', 100),
    mode_file    = start_mode,
    nograph,
    nodisplay
) y c i h pi R g tr div taxrev
  rev_c rev_w rev_k rev_d rev_i spend_ic;
@#endif

@#if STAGE == 2
estimation(
    datafile     = US_Data_Matlab_1961_2013,
    prefilter    = 1,
    lik_init     = 1,
    mh_nblocks   = 4,
    mh_drop      = 0.5,
    mh_replic    = 1000000,
    mh_jscale    = 0.15455,
    mode_compute = 0,
    mode_file    = tax_dsge_est_mode_cal,
    nograph,
    nodisplay
) y c i h pi R g tr div taxrev
  rev_c rev_w rev_k rev_d rev_i spend_ic;
@#endif

// ---------------------------------------------------------------------
// STAGE 5 (2026-10-23): RELOAD the converged chains, sample nothing.
// Run before compute_multipliers and post_estimation.
// ---------------------------------------------------------------------
@#if STAGE == 5
estimation(
    datafile     = US_Data_Matlab_1961_2013,
    prefilter    = 1,
    lik_init     = 1,
    mh_nblocks   = 4,
    mh_drop      = 0.5,
    mh_replic    = 0,
    load_mh_file,
    mode_compute = 0,
    mode_file    = tax_dsge_est_mode_cal,
    nograph,
    nodisplay
) y c i h pi R g tr div taxrev
  rev_c rev_w rev_k rev_d rev_i spend_ic;
@#endif
