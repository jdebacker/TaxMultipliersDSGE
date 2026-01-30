// /************ tax_dsge.mod **************/
// /**  This program sets up the model to be solved by Dynare.
// ***  There are five part to this file:
// ***  1. The preamble, which initializes variables and parameters.
// ***  2. Model, which spells out the model (necessary equations).
// ***  3. Steady State, which gives the steady state values.
// ***  4. Shocks, which defines the shocks.
// ***  5. Computation, which specifies solution method and output.
// ***
// ***  Throughout, I try to use the notation we use in our
// ***  write-up of the model.
// ***/

/* 
*addpath C:\dynare\5.5\matlab // 
*dynare tax_dsge.mod  
*/

// /***************************************************************/
// /* Preamble */
// /***************************************************************/
// /* Endogenous Variables (22) */
// /** Note that there are more variables here than in the paper - this is because it makes the coding easier to define some variables that are functions of others (e.g. marginal cost) **/
// c = consumption (private)
// lambda = marginal utility of consumption
// q = lagrangian multiplier on the law of motion for capital
// l = labor supply
// i = investment
// b = gov't bond holdings
// v = capital utilization rate
// k = capital stock
// k_tau = tax basis of capital stock
// g = gov't consumption
// x = gov't transfers
// d = dividends distributed from intermediate goods producers
// z = total factor productivity (AR(1) process) ** Note that z here is ln(z) from the paper **
// y = output
// p = price level
// w = wage rate
// r_k = capital rental rate
// r = interest rate on gov't bonds
// lambda_p = price markup (AR(1) process) ** Note that lambda_p here is ln(lambda_p) from the paper **
// epsilon_b = stochastic time preference (AR(1) process) ** Note that epsilon_b here is ln(epsilon_b) from the paper **
// p_star = the price chosen by those int goods producers who can change price
// mc = marginal cost of the intermediate goods producer
//
/* Exogenous Variables (3) */
// e_b = shock to stochastic time preference
// e_p = shock to price markup
// e_z = shock to TFP
//
/*  Parameters (30) */
// betta = time preference parameter
// siggma = coefficient of relative risk aversion for utility function
// zetta = Frisch elasticity of labor supply??
// siggma_g = preference parameter for gov't spending (some elasticity)
// chi_g = scale parameter for preferences for gov't spending
// rho_b = persistence of shock to time preference
// gamma = scale paramter in investment adjustment cost function
// delta_0 = SS depreciation rate
// delta_1 = coefficient on linear term in depreciation funciton
// delta_2 = coefficient on quadratic term in depreciation function
// alfa = capital's share of output
// lambda_bar_p = average price markup (BETTER NOTATION!!!!)
// theta = fraction of intermediate goods producers who cannot change price
// rho_z = persistence of shock to TFP
// siggma_z = std dev of shocks to TFP
// siggma_p = std dev of shocks to price markup
// siggma_b = std dev of shocks to time preference
//
// phi_1 = coeff on inflation rate in Taylor rule
// phi_2 = coeff on output gap in Taylor rule
// rho_r = coeff on interest rate difference in Taylor rule
//
// tau_c = consumption tax
// tau_i = tax on interest income
// tau_k = tax on capital income
// tau_l = tax on labor income
// tau_d = tax on dividend income
// tau_ic = investment tax credit
// e_tau = rate of expensing for tax purposes
// delta_tau = rate of depreciation for tax purposes
//
// gamma_x = percentage of gov't spending going to transfers - make this a parameter to help id x and g separately
//
// lbar = SS labor supply (I think we need this to ID model in SS)
/***************************************************************/
var p p_star v z epsilon_b lambda_p r r_k w l k k_tau i y d c lambda q g b x mc ;

varexo e_b e_p e_z ;

parameters betta siggma zetta siggma_g chi_g rho_b gamma delta_0 delta_1 delta_2 alfa lambda_bar_p theta phi_1 phi_2 rho_r rho_z rho_lambda siggma_z siggma_b siggma_p tau_c tau_i tau_k tau_l tau_d tau_ic e_tau delta_tau gamma_x lbar ybar rbar r_k_bar wbar mc_bar ;

betta = 0.95 ;
siggma = 2 ;
zetta = 0.29 ;
siggma_g = 1.1 ;
chi_g = 0.5 ;
rho_b = 0.9 ;
gamma = 2 ;
delta_0 = 0.1 ;
delta_1 = 0.152632 ; //get this from SS, delta_1 = 1/betta-1+delta_0
delta_2 = 0.01 ;
alfa = 0.33 ;
lambda_bar_p = 7 ;
theta = 0.82 ;
rho_z = 0.7 ;
rho_lambda = 0 ;
siggma_z = 0.01 ;
siggma_p = 0.01 ;
siggma_b = 0.01 ;

phi_1 = 0.2 ;
phi_2 = 0.03 ;
rho_r = 0.1 ;

tau_c = 0.1 ;
tau_i = 0.2 ;
tau_k = 0.35 ;
tau_l = 0.25 ;
tau_d = 0.15 ;
tau_ic = 0 ;
e_tau = 0 ;
delta_tau = 0.1 ;
//0.13 ;

gamma_x = 0.2 ;

lbar = 0.5 ;

// Create steady state variables needed for model equations
rbar = (1-betta)/(betta*(1-tau_i)) ;
r_k_bar = (((1/betta)-1+delta_0)/(1-tau_k))*(1-tau_ic -(tau_k*e_tau) - (betta*tau_k*delta_tau*(1-e_tau))) ;
wbar = (1-alfa)*((r_k_bar/alfa)^((-1*alfa)/(1-alfa))) ;
ybar = (((alfa/(1-alfa))*(wbar/r_k_bar))^(alfa))*lbar ;
mc_bar = ((wbar^(1-alfa))*(r_k_bar^alfa))/((alfa^alfa)*((1-alfa)^(1-alfa))) ;



/***************************************************************/
/* Model */
/***************************************************************/
// Equations that define the model:
// 1. HH Budget Constraint -> consumption
// 2. HH FOC, consumption -> lambda (Marg Util Cons)
// 3. HH FOC, labor supply -> labor supply
// 4. HH FOC, capital -> investment
// 4. HH FOC, investment -> q
// 5. HH FOC, bond holdings -> demand for gov't bonds
// 6. HH FOC, capital utilization -> capital utilization
// 7. Law of motion for capital stock -> capital stock
// 8. Law of motion for tax basis of capital stock -> tax basis of cap stock
// 9. Resource constraint -> y
// 10. Calibration of g/y -> g
// 11. Gov't budget constraint -> x
// 12. Taylor Rule for monetary authority -> r
// 13. Int. goods producer FOC, effective capital demand -> (with market clearning condition) r_k
// 14. Int. goods producer FOC, labor demand -> (with market clearing condition) -> w
// 15. Int. goods producr FOCs -> mc (intermediate variable)
// 16. Int. goods producer FOC, price -> p_star
// 17. Int. goods producer profit function -> d
// 18. Calvo pricing rule -> p
// 19. AR(1) process for TFP -> z
// 20. Price markup process -> lambda_p
// 21. AR(1) process for stochastic time preference -> epsilon_b
/***************************************************************/
model ;

// 1. HH Budget Constraint (Eq. 1.2)
c = ((((1+(r*(1-tau_i)))*b(-1))/p) +(((1-tau_l)*w*l)/p) + (((1-tau_k)*r_k*v*k(-1))/p) + (tau_k*delta_tau*k_tau(-1)) + (tau_ic*i) + (tau_k*e_tau*i) + (((1-tau_d)*d)/p) + x - i - (b/p))/(1+tau_c) ;

// 2. HH FOC, consumption (Eq. 1.25)
lambda = exp(epsilon_b)*(zetta*(c^(zetta*(1-siggma)-1))*((1-l)^((1-zetta)*(1-siggma))))/(1+tau_c) ;

// 3. HH FOC, labor supply (Eq. 1.26)
exp(epsilon_b)*((1-zetta)*(c^(zetta*(1-siggma)))*((1-l)^(((1-zetta)*(1-siggma))-1))) = (lambda*w*(1-tau_l))/p ;

// 4. HH FOC, capital (Eq. 1.27)
q = betta*(((lambda(+1)*(1-tau_k)*r_k(+1)*v(+1))/p(+1)) + (q(+1)*(1-delta_0-(delta_1*(v(+1)-1))-((delta_2/2)*((v(+1)-1)^2))))) ;

// 5. HH FOC, investment (Eq. 1.24)
lambda*(1-tau_ic-(tau_k*e_tau)) = (q*(1-((gamma/2)*(((i/i(-1))-1)^2))-(gamma*((i/i(-1))-1)*(i/i(-1))))) + (betta*((lambda(+1)*tau_k*delta_tau*(1-e_tau)) + (q(+1)*gamma*((i(+1)/i)-1)*((i(+1)/i)^2)))) ;

// 6. HH FOC, bond holdings (Eq. 1.28)
lambda/p = ((betta*lambda(+1)*(1+(r(+1)*(1-tau_i))))/p(+1)) ;

// 7. HH FOC, capital utilization (Eq. 1.22)
(lambda*(1-tau_k)*r_k*k(-1))/p = q*(delta_1+(delta_2*(v-1)))*k(-1) ;

// 8. Law of motion for capital stock (Eq. 1.3)
k = ((1-(delta_0+(delta_1*(v-1))+((delta_2/2)*((v-1)^2))))*k(-1)) + (i*(1-((gamma/2)*((i/i(-1))-1)^2))) ;

// 9. Law of motion for tax basis of capital stock (Eq. 1.4)
k_tau = ((1-delta_tau)*k_tau(-1)) + (i*(1-e_tau)) ;

// 10. Resource constraint (Eq. 1.46)
y = c + i + g ;

// 11. Gov't budget constraint (Eq. 1.79)
g = (1-gamma_x)*((tau_c*c) + ((tau_l*l*w)/p) + (b/p) + (tau_k*r_k*v*k(-1)) + ((tau_d*d)/p) - (tau_ic*i) - (tau_k*delta_tau*k_tau(-1)) - (tau_k*e_tau*i) - ((b(-1)*(1+(r*(1-tau_i))))/p)) ;

// 12. Calibrate x is a fraction of gov't spending (Eq. 1.80)
x = (gamma_x)*((tau_c*c) + ((tau_l*l*w)/p) + (b/p) + (tau_k*r_k*v*k(-1)) + ((tau_d*d)/p) - (tau_ic*i) - (tau_k*delta_tau*k_tau(-1)) - (tau_k*e_tau*i) - ((b(-1)*(1+(r*(1-tau_i))))/p)) ;

// 13. Taylor Rule for monetary authority (Eq. 1.82)
// r(+1) = max(((1/betta)*((1+((p-p(-1))/p(-1)))^phi_1)*((y/ybar)^phi_2)*(((1+r)/(1+rbar))^rho_r)-1)/(1-tau_i),0) ;
r(+1) = ((1/betta)*((1+((p-p(-1))/p(-1)))^phi_1)*((y/ybar)^phi_2)*(((1+r)/(1+rbar))^rho_r)-1)/(1-tau_i) ;


// 14. Int. goods producer FOC, effective capital demand (Eq. 1.64)
v*k(-1) = (alfa/r_k)*y*p ;

// 15. Int. goods producer FOC, labor demand (Eq. 1.65)
l = ((1-alfa)/w)*y*p ;

// 16. Int goods producer marginal cost (Eq. T.2.5)
mc = ((w^(1-alfa))*(r_k^alfa))/(p*exp(z)*(alfa^alfa)*((1-alfa)^(1-alfa))) ;

// 17. Int. goods producer FOC, price (Eq. 1.24)
// p_star = ((1 + lambda_bar_p) / (1 - betta * theta)) * mc; // Calvo + CES microfoundations
p_star = (((1+lambda_bar_p)*mc_bar*(1-(betta*theta)))*(((mc-mc_bar)/mc_bar)+((p-1)/1)+((lambda_bar_p/(1+lambda_bar_p))*((lambda_p-lambda_bar_p)/lambda_bar_p)))) + (betta*theta*p_star(+1)) ;
// p_star = ((((1+lambda_bar_p)*mc_bar*(1-(betta*theta)))/(1-((1-lambda_bar_p)*mc_bar*(1-(betta*theta))*(1-theta))))*(((mc-mc_bar)/mc_bar)+((theta*p(-1))-1)+((lambda_bar_p/(1+lambda_bar_p))*((lambda_p-lambda_bar_p)/lambda_bar_p)))) + (betta*theta*p_star(+1)) ;
// the ss state with second equation has larger residuals, if we use this , we need may need to set a proper convergence condition for p ∗(+1) which could involve iterating or adjusting the solution method.

// 18. Int. goods producer profit function (Eq. 1.62)
d = y*((((theta*((p(-1)/p)^((-1*(1+lambda_p))/lambda_p)))*(p(-1)-mc))) + (((1-theta)*((p_star/p)^((-1*(1+lambda_p))/lambda_p)))*(p_star-mc))) ;
// d = y - w * l - r_k * k; used for stable steady state

// 19. Calvo pricing rule (Eq. 1.51)
p = ((theta*(p(-1)^(-1/lambda_p))) + ((1-theta)*(p_star^(-1/lambda_p))))^(-lambda_p) ;

// 20. AR(1) process for TFP (Eq. 1.61)
z = (rho_z*z(-1)) + e_z ;

// 21. Price markup process (Eq. 1.34)
ln(lambda_p) = (rho_lambda*(ln(lambda_p(-1)))) + ((1-rho_lambda)*ln(lambda_bar_p)) + e_p ;

// 22. AR(1) process for stochastic time preference (Eq. 1.33)
epsilon_b = (rho_b*epsilon_b(-1)) + e_b ;


end ;


/***************************************************************/
/* Steady State (symbolic) */
/***************************************************************/

steady_state_model;

// 1. Exogenous variables — steady-state values
z         = 0;                      // Eq. (1.1): TFP shock
epsilon_b = 0;                      // Eq. (1.2): Discount factor shock
lambda_p  = lambda_bar_p;          // Eq. (1.20): Price adjustment cost shock

// 2. Normalizations
p         = 1;                      // Eq. (1.3): Normalize aggregate price level
v         = 1;                      // Normalization: value of intermediate firm

// 3. Prices and policy rates
r         = rbar;                   // Eq. (1.12): Household FOC for bonds
r_k       = r_k_bar;               // Eq. (1.10): Rental rate of capital
w         = wbar;                   // Eq. (1.9): Wage from labor FOC
l         = lbar;                   // Assumed calibration
mc        = mc_bar;                 // Eq. (1.7): Marginal cost from production

// 4. Output and capital
y         = ybar;                   // Eq. (1.4): Output from production function
k         = (alfa / (1 - alfa)) * (w * l / r_k);                       // Eq. (1.10): Capital demand FOC
k_tau     = delta_0 * k * (1 - e_tau) / delta_tau;                    // Eq. (1.18): Tax-depreciation capital stock
i         = delta_0 * k;                                              // Eq. (1.5): Steady-state capital accumulation

// 5. Dividends
d         = y - w * l - r_k * k;   // Eq. (1.6): Profits = Y - WL - RK

// 6. Government spending g — first from budget constraint
g = (                                     // Eq. (1.45): Government budget constraint (used to pin down g)
  tau_c * (y - i)
+ tau_l * w * l / p
+ tau_k * r_k * k
+ tau_d * d / p
- tau_ic * i
- tau_k * delta_tau * k_tau
- tau_k * e_tau * i
) * (1 - gamma_x);

// 7. Calibrate b from Eq. (1.47)
b = (                                     // Eq. (1.47): Intertemporal government budget constraint
  tau_c * (y - i)
+ tau_l * w * l / p
+ tau_k * r_k * k
+ tau_d * d / p
- tau_ic * i
- tau_k * delta_tau * k_tau
- tau_k * e_tau * i
- g
) / (1 + r * (1 - tau_i));

// 8. Transfers
x = (                                     // Eq. (1.46): Transfer rule
  tau_c * (y - i)
+ tau_l * w * l / p
+ b / p
+ tau_k * r_k * k
+ tau_d * d / p
- tau_ic * i
- tau_k * delta_tau * k_tau
- tau_k * e_tau * i
- b * (1 + r * (1 - tau_i)) / p
) * gamma_x;

// 9. Consumption from resource constraint
c = y - i - g;                            // Eq. (1.5): Resource constraint (Y = C + I + G)

// 10. Marginal utility of consumption
lambda = zetta * c^(zetta * (1 - siggma) - 1) *
       (1 - l)^((1 - zetta) * (1 - siggma)) / (1 + tau_c); // (from consumption FOC — derived from utility function)

//lambda = (1 - zetta) * c^(zetta * (1 - siggma)) *
       //(1 - l)^((1 - zetta) * (1 - siggma) - 1) /   // Eq. (1.9): Labor FOC
       //(w * (1 - tau_l));

// 11. Capital Euler equation
q = lambda * (1 - tau_ic - tau_k * e_tau - betta * tau_k * delta_tau * (1 - e_tau));  // Eq. (1.11): FOC w.r.t. capital investment

// 12. Price of adjusting firms
p_star = ((1 + lambda_bar_p) / (1 - betta * theta)) * mc_bar;  // Eq. (1.22): Calvo with CES and markup

end;

/* 
initval;
// 1. Exogenous shocks
z         = 0;
epsilon_b = 0;
lambda_p  = lambda_bar_p;

// 2. Normalizations
p         = 1;
v         = 1;

// 3. Prices and rates
r         = rbar;
r_k       = r_k_bar;
w         = wbar;
l         = lbar;
mc        = mc_bar;

// 4. Output and capital
y         = ybar;
k         = (alfa / (1 - alfa)) * (w * l / r_k);
k_tau     = delta_0 * k * (1 - e_tau) / delta_tau;
i         = delta_0 * k;

// 5. Dividends
d         = y - w * l - r_k * k;

// 6. Government variables
g = (tau_c * (y - i)
 + tau_l * w * l / p
 + tau_k * r_k * k
 + tau_d * d / p
 - tau_ic * i
 - tau_k * delta_tau * k_tau
 - tau_k * e_tau * i) * (1 - gamma_x);

b = (tau_c * (y - i)
 + tau_l * w * l / p
 + tau_k * r_k * k
 + tau_d * d / p
 - tau_ic * i
 - tau_k * delta_tau * k_tau
 - tau_k * e_tau * i
 - g) / (1 + r * (1 - tau_i));

x = (tau_c * (y - i)
 + tau_l * w * l / p
 + b / p
 + tau_k * r_k * k
 + tau_d * d / p
 - tau_ic * i
 - tau_k * delta_tau * k_tau
 - tau_k * e_tau * i
 - b * (1 + r * (1 - tau_i)) / p) * gamma_x;

// 7. Consumption from resource constraint
c = y - i - g;

// 8. Marginal utility
lambda = zetta * c^(zetta * (1 - siggma) - 1) * (1 - l)^((1 - zetta) * (1 - siggma)) / (1 + tau_c);

// 9. Capital Euler equation
q = lambda * (1 - tau_ic - tau_k * e_tau - betta * tau_k * delta_tau * (1 - e_tau));

// 10. Calvo price
p_star = ((1 + lambda_bar_p) / (1 - betta * theta)) * mc_bar;
end;
*/

resid(1);
steady(nocheck);
check;



/***************************************************************/
/* Shocks */
/***************************************************************/
// The model has three shocks:
// 1. Shock to time preference (value of consumption across periods)
// 2. Shock to TFP
// 3. Shock to price markup
/***************************************************************/
shocks;
var e_b = siggma_b^2 ;
var e_z = siggma_z^2 ;
var e_p = siggma_p^2 ;
end;

/***************************************************************/
/* Computation */
/***************************************************************/
// Solve the model with second order taylor approximation.
// Get all of output.
// Specify seed to can confirm shock working.
stoch_simul(order=2, irf=20);











