// /**** TO DO IN SIMPLE VERSION ***/
// 3. Add Taylor rule for monetary authority
// 4. Add law of motion for capital stock (will include e_tau and delta_tau)


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


// /***************************************************************/
// /* Preamble */
// /***************************************************************/
// /* Endogenous Variables (20????) */
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
/*  Parameters (28) */
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
var z v epsilon_b r b r_k l k i y w d c x g lambda q p_star lambda_p mc p k_tau ;

varexo e_b e_z e_p ;

parameters betta siggma zetta siggma_g chi_g rho_b gamma delta_0 delta_1 delta_2 alfa lambda_bar_p theta phi_1 phi_2 rho_r rho_z rho_lambda siggma_z siggma_b siggma_p tau_c tau_i tau_k tau_l tau_d tau_ic e_tau delta_tau gamma_x lbar ybar rbar r_k_bar wbar mc_bar ;

betta = 0.9883 ;
siggma = 2 ;
zetta = 0.29 ;
siggma_g = 1.1 ;
chi_g = 0.5 ;
rho_b = 0.94 ;
gamma = 2.5 ;
delta_0 = 0.025 ;
delta_1 = 1/betta-1+delta_0 ; // get this from SS condition
psi = 0.32 ;
delta_2 = (psi/(1-psi))*delta_1 ;
alfa = 0.33 ;
lambda_bar_p = 8 ;
theta = 0.82 ;
rho_z = 0.89 ;
rho_lambda = 0 ;
siggma_z = 0 ; //0.64 ;
siggma_p = 0 ; //0.19 ;
siggma_b = 0 ; //2.4 ;

phi_1 = 0.280 ;
phi_2 = -0.006 ;
rho_r = 0.857 ;

tau_c = 0.022 ;
tau_i = 0.108 ;
tau_k = 0.214 ;
tau_l = 0.178 ;
tau_d = 0.108 ;
tau_ic = 0 ;
e_tau = 0 ;
delta_tau = 0.031 ;

gamma_x = 0.563 ;

lbar = 0.3333 ;

// Create steady state variables needed for model equations
rbar = (1-betta)/betta ;
r_k_bar = (1/betta)-1+delta_0 ;
wbar = (1-alfa)*((r_k_bar/alfa)^((-1*alfa)/(1-alfa))) ;
kbar = (alfa/(1-alfa))*wbar*(lbar/r_k_bar) ;
ybar = (kbar^alfa)*(lbar^(1-alfa)) ;
mc_bar = ((wbar^(1-alfa))*(r_k_bar^alfa))/((alfa^alfa)*((1-alfa)^(1-alfa))) ;



/***************************************************************/
/* Model */
/***************************************************************/
// Equations that define the model:
// 1. HH Budget Constraint -> consumption
// 2. HH FOC, consumption -> lambda (Marg Util Cons)
// 3. HH FOC, labor supply -> labor supply
// 4. HH FOC, capital -> investment
// 5. HH FOC, investment -> q
// 6. HH FOC, bond holdings -> demand for gov't bonds
// 7. HH FOC, capital utilization -> capital utilization
// 8. Law of motion for capital stock -> capital stock
// 9. Law of motion for tax basis of capital stock -> tax basis of cap stock
// 10. Resource constraint -> y
// 11. Calibration of g/y -> g
// 12. Gov't budget constraint -> x
// 13. Taylor Rule for monetary authority -> r
// 14. Int. goods producer FOC, effective capital demand -> (with market clearning condition) r_k
// 15. Int. goods producer FOC, labor demand -> (with market clearing condition) -> w
// 16. Int. goods producr FOCs -> mc (intermediate variable)
// 17. Int. goods producer FOC, price -> p_star
// 18. Int. goods producer profit function -> d
// 19. Calvo pricing rule -> p
// 20. AR(1) process for TFP -> z
// 21. Price markup process -> lambda_p
// 22. AR(1) process for stochastic time preference -> epsilon_b
/***************************************************************/
model ;

// 1. HH Budget Constraint
c + i + (b/p) = (r_k*v*k(-1)) + ((1-tau_l)*w*l) + d + (((1+r)*b(-1))/p) + x  ;


// 2. HH FOC, consumption
lambda = exp(epsilon_b)*(zetta*(c^(zetta*(1-siggma)-1))*((1-l)^((1-zetta)*(1-siggma)))) ;

// 3. HH FOC, labor supply
exp(epsilon_b)*((1-zetta)*(c^(zetta*(1-siggma)))*((1-l)^(((1-zetta)*(1-siggma))-1))) = (lambda*w*(1-tau_l)) ;

// 4. HH FOC, capital
q = (betta*lambda(+1)*r_k(+1)*v(+1)) + (betta*q(+1)*(1-(delta_0+(delta_1*(v(+1)-1))+((delta_2/2)*((v(+1)-1)^2))))) ;


// 5. HH FOC, investment
lambda = (q*(1-((gamma/2)*(((i/i(-1))-1)^2))-(gamma*((i/i(-1))-1)*(i/i(-1))))) - (betta*((q(+1)*gamma*((i(+1)/i)-1)*((i(+1)/i)^2)))) ;

// 6. HH FOC, bond holdings
lambda(-1)/p(-1) = (betta*lambda*(1+r))/p ;

// 7. HH FOC, capital utilization
(lambda*r_k) = q*(delta_1+(delta_2*(v-1))) ;

// 8. Law of motion for capital stock
k = ((1-(delta_0+(delta_1*(v-1))+((delta_2/2)*((v-1)^2))))*k(-1)) + (i*(1-((gamma/2)*((i/i(-1))-1)^2))) ;

// 9. Law of motion for tax basis of capital stock
k_tau = ((1-delta_tau)*k_tau(-1)) + (i*(1-e_tau)) ;

// 10. Resource constraint
// y = c + k - (1-delta_0)*k(-1) ;

// 11. Gov't budget constraint
b = 0 ;
g = ((1-gamma_x)*x)/gamma_x ;


// 12. Calibrate x as a fraction of gov't spending
x = gamma_x*(tau_l*w*l) ;


// 13. Taylor Rule for monetary authority
//r(+1) = max(((1/betta)*((1+((p-p(-1))/p(-1)))^phi_1)*((y/ybar)^phi_2)*(((1+r)/(1+rbar))^rho_r)-1)/(1-tau_i),0) ;
//r(+1) = ((1/betta)*((1+((p-p(-1))/p(-1)))^phi_1)*((y/ybar)^phi_2)*(((1+r)/(1+rbar))^rho_r)-1) ;


// 14. Int. goods producer FOC, effective capital demand
k(-1) = (r_k/(alfa*exp(z)*(v^alfa)*(l^(1-alfa))))^(1/(alfa-1)) ;


y = exp(z)*((v*k(-1))^alfa)*(l^(1-alfa)) ;

// 15. Int. goods producer FOC, labor demand
w = (1-alfa)*(y/l) ;

// 16. Int goods producer marginal cost
mc = ((w^(1-alfa))*(r_k^alfa))/(p*exp(z)*(alfa^alfa)*((1-alfa)^(1-alfa))) ;


// 17. Int. goods producer FOC, price
//p_star = 0 ;
p_star = ((((1+lambda_bar_p)/lambda_bar_p)*mc_bar*(1-(betta*theta)))*(((mc-mc_bar)/mc_bar)+((p-1)/1)+((lambda_bar_p/(1+lambda_bar_p))*((lambda_p-lambda_bar_p)/lambda_bar_p)))) + (betta*theta*p_star(+1));
//p_star = ((((1+lambda_bar_p)*mc_bar*(1-(betta*theta)))/(1-((1-lambda_bar_p)*mc_bar*(1-(betta*theta))*(1-theta))))*(((mc-mc_bar)/mc_bar)+((theta*p(-1))-1)+((lambda_bar_p/(1+lambda_bar_p))*((lambda_p-lambda_bar_p)/lambda_bar_p)))) + (betta*theta*p_star(+1)) ;


// 18. Int. goods producer profit function
d = y - (r_k*k(-1)*v) - (w*l) ;

// 19. Calvo pricing rule
p = ((theta*(p(-1)^(-1/lambda_p))) + ((1-theta)*(exp(p_star)^(-1/lambda_p))))^(-lambda_p) ;

// 20. AR(1) process for TFP
z = (rho_z*z(-1)) + e_z ;

// 21. Price markup process
ln(lambda_p) = (rho_lambda*(ln(lambda_p(-1)))) + ((1-rho_lambda)*ln(lambda_bar_p)) + e_p ;

// 22. AR(1) process for stochastic time preference
epsilon_b = (rho_b*epsilon_b(-1)) + e_b ;


end ;


/**************************************************************/
/* Steady State */
/***************************************************************/
// Set steady state values.  Calls matlab program to compute.
/***************************************************************/
initval;
//tax_dsge2_steadystate ;

z = 0 ;
v = 1 ;
epsilon_b = 0 ;
r = (1/betta)-1 ;
r_k = (1/betta)-1+delta_0 ;
l = lbar ;
k = 8.7912 ;
i = 8.7912*delta_0 ;
y = 2.9444 ;
w = 0.2 ;
d = 0 ;
c = 2.2850 ;
x = 0.1 ;
g = 0.8 ;
lambda = 0.1915 ;
q = 0.1915 ;
p_star = 0 ;
lambda_p = lambda_bar_p ;
mc = 1 ;
p = ((1+lambda_bar_p)/lambda_bar_p)*mc ;
b = 0 ;

e_b = 0 ;
e_p = 0 ;
e_z = 0 ;

end ;

steady ;

check ;


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
/***************************************************************/
//stoch_simul ;
stoch_simul(order=2) ;












