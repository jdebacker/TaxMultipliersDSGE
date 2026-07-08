function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 112);

T = tax_dsge_est.dynamic_g1_tt(T, y, x, params, steady_state, it_);

T(101) = getPowerDeriv(T(48),1-params(6),2);
T(102) = T(81)*T(46)*getPowerDeriv(1-y(27),1-T(45),2)+T(80)*T(80)*T(101);
T(103) = T(47)*T(87)*T(47)*T(87)*T(101)+T(81)*T(47)*getPowerDeriv(y(29),T(45),2);
T(104) = getPowerDeriv(y(26)*y(3),params(3),2);
T(105) = getPowerDeriv(y(27),1-params(3),2);
T(106) = getPowerDeriv(y(26)*y(3),params(3)-1,2);
T(107) = (-((-(y(24)*y(59)))*(y(2)+y(2))))/(y(2)*y(2)*y(2)*y(2));
T(108) = (-(T(9)*(T(65)*T(107)+T(64)*2*T(64))));
T(109) = (-((-y(84))*(y(33)+y(33))))/(y(33)*y(33)*y(33)*y(33));
T(110) = (-((-((-(y(9)*y(48)*(y(8)-1)))*(y(38)+y(38))))/(y(38)*y(38)*y(38)*y(38))));
T(111) = (-(T(85)*T(85)))/(T(56)*T(56));
T(112) = (-(1/(params(20)*T(43))*1/(params(20)*T(43))))/(T(54)*T(54));

end
