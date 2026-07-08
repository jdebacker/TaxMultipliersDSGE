function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 100);

T = tax_dsge_est.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(63) = 1/(params(19)*T(43));
T(64) = (-(y(24)*y(59)))/(y(2)*y(2));
T(65) = 2*T(10);
T(66) = y(59)/y(2);
T(67) = (-y(82))/(y(24)*y(24));
T(68) = 2*T(14);
T(69) = params(7)*y(95)*T(67)*T(68);
T(70) = (-(y(95)*y(82)))/(y(24)*y(24));
T(71) = 1/y(24);
T(72) = params(7)*y(95)*T(68)*T(71);
T(73) = y(95)/y(24);
T(74) = T(17)*T(72)+T(16)*T(73);
T(75) = getPowerDeriv(y(26)*y(3),params(3),1);
T(76) = getPowerDeriv(y(26)*y(3),params(3)-1,1);
T(77) = T(51)*getPowerDeriv(y(26),1+params(8),1);
T(78) = (1-y(90))*y(86)-T(51)*getPowerDeriv(y(83),1+params(8),1);
T(79) = (-(getPowerDeriv(1-y(27),1-T(45),1)));
T(80) = T(46)*T(79);
T(81) = getPowerDeriv(T(48),1-params(6),1);
T(82) = T(80)*T(81);
T(83) = getPowerDeriv(y(27),1-params(3),1);
T(84) = getPowerDeriv(y(27),(-params(3)),1);
T(85) = 1/T(43);
T(86) = T(85)/T(56);
T(87) = getPowerDeriv(y(29),T(45),1);
T(88) = T(81)*T(47)*T(87);
T(89) = y(29)*(1+y(46))*y(29)*(1+y(46));
T(90) = params(1)*(-y(84))/(y(33)*y(33));
T(91) = params(1)*1/y(33);
T(92) = 1/(1+T(1));
T(93) = 1/(params(20)*T(43))/T(54);
T(94) = 1/T(44);
T(95) = 1/params(15);
T(96) = 1/params(16);
T(97) = 1/params(40);
T(98) = 1/params(48);
T(99) = 1/params(52);
T(100) = y(24)/y(2);

end
