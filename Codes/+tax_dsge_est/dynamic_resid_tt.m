function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 62);

T(1) = (params(17)/params(1)-1)/(1-params(48));
T(2) = (y(26)*y(3))^params(3);
T(3) = y(27)^(1-params(3));
T(4) = params(9)/2*(y(38)-params(17))^2;
T(5) = y(27)^(-params(3));
T(6) = (y(26)*y(3))^(params(3)-1);
T(7) = params(3)*y(57)*y(37)*T(6);
T(8) = params(10)/2*(y(39)-params(17))^2;
T(9) = params(7)/2;
T(10) = y(24)*y(59)/y(2)-1;
T(11) = 1-T(9)*T(10)^2;
T(12) = params(7)*y(24)*y(59)/y(2);
T(13) = T(11)-T(10)*T(12);
T(14) = y(82)/y(24);
T(15) = T(14)^2;
T(16) = params(7)*y(95)*T(15);
T(17) = y(95)*y(82)/y(24)-1;
T(18) = T(16)*T(17);
T(19) = params(1)*params(10)*y(84)/y(33);
T(20) = params(1)*y(84)/y(33);
T(21) = params(11)*y(93)+params(12)/(params(12)-1)*(1-y(87)-y(93));
T(22) = params(13)*y(94)+params(14)/(params(14)-1)*(1-y(87)-y(94));
T(23) = y(84)/y(33)*params(1)*params(9);
T(24) = y(44)/params(15);
T(25) = y(11)/params(15);
T(26) = y(45)/params(16);
T(27) = y(12)/params(16);
T(28) = y(38)/params(17);
T(29) = y(46)/params(40);
T(30) = y(13)/params(40);
T(31) = y(48)/params(48);
T(32) = y(15)/params(48);
T(33) = y(49)/params(52);
T(34) = y(16)/params(52);
T(35) = y(40)/(1+T(1));
T(36) = y(8)/(1+T(1));
T(37) = (1/params(1)-1+params(2))*(1-params(44)-params(15)*params(56)-params(1)*params(15)*params(60)/(1-params(1)*(1-params(60)))*(1-params(56)));
T(38) = T(37)/(1-params(15));
T(39) = (1-params(61))*params(2)*params(3)/T(38);
T(40) = 1+1/(params(4)*((1-params(11))*(1-params(19)-T(39))/(params(1)*params(11)*(params(12)-1)/(params(1)*params(12)-1)-1)+params(19)*(1-params(13))/(params(1)*params(13)*(params(14)-1)/(params(1)*params(14)-1)-1)-T(39)));
T(41) = (T(38)/(params(3)*T(40)))^(1/(params(3)-1));
T(42) = (1-params(3))*T(40)*(T(41))^params(3);
T(43) = (params(18)*T(42)+T(38)*params(18)*T(41))/(1-params(61));
T(44) = params(40)*(T(43)-params(2)*params(18)*T(41)-params(19)*T(43))+params(18)*params(16)*T(42)+params(15)*(T(38)*params(18)*T(41)-(1-params(56))*params(2)*params(18)*T(41)-params(56)*params(2)*params(18)*T(41))+params(20)*T(43)*params(48)*T(1)/params(17)+params(52)*params(61)*T(43)-params(44)*params(2)*params(18)*T(41)+params(20)*T(43)*(1-(1+T(1))/params(17))-params(19)*T(43);
T(45) = (1+params(40))*(1-params(11))*(T(43)-params(2)*params(18)*T(41)-params(19)*T(43))/((1+params(40))*(1-params(11))*(T(43)-params(2)*params(18)*T(41)-params(19)*T(43))+(1-params(18))*T(42)/(params(5)/((params(5)-1)*(1-params(16)))));
T(46) = y(29)^T(45);
T(47) = (1-y(27))^(1-T(45));
T(48) = T(46)*T(47);
T(49) = T(48)^(1-params(6));
T(50) = y(57)*T(2)*T(3)-(params(18)^(1-params(3))*(params(18)*T(41))^params(3)-T(43));
T(51) = T(37)/(1+params(8));
T(52) = (y(26)^(1+params(8))-1)*T(51);
T(53) = (1-params(2))*y(85)+(1-y(90))*y(86)*y(83)-(y(83)^(1+params(8))-1)*T(51);
T(54) = y(9)/(params(20)*T(43));
T(55) = log(T(54));
T(56) = y(4)/T(43);
T(57) = log(T(56));
T(58) = y(23)/(params(19)*T(43));
T(59) = y(1)/(params(19)*T(43));
T(60) = y(42)/T(44);
T(61) = y(10)/T(44);
T(62) = y(28)/T(43);

end
