function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 27);

T(1) = (1-params(44)-params(15)*params(56)-params(1)*params(15)*params(60)/(1-params(1)*(1-params(60)))*(1-params(56)))*(1/params(1)-1+params(2));
T(2) = T(1)/(1-params(15));
T(3) = params(2)*params(3)/T(2)*(1-params(61));
T(4) = 1+1/(params(4)*((1-params(19)-T(3))*(1-params(11))/(params(1)*params(11)*(params(12)-1)/(params(1)*params(12)-1)-1)+params(19)*(1-params(13))/(params(1)*params(13)*(params(14)-1)/(params(1)*params(14)-1)-1)-T(3)));
T(5) = (T(2)/(params(3)*T(4)))^(1/(params(3)-1));
T(6) = T(4)*(1-params(3))*(T(5))^params(3);
T(7) = (params(18)*T(6)+T(2)*params(18)*T(5))/(1-params(61));
T(8) = (params(17)/params(1)-1)/(1-params(48));
T(9) = (T(7)-params(2)*params(18)*T(5)-params(19)*T(7))*params(40)+params(18)*T(6)*params(16)+params(15)*(T(2)*params(18)*T(5)-(1-params(56))*params(2)*params(18)*T(5)-params(56)*params(2)*params(18)*T(5))+T(7)*params(20)*params(48)*T(8)/params(17)+params(61)*T(7)*params(52)-params(44)*params(2)*params(18)*T(5)+T(7)*params(20)*(1-(1+T(8))/params(17))-params(19)*T(7);
T(10) = (1-params(11))*(T(7)-params(2)*params(18)*T(5)-params(19)*T(7))*(1+params(40))/((1-params(11))*(T(7)-params(2)*params(18)*T(5)-params(19)*T(7))*(1+params(40))+(1-params(18))*T(6)/(params(5)/((params(5)-1)*(1-params(16)))));
T(11) = y(8)^T(10);
T(12) = (1-y(6))^(1-T(10));
T(13) = T(11)*T(12);
T(14) = T(13)^(1-params(6));
T(15) = (y(5)*y(4))^params(3);
T(16) = y(6)^(1-params(3));
T(17) = y(36)*T(15)*T(16)-(params(18)^(1-params(3))*(params(18)*T(5))^params(3)-T(7));
T(18) = params(9)/2*(y(17)-params(17))^2;
T(19) = y(6)^(-params(3));
T(20) = (y(5)*y(4))^(params(3)-1);
T(21) = params(3)*y(36)*y(16)*T(20);
T(22) = params(10)/2*(y(18)-params(17))^2;
T(23) = (y(5)^(1+params(8))-1)*T(1)/(1+params(8));
T(24) = 1-params(7)/2*(y(38)-1)^2;
T(25) = log(y(20)/(T(7)*params(20)));
T(26) = log(y(7)/T(7));
T(27) = log(y(21)/T(9));

end
