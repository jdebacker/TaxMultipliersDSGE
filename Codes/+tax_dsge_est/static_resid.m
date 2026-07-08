function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = tax_dsge_est.static_resid_tt(T, y, x, params);
end
residual = zeros(60, 1);
lhs = y(8);
rhs = y(1)-params(11)*y(10);
residual(1) = lhs - rhs;
lhs = y(9);
rhs = y(2)-params(13)*y(11);
residual(2) = lhs - rhs;
lhs = y(10);
rhs = params(12)*y(10)+y(1)*(1-params(12));
residual(3) = lhs - rhs;
lhs = y(11);
rhs = params(14)*y(11)+y(2)*(1-params(14));
residual(4) = lhs - rhs;
lhs = y(12);
rhs = T(10)*y(37)*T(14)/(y(8)*(1+y(25)));
residual(5) = lhs - rhs;
lhs = T(14)*(1-T(10))*y(37)/(1-y(6));
rhs = y(12)*y(15)/y(39);
residual(6) = lhs - rhs;
lhs = y(7);
rhs = T(17)-T(18);
residual(7) = lhs - rhs;
lhs = y(15);
rhs = T(15)*(1-params(3))*y(36)*y(16)*T(19);
residual(8) = lhs - rhs;
lhs = y(14);
rhs = T(16)*T(21);
residual(9) = lhs - rhs;
lhs = T(17);
rhs = y(15)*T(22)+T(18)+y(1)+y(2)+y(3)+y(4)*T(23);
residual(10) = lhs - rhs;
lhs = y(4);
rhs = y(4)*(1-params(2))+y(3)*T(24);
residual(11) = lhs - rhs;
lhs = y(29);
rhs = (1-params(60))*y(29)+y(3)*(1-y(32));
residual(12) = lhs - rhs;
lhs = y(12)*y(30);
rhs = params(1)*y(12)*((1-params(60))*y(30)+params(60)*y(23));
residual(13) = lhs - rhs;
lhs = y(14)*(1-y(23));
rhs = T(1)*y(5)^params(8);
residual(14) = lhs - rhs;
lhs = y(12)*y(13);
rhs = params(1)*y(12)*((1-params(2))*y(13)+y(5)*y(14)*(1-y(23))-T(23));
residual(15) = lhs - rhs;
lhs = y(12)*(1-y(26)-y(32)*y(23));
rhs = y(12)*y(13)*(T(24)-(y(38)-1)*params(7)*y(38))+(y(38)-1)*params(7)*y(38)*params(1)*y(12)*y(13)+(1-y(32))*y(12)*y(30);
residual(16) = lhs - rhs;
lhs = y(12);
rhs = y(12)*params(1)*(1+(y(19)-1)*(1-y(27)))/y(17);
residual(17) = lhs - rhs;
lhs = y(18);
rhs = y(17);
residual(18) = lhs - rhs;
lhs = y(6)*(params(5)-1)*(1-y(24))+(y(18)-params(17))*params(10)*y(18)-params(5)*y(6)/y(39);
rhs = (y(18)-params(17))*y(18)*params(1)*params(10);
residual(19) = lhs - rhs;
residual(20) = 1-y(16)-y(35);
lhs = (1-y(16)-y(33))/(params(12)-1);
rhs = params(1)*(params(11)*y(33)+(1-y(16)-y(33))*params(12)/(params(12)-1));
residual(21) = lhs - rhs;
lhs = (1-y(16)-y(34))/(params(14)-1);
rhs = params(1)*(params(13)*y(34)+(1-y(16)-y(34))*params(14)/(params(14)-1));
residual(22) = lhs - rhs;
lhs = params(4)*(y(8)*y(33)+y(9)*y(34)+y(35)*(y(7)-y(1)-y(2)))+(y(17)-params(17))*params(9)*y(17)-y(7);
rhs = (y(17)-params(17))*y(17)*params(1)*params(9);
residual(23) = lhs - rhs;
lhs = y(31);
rhs = y(7)-y(6)*y(15)-y(4)*y(5)*y(14);
residual(24) = lhs - rhs;
lhs = y(20);
rhs = y(2)+y(19)*y(20)/y(17)+y(21)-y(22)+x(14);
residual(25) = lhs - rhs;
lhs = y(22);
rhs = y(1)*y(25)+y(6)*y(15)*y(24)+y(23)*(y(4)*y(5)*y(14)-params(60)*y(29)-y(3)*y(32))+y(20)*(y(19)-1)*y(27)/y(17)+y(31)*y(28)-y(3)*y(26);
residual(26) = lhs - rhs;
lhs = log(y(23)/params(15));
rhs = x(3)+log(y(23)/params(15))*params(33)+params(34)*T(25)+params(35)*T(26);
residual(27) = lhs - rhs;
lhs = log(y(24)/params(16));
rhs = x(4)+x(3)*params(39)+log(y(24)/params(16))*params(36)+T(25)*params(37)+T(26)*params(38);
residual(28) = lhs - rhs;
lhs = log(y(2)/(params(19)*T(7)));
rhs = x(1)+log(y(2)/(params(19)*T(7)))*params(27)-T(25)*params(28)+T(26)*params(29);
residual(29) = lhs - rhs;
lhs = T(27);
rhs = x(2)+T(27)*params(30)-T(25)*params(31)+T(26)*params(32);
residual(30) = lhs - rhs;
lhs = log(y(19)/(1+T(8)));
rhs = x(8)+log(y(19)/(1+T(8)))*params(21)+(1-params(21))*(params(22)*log(y(17)/params(17))+T(26)*params(23));
residual(31) = lhs - rhs;
lhs = log(y(25)/params(40));
rhs = x(9)+log(y(25)/params(40))*params(41)+T(25)*params(42)+T(26)*params(43);
residual(32) = lhs - rhs;
lhs = y(26)-params(44);
rhs = x(10)+(y(26)-params(44))*params(45)+T(25)*params(46)+T(26)*params(47);
residual(33) = lhs - rhs;
lhs = log(y(27)/params(48));
rhs = x(11)+log(y(27)/params(48))*params(49)+T(25)*params(50)+T(26)*params(51);
residual(34) = lhs - rhs;
lhs = log(y(28)/params(52));
rhs = x(12)+log(y(28)/params(52))*params(53)+T(25)*params(54)+T(26)*params(55);
residual(35) = lhs - rhs;
lhs = y(32)-params(56);
rhs = x(13)+(y(32)-params(56))*params(57)+T(25)*params(58)+T(26)*params(59);
residual(36) = lhs - rhs;
lhs = log(y(36));
rhs = log(y(36))*params(24)+x(5);
residual(37) = lhs - rhs;
lhs = log(y(37));
rhs = log(y(37))*params(25)+x(6);
residual(38) = lhs - rhs;
lhs = log(y(38));
rhs = log(y(38))*params(26)+x(7);
residual(39) = lhs - rhs;
lhs = y(40);
rhs = y(1)*y(25);
residual(40) = lhs - rhs;
lhs = y(41);
rhs = y(6)*y(15)*y(24);
residual(41) = lhs - rhs;
lhs = y(42);
rhs = y(23)*(y(4)*y(5)*y(14)-params(60)*y(29)-y(3)*y(32));
residual(42) = lhs - rhs;
lhs = y(43);
rhs = y(20)*(y(19)-1)*y(27)/y(17);
residual(43) = lhs - rhs;
lhs = y(44);
rhs = y(31)*y(28);
residual(44) = lhs - rhs;
lhs = y(45);
rhs = y(3)*y(26);
residual(45) = lhs - rhs;
lhs = y(46);
rhs = y(44)+y(43)+y(42)+y(40)+y(41)-y(45);
residual(46) = lhs - rhs;
lhs = y(47);
rhs = 100*(y(1)-(T(7)-params(2)*params(18)*T(5)-params(19)*T(7)))/(T(7)-params(2)*params(18)*T(5)-params(19)*T(7));
residual(47) = lhs - rhs;
lhs = y(48);
rhs = 100*(y(3)-params(2)*params(18)*T(5))/(params(2)*params(18)*T(5));
residual(48) = lhs - rhs;
lhs = y(49);
rhs = 100*(y(2)-params(19)*T(7))/(params(19)*T(7));
residual(49) = lhs - rhs;
lhs = y(50);
rhs = y(17)-params(17);
residual(50) = lhs - rhs;
lhs = y(51);
rhs = y(19)-(1+T(8));
residual(51) = lhs - rhs;
lhs = y(52);
rhs = (y(23)-params(15))/params(15);
residual(52) = lhs - rhs;
lhs = y(53);
rhs = (y(24)-params(16))/params(16);
residual(53) = lhs - rhs;
lhs = y(54);
rhs = 100*(y(20)-T(7)*params(20))/(T(7)*params(20));
residual(54) = lhs - rhs;
lhs = y(55);
rhs = 100*(y(21)-T(9))/T(9);
residual(55) = lhs - rhs;
lhs = y(56);
rhs = (y(25)-params(40))/params(40);
residual(56) = lhs - rhs;
lhs = y(57);
rhs = (y(28)-params(52))/params(52);
residual(57) = lhs - rhs;
lhs = y(60);
rhs = (y(27)-params(48))/params(48);
residual(58) = lhs - rhs;
lhs = y(58);
rhs = y(26)-params(44);
residual(59) = lhs - rhs;
lhs = y(59);
rhs = y(32)-params(56);
residual(60) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
