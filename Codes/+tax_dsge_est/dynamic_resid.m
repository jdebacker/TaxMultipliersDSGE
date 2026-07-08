function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = tax_dsge_est.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(60, 1);
lhs = y(29);
rhs = y(22)-params(11)*y(5);
residual(1) = lhs - rhs;
lhs = y(30);
rhs = y(23)-params(13)*y(6);
residual(2) = lhs - rhs;
lhs = y(31);
rhs = params(12)*y(5)+y(22)*(1-params(12));
residual(3) = lhs - rhs;
lhs = y(32);
rhs = params(14)*y(6)+y(23)*(1-params(14));
residual(4) = lhs - rhs;
lhs = y(33);
rhs = y(58)*T(45)*T(49)/(y(29)*(1+y(46)));
residual(5) = lhs - rhs;
lhs = T(49)*y(58)*(1-T(45))/(1-y(27));
rhs = y(33)*y(36)/y(60);
residual(6) = lhs - rhs;
lhs = y(28);
rhs = T(50)-T(4);
residual(7) = lhs - rhs;
lhs = y(36);
rhs = T(2)*(1-params(3))*y(57)*y(37)*T(5);
residual(8) = lhs - rhs;
lhs = y(35);
rhs = T(3)*T(7);
residual(9) = lhs - rhs;
lhs = T(50);
rhs = y(36)*T(8)+T(4)+y(22)+y(23)+y(24)+y(3)*T(52);
residual(10) = lhs - rhs;
lhs = y(25);
rhs = y(3)*(1-params(2))+y(24)*T(11);
residual(11) = lhs - rhs;
lhs = y(50);
rhs = (1-params(60))*y(17)+y(24)*(1-y(53));
residual(12) = lhs - rhs;
lhs = y(33)*y(51);
rhs = params(1)*y(84)*((1-params(60))*y(92)+params(60)*y(90));
residual(13) = lhs - rhs;
lhs = y(35)*(1-y(44));
rhs = y(26)^params(8)*T(37);
residual(14) = lhs - rhs;
lhs = y(33)*y(34);
rhs = params(1)*y(84)*T(53);
residual(15) = lhs - rhs;
lhs = y(33)*(1-y(47)-y(53)*y(44));
rhs = y(33)*y(34)*T(13)+params(1)*y(84)*y(85)*T(18)+(1-y(53))*y(33)*y(51);
residual(16) = lhs - rhs;
lhs = y(33);
rhs = y(84)*params(1)*(1+(y(40)-1)*(1-y(91)))/y(88);
residual(17) = lhs - rhs;
lhs = y(39);
rhs = y(38)*y(36)/y(7);
residual(18) = lhs - rhs;
lhs = y(27)*(params(5)-1)*(1-y(45))+(y(39)-params(17))*params(10)*y(39)-params(5)*y(27)/y(60);
rhs = T(19)*y(89)*(y(89)-params(17));
residual(19) = lhs - rhs;
residual(20) = 1-y(37)-y(56);
lhs = (1-y(37)-y(54))/(params(12)-1);
rhs = T(20)*T(21);
residual(21) = lhs - rhs;
lhs = (1-y(37)-y(55))/(params(14)-1);
rhs = T(20)*T(22);
residual(22) = lhs - rhs;
lhs = params(4)*(y(29)*y(54)+y(30)*y(55)+y(56)*(y(28)-y(22)-y(23)))+(y(38)-params(17))*params(9)*y(38)-y(28);
rhs = y(88)*T(23)*(y(88)-params(17));
residual(23) = lhs - rhs;
lhs = y(52);
rhs = y(28)-y(27)*y(36)-y(3)*y(26)*y(35);
residual(24) = lhs - rhs;
lhs = y(41);
rhs = y(23)+y(8)*y(9)/y(38)+y(42)-y(43)+x(it_, 14);
residual(25) = lhs - rhs;
lhs = y(43);
rhs = y(22)*y(46)+y(27)*y(36)*y(45)+y(44)*(y(3)*y(26)*y(35)-params(60)*y(17)-y(24)*y(53))+y(9)*y(48)*(y(8)-1)/y(38)+y(52)*y(49)-y(24)*y(47);
residual(26) = lhs - rhs;
lhs = log(T(24));
rhs = x(it_, 3)+params(33)*log(T(25))+params(34)*T(55)+params(35)*T(57);
residual(27) = lhs - rhs;
lhs = log(T(26));
rhs = x(it_, 4)+x(it_, 3)*params(39)+params(36)*log(T(27))+params(37)*T(55)+params(38)*T(57);
residual(28) = lhs - rhs;
lhs = log(T(58));
rhs = x(it_, 1)+params(27)*log(T(59))-params(28)*T(55)+params(29)*T(57);
residual(29) = lhs - rhs;
lhs = log(T(60));
rhs = x(it_, 2)+params(30)*log(T(61))-params(31)*T(55)+params(32)*T(57);
residual(30) = lhs - rhs;
lhs = log(T(35));
rhs = x(it_, 8)+params(21)*log(T(36))+(1-params(21))*(params(22)*log(T(28))+params(23)*log(T(62)));
residual(31) = lhs - rhs;
lhs = log(T(29));
rhs = x(it_, 9)+params(41)*log(T(30))+params(42)*T(55)+params(43)*T(57);
residual(32) = lhs - rhs;
lhs = y(47)-params(44);
rhs = x(it_, 10)+params(45)*(y(14)-params(44))+params(46)*T(55)+params(47)*T(57);
residual(33) = lhs - rhs;
lhs = log(T(31));
rhs = x(it_, 11)+params(49)*log(T(32))+params(50)*T(55)+params(51)*T(57);
residual(34) = lhs - rhs;
lhs = log(T(33));
rhs = x(it_, 12)+params(53)*log(T(34))+params(54)*T(55)+params(55)*T(57);
residual(35) = lhs - rhs;
lhs = y(53)-params(56);
rhs = x(it_, 13)+params(57)*(y(18)-params(56))+params(58)*T(55)+params(59)*T(57);
residual(36) = lhs - rhs;
lhs = log(y(57));
rhs = params(24)*log(y(19))+x(it_, 5);
residual(37) = lhs - rhs;
lhs = log(y(58));
rhs = params(25)*log(y(20))+x(it_, 6);
residual(38) = lhs - rhs;
lhs = log(y(59));
rhs = params(26)*log(y(21))+x(it_, 7);
residual(39) = lhs - rhs;
lhs = y(61);
rhs = y(22)*y(46);
residual(40) = lhs - rhs;
lhs = y(62);
rhs = y(27)*y(36)*y(45);
residual(41) = lhs - rhs;
lhs = y(63);
rhs = y(44)*(y(3)*y(26)*y(35)-params(60)*y(17)-y(24)*y(53));
residual(42) = lhs - rhs;
lhs = y(64);
rhs = y(9)*y(48)*(y(8)-1)/y(38);
residual(43) = lhs - rhs;
lhs = y(65);
rhs = y(52)*y(49);
residual(44) = lhs - rhs;
lhs = y(66);
rhs = y(24)*y(47);
residual(45) = lhs - rhs;
lhs = y(67);
rhs = y(65)+y(64)+y(63)+y(61)+y(62)-y(66);
residual(46) = lhs - rhs;
lhs = y(68);
rhs = 100*(y(22)-(T(43)-params(2)*params(18)*T(41)-params(19)*T(43)))/(T(43)-params(2)*params(18)*T(41)-params(19)*T(43));
residual(47) = lhs - rhs;
lhs = y(69);
rhs = 100*(y(24)-params(2)*params(18)*T(41))/(params(2)*params(18)*T(41));
residual(48) = lhs - rhs;
lhs = y(70);
rhs = 100*(y(23)-params(19)*T(43))/(params(19)*T(43));
residual(49) = lhs - rhs;
lhs = y(71);
rhs = y(38)-params(17);
residual(50) = lhs - rhs;
lhs = y(72);
rhs = y(40)-(1+T(1));
residual(51) = lhs - rhs;
lhs = y(73);
rhs = (y(44)-params(15))/params(15);
residual(52) = lhs - rhs;
lhs = y(74);
rhs = (y(45)-params(16))/params(16);
residual(53) = lhs - rhs;
lhs = y(75);
rhs = 100*(y(41)-params(20)*T(43))/(params(20)*T(43));
residual(54) = lhs - rhs;
lhs = y(76);
rhs = 100*(y(42)-T(44))/T(44);
residual(55) = lhs - rhs;
lhs = y(77);
rhs = (y(46)-params(40))/params(40);
residual(56) = lhs - rhs;
lhs = y(78);
rhs = (y(49)-params(52))/params(52);
residual(57) = lhs - rhs;
lhs = y(81);
rhs = (y(48)-params(48))/params(48);
residual(58) = lhs - rhs;
lhs = y(79);
rhs = y(47)-params(44);
residual(59) = lhs - rhs;
lhs = y(80);
rhs = y(53)-params(56);
residual(60) = lhs - rhs;

end
