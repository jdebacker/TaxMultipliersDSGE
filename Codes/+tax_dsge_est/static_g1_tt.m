function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
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

assert(length(T) >= 44);

T = tax_dsge_est.static_resid_tt(T, y, x, params);

T(28) = 1/(params(19)*T(7))/(y(2)/(params(19)*T(7)));
T(29) = getPowerDeriv(y(5)*y(4),params(3),1);
T(30) = getPowerDeriv(y(5)*y(4),params(3)-1,1);
T(31) = T(1)/(1+params(8))*getPowerDeriv(y(5),1+params(8),1);
T(32) = getPowerDeriv(T(13),1-params(6),1);
T(33) = T(11)*(-(getPowerDeriv(1-y(6),1-T(10),1)))*T(32);
T(34) = getPowerDeriv(y(6),1-params(3),1);
T(35) = 1/T(7)/(y(7)/T(7));
T(36) = T(32)*T(12)*getPowerDeriv(y(8),T(10),1);
T(37) = 1/(1+T(8))/(y(19)/(1+T(8)));
T(38) = 1/(T(7)*params(20))/(y(20)/(T(7)*params(20)));
T(39) = 1/T(9)/(y(21)/T(9));
T(40) = 1/params(15)/(y(23)/params(15));
T(41) = 1/params(16)/(y(24)/params(16));
T(42) = 1/params(40)/(y(25)/params(40));
T(43) = 1/params(48)/(y(27)/params(48));
T(44) = 1/params(52)/(y(28)/params(52));

end
