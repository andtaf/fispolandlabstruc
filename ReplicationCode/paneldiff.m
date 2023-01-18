%
%function [output] =  paneldiff(data,I);
%
% ============================================================
% Takes the first difference of data which is in panel format.
% Required input data includes:
%  - data: a matrix of dimension (time*I -by- variables)
%  - I:    number of cross sectional units
% 
% the output is the first differenced data of dimension
% (time*I -by- variables)
%
% ------------------------------------------------------------
% Calls on the function panellag
% ============================================================
% See also panellag, paneldetrend
%
% Written by Pascal Towbin and Sebastian Weber, 20.02.2009 
%

function [output] =  paneldiff(data,I,lag);

output=data-panellag(data,I,lag);

%=============================================================