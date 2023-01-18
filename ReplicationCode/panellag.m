%
%function [output] =  panellag(data,I,length);
%
% ==================================================================
% Lags the whole Panel data by the wished number of lags.  
% Required inputs are
%  - data:   The panel data with dimenstion (time*I - by- variables)
%  - I:      number of cross sectional units
%  - length: the number of lags that should be taken from the data
%
% Output
%  - output: a matrix of dimension (time*I - by- variables) with the 
%            lagged data
% ===================================================================
% See also paneldiff, paneldetrend
%
% Written by Pascal Towbin and Sebastian Weber - 20.02.2009
%

function [output] =  panellag(data,I,length);

%
[T K]=size(data);
time = size(data,1)/I;
if time ~= floor(time);
    display('Error: The time dimension is not identical across individuals!!')
end
for actual_country=1:I
    if length > 0
      start=1+(actual_country-1)*time;
      ende=actual_country*time;
      output(start:start-1+length,1:K)=nan;
      output(length+start:ende,:)=data(start:ende-length,:);
    else
      start=1+(actual_country-1)*time;
      ende=actual_country*time;
      output(start:start-1+length,1:K)=nan;
      output(start-length:ende,:)=data(start-length:ende,:);
    end
end;

%=====================================================================