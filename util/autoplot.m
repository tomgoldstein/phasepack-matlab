%                             autoplot.m
%
% This function generates a line plot from x-axis and y-axis data.  It
% automatically selects either a linear or log axis, and it chooses colors
% and markers for the lines to make them easily distinguisable.
%
% inputs:
%       xvals:  a column vector of horizontal axis data
%       yvals:  a matrix, each column of which contains the data to
%                   generate a line of the plot
%       curveNames:  a cell array of string containing the name of each
%               curve to appear in the legend.


function autoplot( xvals, yvals, curveNames)

%% figure out what kind of axis to use
% Measure correlation between x/y values and a counter.  Also measure
% correlation with the log-values.
count = (1:numel(xvals))';
xLinCorr = corrcoef([count xvals(:)]);       xLinCorr = min(abs(xLinCorr(:)));
xLogCorr = corrcoef([count log(xvals(:))]);  xLogCorr = min(abs(xLogCorr(:)));
yLinCorr = corrcoef([count yvals]);       yLinCorr = min(abs(yLinCorr(:)));
yLogCorr = corrcoef([count log(yvals)]);  yLogCorr = min(abs(yLogCorr(:)));
% Pick the type of axis that produces the largest correlation
useLogY = min(yvals(:))>0 && yLogCorr>yLinCorr;
useLogX = min(xvals(:))>0 && xLogCorr>xLinCorr;  

%% Choose plotting function with appropriate axes
if useLogY && useLogX
    myplot = @(x,y,optName,optVal) loglog(x,y,optName,optVal);
elseif useLogY && ~useLogX
    myplot = @(x,y,optName,optVal) semilogy(x,y,optName,optVal);
elseif ~useLogY && useLogX
    myplot = @(x,y,optName,optVal) semilogx(x,y,optName,optVal);
else
    myplot = @(x,y,optName,optVal) plot(x,y,optName,optVal);
end

%% Plot stuff
markers = {'.','o','x','*','s','d','v','p','h','+'}; 
lineTypes = {'-','-.','--'};

nm = numel(markers);
nl = numel(lineTypes);

figure;grid on;
for k=1:length(curveNames)
    h = myplot(xvals, yvals(:,k), 'DisplayName', curveNames{k});
    set(h,'linewidth',1.75, 'LineStyle', lineTypes{mod(k-1,nl)+1}, 'Marker', markers{mod(k-1,nm)+1}); 
    hold on;grid on;
end


end

