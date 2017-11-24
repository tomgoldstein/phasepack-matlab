% This function automatically chooses the axis types (linear or log) for
% your plot, and then plots the curves.

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

