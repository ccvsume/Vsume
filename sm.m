function [fitresult,gof] = sm(ind, temzfor)
[xData, yData] = prepareCurveData( ind, temzfor );

ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.7;

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft, opts );