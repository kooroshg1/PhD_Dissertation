function out = calcRMSE(y, yHat)

y(isnan(y)) = [];
yHat(isnan(yHat)) = [];

out = (sum((yHat - y).^2) / length(y)) / (max(y) - min(y));