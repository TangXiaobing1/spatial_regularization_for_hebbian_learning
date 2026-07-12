function [pdTable, bestPd] = compute_ISI_models(ISI, plotFlag)
% COMPUTE_ISI_MODELS_V3 Robust multi-distribution fitting for ISI data
% (v2) Adds RMSE calculation and modifies plot style to hollow circles.
% -------------------------------------------------------------------------
% COMPUTE_ISI_MODELS Robust multi‑distribution fitting for ISI data
% [pdTable, bestPd] = COMPUTE_ISI_MODELS(ISI, plotFlag) fits four
% candidate distributions (Weibull, Lognormal, Gamma, Exponential) to
% positive inter‑spine interval values and returns a summary table plus
% the best distribution (smallest AIC). If plotFlag (default true) is
% enabled, the function draws a normalised histogram overlaid with the
% pdf of bestPd.
% -------------------------------------------------------------------------
% • Filters out non‑positive ISI values (0 or negative) that break
% maximum‑likelihood fits. A warning is issued if values are removed.
% • Gracefully skips models whose fitting fails and sets their AIC/BIC
% to Inf.
% • Calculates AIC = 2k − 2logL and BIC = k·log(n) − 2logL.
% • Returns a ranked pdTable for easy inspection.
% -------------------------------------------------------------------------
% Example
% ISI = diff(sort(Pvec_tree(tree)(spine_idx)));
% [pdTable, bestPd] = compute_ISI_models(ISI);
% -------------------------------------------------------------------------
% painting option: when plotFlag = 1(default), draw all 4 models on plot.
% when plotFlag = 0, draw only best fit model.
% -------------------------------------------------------------------------

if nargin < 2, plotFlag = true; end

ISI = double(ISI(:));
negOrZero = ISI <= 0.016 | ISI > 30; 
if any(negOrZero)
    warning('移除了 %d 个非正数的ISI值。', sum(negOrZero));
    ISI = ISI(~negOrZero);
end

n = numel(ISI);
if n < 2, error('有效数据点少于2个，无法进行拟合。'); end

binWidth = 0.2;
[y_empirical, x_edges] = histcounts(ISI, 'Normalization','pdf','BinWidth', binWidth);
x_empirical = x_edges(1:end-1) + diff(x_edges)/2;


models = {'Weibull','Lognormal','Gamma','Exponential'};
pdObjs = cell(size(models));
logL   = NaN(size(models));
params = cell(size(models));
pValue = NaN(size(models));
RMSE   = NaN(size(models)); 

for i = 1:numel(models)
    try
        pdObjs{i} = fitdist(ISI, models{i});
        logL(i)   = -pdObjs{i}.NLogL;
        
        switch models{i}
            case 'Weibull',   params{i} = [pdObjs{i}.B, pdObjs{i}.A];
            case 'Lognormal', params{i} = [pdObjs{i}.mu, pdObjs{i}.sigma];
            case 'Gamma',     params{i} = [pdObjs{i}.a, pdObjs{i}.b];
            case 'Exponential', params{i} = pdObjs{i}.mu;
        end
        
        [~, pValue(i)] = kstest(ISI, 'CDF', pdObjs{i});
        
        % --- ComputeRMSE ---
        y_predicted = pdf(pdObjs{i}, x_empirical);
        squared_error = (y_empirical - y_predicted).^2;
        RMSE(i) = sqrt(mean(squared_error));

    catch ME
        warning('无法拟合 %s 模型: %s', models{i}, ME.message);
    end
end

k = cellfun(@(p) numel(p), params);
AIC = 2.*k - 2.*logL;
BIC = k.*log(n) - 2.*logL;

AIC(isnan(AIC)) = Inf; 
BIC(isnan(BIC)) = Inf;
pValue(isnan(pValue)) = 0;
RMSE(isnan(RMSE)) = Inf;

pdTable = table(models', params', logL', AIC', BIC', pValue', RMSE', ...
    'VariableNames', {'Model','Params','LogL','AIC','BIC','PValue','RMSE'});

[~, bestIdx] = min(AIC);
bestPd = pdObjs{bestIdx};

if plotFlag && ~isempty(bestPd)
    figure; hold on;
    
    plot(x_empirical, y_empirical, 'o', ...
        'MarkerFaceColor', [0.3 0.3 0.3], ...
        'MarkerEdgeColor', [0 0 1], ...
        'MarkerSize', 6, ...
        'DisplayName', 'ISI Data (Empirical PDF)');
    
    x_curve = linspace(min(ISI), max(ISI), 200);
    colors = lines(numel(models)); 
    
    for i = 1:numel(models)
        if ~isempty(pdObjs{i})
            y_curve = pdf(pdObjs{i}, x_curve);
            
            lineWidth = 1.5;
            if i == bestIdx
                lineWidth = 2.5;
            end
            
            displayName = sprintf('%s (AIC: %.2f)', models{i}, AIC(i));
            
            plot(x_curve, y_curve, 'LineWidth', lineWidth, 'Color', colors(i,:), ...
                 'DisplayName', displayName);
        end
    end
    
    legend('show', 'Location', 'northeast');
    xlabel('Inter-Spine Interval (\mum)'); 
    ylabel('Probability Density');
    title('ISI Distribution and All Fitted Models');
    grid on;
    hold off;
end


if plotFlag == 0 && ~isempty(bestPd)
    figure; hold on;
    histogram(ISI, 'Normalization','pdf','BinWidth',0.2);
    x = linspace(min(ISI), max(ISI), 200);
    y = pdf(bestPd, x);
    plot(x, y, 'LineWidth',2);
    legend('ISI data', sprintf('Best fit: %s', models{bestIdx}));
    xlabel('Inter‑Spine Interval (\mum)'); ylabel('Probability Density');
    title('ISI distribution and best‑fit model');
end


[~, ord] = sort(pdTable.AIC);
pdTable = pdTable(ord, :);

end