function metrics = runAdaptiveVECMModel(prices, baseWindow, volMult, shrinkFact, growFact, stepSize, minWindow, maxWindow)
% Adaptive-window VECM regime modeling based on volatility
% Inputs:
%   prices     : table or numeric matrix of price data
%   baseWindow : initial rolling window size
%   volMult    : factor to compare window volatility to global average
%   shrinkFact : scale factor to shrink the window if volatility is high
%   growFact   : scale factor to grow the window if volatility is below avg
%   stepSize   : move forward in time by 'stepSize' rows each iteration
%   minWindow  : lower bound on window size
%   maxWindow  : upper bound on window size
%
% Output:
%   metrics : struct with fields:
%       .regimeStabilityIndex
%       .rankInstability
%       .normalizedMSE

    nObs = height(prices);
    if nObs < (baseWindow + stepSize + 1)
        warning('Dataset too small: baseWindow + stepSize exceeds data length.');
        metrics.regimeStabilityIndex = NaN;
        metrics.rankInstability      = NaN;
        metrics.normalizedMSE        = NaN;
        return;
    end

    % Overall dataset volatility for reference
    allReturns = tick2ret(prices{:,:}, 'Method','continuous');
    avgVol     = std(allReturns(:));

    total_error  = 0;
    rank_changes = 0;
    prev_rank    = -1;
    errors       = [];
    count        = 0;
    start_idx    = 1;

    while (start_idx + baseWindow) <= (nObs - stepSize)
        w       = round(baseWindow);
        end_idx = start_idx + w - 1;
        if end_idx >= nObs
            break;
        end

        train_data = prices{start_idx:end_idx, :};
        test_data  = prices{end_idx+1 : end_idx+stepSize, :};

        % --- Check volatility in the training window ---
        train_rets = tick2ret(train_data, 'Method','continuous');
        windowVol   = std(train_rets(:));

        % --- Adapt baseWindow based on volatility ---
        if (windowVol > volMult * avgVol) && (baseWindow > minWindow)
            baseWindow = max(minWindow, baseWindow * shrinkFact);
        elseif (windowVol < avgVol) && (baseWindow < maxWindow)
            baseWindow = min(maxWindow, baseWindow * growFact);
        end

        % --- Johansen test + Bartlett correction ---
        [h, ~, stat, cValue, mles] = jcitest(train_data,'Model','H1','Lags',4);

        % Identify numeric 'r#' columns
        testedRanks     = str2double(erase(h.Properties.VariableNames(1:end-4), 'r'));
        bartlettFactors = computeBartlett(mles, testedRanks, train_data);

        % Adjust the 'r#' columns in 'stat' with bartlettFactors
        adjustedStat = stat;
        varNamesStat = adjustedStat.Properties.VariableNames;
        isRankVarStat = startsWith(varNamesStat, "r");  % e.g. r0, r1, r2...

        % Apply Bartlett to the numeric columns
        rankCols = varNamesStat(isRankVarStat); 
        for iRank = 1:length(rankCols)
            colName = rankCols{iRank};
            adjustedStat{1, colName} = stat{1, colName} / bartlettFactors(iRank);
        end
        
        % Keep only 'r#' columns in adjustedStat
        adjustedStat = adjustedStat(:, isRankVarStat);

        % Also keep only 'r#' columns in cValue
        varNamesCV  = cValue.Properties.VariableNames;
        isRankVarCV = startsWith(varNamesCV, "r");
        cValue      = cValue(:, isRankVarCV);

        % Convert to numeric arrays for rank comparison
        rankCandidate = adjustedStat{1, :};        
        critVals      = table2array(cValue(1, :));
        rank          = sum(rankCandidate > critVals);

        if rank ~= 0
            count = count + 1;
            if (prev_rank ~= -1) && (rank ~= prev_rank)
                rank_changes = rank_changes + 1;
            end
            prev_rank = rank;

            % --- VECM estimation (no HAC) ---
            vecmMdl = vecm(size(train_data,2), rank, 1);
            vecm_est= estimate(vecmMdl, train_data);

            % Optionally, if you had E = infer(...) for something else, comment out:
            % E = infer(vecm_est, train_data);
            
            % Forecast + MSE
            [fcast, ~] = forecast(vecm_est, stepSize, test_data);
            MSE         = nanmean((fcast - test_data).^2, 'all');
            total_error = total_error + MSE;
            errors(end+1) = MSE;  %#ok<AGROW>
        end

        % Move to next iteration
        start_idx = start_idx + stepSize;
    end

    % --- Final metric calculation ---
    if count > 0
        rank_instability       = rank_changes / count;
        error_volatility       = std(errors) / mean(errors);
        normalized_MSE         = total_error / count;
        regime_stability_index = 0.7 * rank_instability + 0.3 * error_volatility;
    else
        % No windows had rank>0 => no valid cointegration or no iteration
        rank_instability       = 0; 
        normalized_MSE         = 0;
        regime_stability_index = 0;
        warning('No valid (rank>0) windows found in runAdaptiveVECMModel. All metrics set to 0.');
    end

    metrics.regimeStabilityIndex = regime_stability_index;
    metrics.rankInstability      = rank_instability;
    metrics.normalizedMSE        = normalized_MSE;
end
