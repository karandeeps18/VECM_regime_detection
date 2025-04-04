function metrics = runFixedVECMModel(prices, windowSize, stepSize)
% Fixed rolling-window VECM regime modeling
% Inputs:
%   prices     : table or numeric matrix of price data
%   windowSize : integer, size of rolling window
%   stepSize   : integer, how often to move the window
%
% Output:
%   metrics : struct with fields:
%       .regimeStabilityIndex
%       .rankInstability
%       .normalizedMSE

    % --- Basic dimension check ---
    nObs = height(prices);
    if nObs < (windowSize + stepSize + 1)
        warning('Dataset too small: windowSize + stepSize exceeds data length.');
        metrics.regimeStabilityIndex = NaN;
        metrics.rankInstability      = NaN;
        metrics.normalizedMSE        = NaN;
        return;
    end

    total_error  = 0;
    rank_changes = 0;
    prev_rank    = -1;
    errors       = [];
    count        = 0;

    % Start the rolling loop
    for start_idx = 1 : stepSize : (nObs - windowSize - stepSize)
        end_idx   = start_idx + windowSize - 1;

        train_data = prices{start_idx:end_idx, :};
        test_data  = prices{end_idx+1 : end_idx+stepSize, :};

        % ===== Johansen test with Bartlett correction =====
        [h, ~, stat, cValue, mles] = jcitest(train_data,'Model','H1','Lags',4);

        % 1) Identify numeric 'r#' columns from test
        testedRanks = str2double(erase(h.Properties.VariableNames(1:end-4), 'r'));
        bartlettFactors = computeBartlett(mles, testedRanks, train_data);

        % 2) Adjust the 'r#' columns in 'stat' with bartlettFactors
        adjustedStat = stat;
        varNamesStat = adjustedStat.Properties.VariableNames;
        isRankVarStat = startsWith(varNamesStat, "r");  % e.g. r0, r1, r2...

        % Only apply Bartlett to those numeric columns
        rankCols = varNamesStat(isRankVarStat); % cell array of 'r0','r1',...
        for iRank = 1:length(rankCols)
            colName = rankCols{iRank};
            adjustedStat{1, colName} = stat{1, colName} / bartlettFactors(iRank);
        end

        % Keep only the numeric 'r#' columns in adjustedStat
        adjustedStat = adjustedStat(:, isRankVarStat);

        % Do the same "keep only r# columns" in cValue
        varNamesCV = cValue.Properties.VariableNames;
        isRankVarCV = startsWith(varNamesCV, "r");
        cValue = cValue(:, isRankVarCV);

        % Now get numeric arrays for the rank comparison
        rankCandidate = adjustedStat{1, :};         % e.g. [r0_val, r1_val, ...]
        critVals      = table2array(cValue(1, :));  % e.g. [crit0, crit1, ...]
        rank          = sum(rankCandidate > critVals);

        % ===== If we have cointegration => do VECM, forecast, etc. =====
        if rank ~= 0
            count = count + 1;
            if (prev_rank ~= -1) && (rank ~= prev_rank)
                rank_changes = rank_changes + 1;
            end
            prev_rank = rank;

            % VECM + HAC
            vecmMdl = vecm(size(train_data,2), rank, 1);
            vecm_est= estimate(vecmMdl, train_data);
           % E       = infer(vecm_est, train_data);
           % max_lags= floor(4 * (size(train_data,1) / 100)^(2/9));
           % hac_cov = hacCustom(E, max_lags);
           % vecm_est.Covariance = hac_cov;

            % Forecast + MSE
            [fcast, ~] = forecast(vecm_est, stepSize, test_data);
            MSE         = nanmean((fcast - test_data).^2, 'all');
            total_error = total_error + MSE;
            errors(end+1) = MSE;  %#ok<AGROW>
        end
    end

    % ===== Final metric calculation =====
    if count > 0
        rank_instability       = rank_changes / count;
        error_volatility       = std(errors) / mean(errors);
        normalized_MSE         = total_error / count;
        regime_stability_index = 0.7 * rank_instability + 0.3 * error_volatility;
    else
        % No windows had rank>0 => no valid VECM or no iteration
        rank_instability       = 0; 
        normalized_MSE         = 0;
        regime_stability_index = 0;
        warning('No valid (rank>0) windows found in runFixedVECMModel. All metrics set to 0.');
    end

    % Return metrics as a struct
    metrics.regimeStabilityIndex = regime_stability_index;
    metrics.rankInstability      = rank_instability;
    metrics.normalizedMSE        = normalized_MSE;
end
