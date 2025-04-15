# Adaptive Regime Detection in Rolling-Window VECM Modeling
## 1. Affiliation
Project is conducted in collaboration with Professor Michael Robbins, Columbia University

## 2. Overview
This Project develops and tests a rolling-window Vector Error Correction Model (VECM) framework enhanced with threshold and Markov-switching capabilities to adaptively define the regime changes in extensive time-series data.

Johansen's cointegration with finite-sample correction method (bartlett correction) is being used to identify the long-term and stable relationships among non-stationary assets. The model then adjusts how these relationships are corrected (error-correction) depending on current regime.

We are testing threshold-VECMs (TVECMs) to see whether partitioning the error-correction term (ECT) space can capture nonlinear adjustments around the transaction cost boundary. Exploring the Markov-switching VECMs (MS-VECMs) to handle more gradual market transition through hidden state variables. Gridsearch parameter optimization is used to determine optimal, training-window length, reestimation frequency, cointegration rank and lag structure.

---------------------------------------------------------------------------------------------------

## License
Project affiliated with Professor Michael Robbins (Columbia University).

**No Use License**  
All rights to this project are reserved. You may not use, copy, modify, merge, publish, distribute, sublicense, or sell copies of this project or any portion thereof, unless explicitly authorized by the project owner. For inquiries about licensing, please contact michaelprobbins@aol.com.


### Key Files
- `runFixedVECMModel.m`  
  Implements a standard rolling-window approach:
  - windowSize = number of observations in each training window
  - stepSize = how far we move the window forward each iteration
  - Returns metrics:
    - `normalizedMSE`
    - `rankInstability`
    - `regimeStabilityIndex`

- `runAdaptiveVECMModel.m`  
  Implements an adaptive approach where the window size shrinks or grows based on observed volatility in the current window. Key parameters:
  - baseWindow, volMult, shrinkFact, growFact, etc.
  - Returns the same metrics as the fixed approach.

- `computeBartlett.m`  
  Provides small-sample corrections to the Johansen test statistics.

- `Experiment1Initialization2.mlx`  
  Loads `partition5082.mat` and prepares the `prices` variable for the experiment.

- `Experiment1Function1.mlx`  
  The main “training function” used by the Experiment Manager. It selects whether to call `runFixedVECMModel` or `runAdaptiveVECMModel` based on the parameter `modelType`.

---

## 2. Overview of the Adaptive Logic

- Measure Volatility in each training chunk:
  - If volatility exceeds `volMult * avgVol`, shrink the window (`baseWindow *= shrinkFact`).
  - If volatility is below `avgVol`, grow the window (`baseWindow *= growFact`).
- This attempts to shorten the window in high-volatility regimes for fresher data and lengthen it in calmer periods to capture a broader history.

---

## 3. Experiment Manager Usage

1. Open MATLAB and ensure all `.m` and `.mlx` files are in your MATLAB path.
2. Launch the Experiment Manager by typing:
   ```matlab
   experimentManager

