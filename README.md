# Modeling Bitcoin Portfolio Risk & Hedging Effectiveness using Regime-Aware DCC-MIDAS-X Models

## 1. Project Name

Bachelor Graduation Thesis - June 2026 - Foreign Trade University, HCMC Campus 

Evaluation Score: Excellent (4.0/4.0)

## 2. Domain

* Financial Risk Management
* Advanced Econometrics Modelling
* Portfolio Hedging
* Statistic, Linear Algebra, R Statitical Programming

## 3. Served Stakeholder

Personal project for Bachelor Graduation

Academic Supervisor: Le Trung Thanh, PhD - Foreign Trade University, HCMC Campus

## 4. Timeline

Mar - June 2026

## 5. Scope of Work

* Model Bitcoin volatility and dynamic correlations across 3 asset classes: MSCI World, S&P GSCI, and PIMCO Investment Grade Corporate Bond Index.
* Examine the macro-financial impact of Global Economic Policy Uncertainty (GEPU).
* Develop advanced GARCH-MIDAS-X and DCC-MIDAS-X models incorporating structural breaks.
* Construct time-varying Optimal Hedge Ratios (OHR) for dynamic portfolio risk management.
* Evaluate hedging performance through in-sample, out-of-sample, and transaction-cost stress testing.

## 6. Process and Approach

Input: Daily asset prices and monthly GEPU data (Oct 2016 – Dec 2025).  
Output: Dynamic correlation estimates, time-varying hedge ratios, and hedging performance metrics.

Methodology and process:
* **Data Preparation:** Process daily financial log returns and log-differenced GEPU policy shock series.
* **Structural Breaks:** Apply Bai–Perron test to detect endogenous regime changes in policy uncertainty.
* **Volatility Modelling:** Estimate short- and long-run volatility components using GARCH-MIDAS-X.
* **Correlation Modelling:** Capture dynamic Bitcoin–asset joint distributions via regime-aware DCC-MIDAS-X.
* **Dynamic Hedging:** Generate time-varying Optimal Hedge Ratios (OHR) and hedged portfolio returns.
* **Model Evaluation:** Compare model selection and predictive superiority via AIC, BIC, Likelihood-Ratio, Diebold–Mariano, and Clark–West tests.
* **Robustness Testing:** Conduct rolling out-of-sample backtests, alternative train-test splits, and transaction-cost analysis (up to 50 bps).

## 7. Outcome

* Identified 3 significant GEPU structural break regimes: May 2019, Jan 2021, and Jul 2024.
* Quantified Bitcoin's volatility dynamics across turbulent macro-financial regimes and mapped how its dynamic financial relationships with traditional assets evolve in response to specific macroeconomic events. 
* Proved that Bitcoin's hedging relationship with traditional assets is strictly asset-specific and regime-dependent.
* Structural-break models consistently superior in model fit and out-of-sample forecasting accuracy.
* The full DCC-MIDAS-X + Structural Break model achieves the highest variance reduction in hedged portfolios:
  * BTC – MSCI World: 17.05%
  * BTC – S&P GSCI: 2.18%
  * BTC – PIMCO Bond: 3.03%
* Out-of-sample forecasting improvements are statistically significant for BTC–MSCI and BTC–GSCI pairs.
* Hedging performance remains robust after accounting for market transaction costs up to 50 bps.
* Demonstrates empirically that Bitcoin acts as a conditional hedge/diversifier rather than a universal safe haven.

## 8. My Role

* Project owner, Quantitative researcher & Engineer
* Conducting 100% all tasks independently (Research design, data processing, R numerical optimization, econometric modelling, dynamic hedging analysis, and thesis manuscript)

## 9. Lessons Learned

* **Econometric & Financial:** Advanced GARCH/DCC/MIDAS modeling, Quasi-Maximum Likelihood Estimation (QMLE), structural break testing, and dynamic portfolio backtesting
* **Engineering & Software:** Custom numerical optimization routines in R, end-to-end data pipelines, and robust statistical testing frameworks

## Prerequisites

R 4.5.x

## Disclaimer

This project is intended for educational purposes only.
