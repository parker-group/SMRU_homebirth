# Predicted probabilities of home birth from the mixed-effects model

This document describes how predicted probabilities of home birth were generated from the final mixed-effects logistic regression model used in the analysis.

---

## 1. Model framework

Let `Y_ij` denote the binary outcome for woman `i` at outreach site `j`, where  
`Y_ij = 1` indicates home birth and `Y_ij = 0` indicates facility birth.

The fitted mixed-effects logistic regression model is:

    logit( P(Y_ij = 1) ) = η_ij = x_ij^T β + u_j

where:

- `x_ij` is the vector of fixed-effect covariates for woman `i` at site `j`
  (including number of ANC visits, gestational age at enrollment, parity indicator, and other covariates),
- `β` is the vector of fixed-effect coefficients,
- `u_j ~ Normal(0, σ^2)` is the site-level random intercept,
- `σ^2` is the estimated between-site variance.

---

## 2. Prediction target

We aim to estimate the **population-averaged predicted probability** of home birth as a function of the number of antenatal care (ANC) visits, marginalizing over outreach-site effects and holding other covariates at their observed values.

Predictions are generated:

- across the observed range of ANC visits,
- for the full population,
- and stratified by parity (primiparous vs multiparous).

---

## 3. ANC visit grid and standardization

Let `a` denote a raw ANC visit value.

Predictions are generated on a grid:

    a ∈ {a_min, a_min + 1, ..., a_max}

where `a_min` and `a_max` are the minimum and maximum observed ANC visits among analytic rows used in model fitting.

Because ANC visits were standardized before model fitting, raw grid values are converted to standardized units:

    z(a) = (a - μ_ANC) / s_ANC

where `μ_ANC` and `s_ANC` are the mean and standard deviation of ANC visits in the analytic dataset.

---

## 4. Marginal standardization over covariates

Let `x_i(a)` denote the covariate vector for individual `i`, with:

- the ANC component replaced by `z(a)`,
- all other covariates held at their observed values.

The fixed-effect linear predictor for individual `i` at ANC value `a` is:

    η_i(a) = x_i(a)^T β

This corresponds to **marginal standardization** (predictive margins), because predictions are averaged over the empirical distribution of all other covariates.

---

## 5. Marginalization over site-level random effects

To obtain population-averaged probabilities, predictions are marginalized over the random intercept:

    P_i(a) = E_u [ logit^{-1}( η_i(a) + u ) ],   u ~ Normal(0, σ^2)

This expectation has no closed-form expression and is approximated by Monte Carlo integration.

---

## 6. Monte Carlo approximation

For each ANC value `a` and individual `i`:

1. Draw `K` realizations of the random intercept:

       u^(k) ~ Normal(0, σ^2),   k = 1, ..., K

2. Compute the predicted probability for each draw:

       p_i^(k)(a) = logit^{-1}( η_i(a) + u^(k) )

3. Average across draws:

       P̂_i(a) = (1/K) * sum_{k=1}^K p_i^(k)(a)

---

## 7. Averaging over individuals

The predicted probability for ANC value `a` is then obtained by averaging across all `N` individuals in the prediction dataset:

    P̂(a) = (1/N) * sum_{i=1}^N P̂_i(a)

This yields a **population-averaged predicted probability curve** for home birth as a function of ANC attendance.

---

## 8. Parity-stratified predictions

To generate parity-specific curves:

- the same procedure is applied using only rows for primiparous women,
- and separately using only rows for multiparous women.

All covariates other than ANC visits are again held at their observed values within each subgroup.

---

## 9. Uncertainty estimation

Uncertainty in predicted probabilities is quantified by jointly propagating uncertainty in:

- the fixed-effect coefficients `β`, and
- the site-level random-intercept variance `σ^2`.

---

### 9.1 Fixed-effect uncertainty

Fixed-effect coefficients are drawn from a multivariate normal distribution:

    β^(b) ~ Normal( β, V_β )

where `V_β` is the estimated variance–covariance matrix of `β`.

---

### 9.2 Random-effect uncertainty

When available, uncertainty in the site-level standard deviation `σ` is approximated by drawing from a log-normal distribution calibrated to the profile-likelihood confidence interval for `σ`. When a valid profile interval cannot be obtained, `σ` is held fixed at its maximum-likelihood estimate.

---

### 9.3 Joint simulation

For each simulation draw `b = 1, ..., B`:

1. Draw `β^(b)` and `σ^(b)`.
2. Recompute predicted probabilities `P̂^(b)(a)` using the same Monte Carlo marginalization over random effects.

Pointwise 95% confidence intervals are then obtained as:

    CI_0.95(a) = [
      quantile_0.025 { P̂^(b)(a) },
      quantile_0.975 { P̂^(b)(a) }
    ]

---

## 10. Summary

Predicted probabilities of home birth are:

- population-averaged with respect to outreach-site effects,
- standardized over the empirical distribution of all other covariates,
- generated across the observed range of ANC visits, and
- accompanied by simulation-based confidence intervals that account for uncertainty in both fixed and random effects.

This approach yields interpretable risk curves that reflect both individual-level predictors and residual between-site heterogeneity.
