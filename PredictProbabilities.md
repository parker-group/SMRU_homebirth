# Predicted probabilities of home birth from the mixed-effects model

This document describes how predicted probabilities of home birth were generated from the final mixed-effects logistic regression model used in the analysis.

---

## 1. Model framework

Let $Y_{ij}$ denote the binary outcome for woman $i$ at outreach site $j$, where  
$Y_{ij} = 1$ indicates home birth and $Y_{ij} = 0$ indicates facility birth.

The fitted mixed-effects logistic regression model is

$$
\operatorname{logit}\!\left[ P(Y_{ij} = 1) \right]
\;=\;
\eta_{ij}
\;=\;
\mathbf{x}_{ij}^{\top} \boldsymbol{\beta} + u_j,
$$

where:

- $\mathbf{x}_{ij}$ is the vector of fixed-effect covariates for woman $i$ at site $j$  
  (including number of ANC visits, gestational age at enrollment, parity indicator, and other covariates),
- $\boldsymbol{\beta}$ is the vector of fixed-effect coefficients,
- $u_j \sim \mathcal{N}(0, \sigma^2)$ is the site-level random intercept,
- $\sigma^2$ is the estimated between-site variance.

---

## 2. Prediction target

We aim to estimate the **population-averaged predicted probability** of home birth as a function of the number of antenatal care (ANC) visits, marginalizing over outreach-site effects and holding other covariates at their observed values.

Predictions are generated:

- across the observed range of ANC visits,  
- for the full population,  
- and stratified by parity (primiparous vs multiparous).

---

## 3. ANC visit grid and standardization

Let $a$ denote a raw ANC visit value.

Predictions are generated on a grid:

$$
a \in \{ a_{\min}, a_{\min}+1, \dots, a_{\max} \},
$$

where $a_{\min}$ and $a_{\max}$ are the minimum and maximum observed ANC visits among analytic rows used in model fitting.

Because ANC visits were standardized before model fitting, raw grid values are converted to standardized units:

$$
z(a) \;=\; \frac{a - \mu_{\text{ANC}}}{s_{\text{ANC}}},
$$

where $\mu_{\text{ANC}}$ and $s_{\text{ANC}}$ are the mean and standard deviation of ANC visits in the analytic dataset.

---

## 4. Marginal standardization over covariates

Let $\mathbf{x}_i(a)$ denote the covariate vector for individual $i$, with:

- the ANC component replaced by $z(a)$,  
- all other covariates held at their observed values.

The fixed-effect linear predictor for individual $i$ at ANC value $a$ is:

$$
\eta_i(a) \;=\; \mathbf{x}_i(a)^{\top} \boldsymbol{\beta}.
$$

This corresponds to **marginal standardization** (predictive margins), because predictions are averaged over the empirical distribution of all other covariates.

---

## 5. Marginalization over site-level random effects

To obtain population-averaged probabilities, predictions are marginalized over the random intercept:

$$
P_i(a)
\;=\;
\mathbb{E}_{u}\!\left[
\operatorname{logit}^{-1}\!\left( \eta_i(a) + u \right)
\right],
\quad
u \sim \mathcal{N}(0, \sigma^2).
$$

This expectation has no closed-form expression and is approximated by Monte Carlo integration.

---

## 6. Monte Carlo approximation

For each ANC value $a$ and individual $i$:

1. Draw $K$ realizations of the random intercept:

$$
u^{(k)} \sim \mathcal{N}(0, \sigma^2),
\quad k = 1, \dots, K.
$$

2. Compute the predicted probability for each draw:

$$
p_i^{(k)}(a)
\;=\;
\operatorname{logit}^{-1}\!\left( \eta_i(a) + u^{(k)} \right).
$$

3. Average across draws:

$$
\widehat{P}_i(a)
\;=\;
\frac{1}{K} \sum_{k=1}^{K} p_i^{(k)}(a).
$$

---

## 7. Averaging over individuals

The predicted probability for ANC value $a$ is then obtained by averaging across all $N$ individuals in the prediction dataset:

$$
\widehat{P}(a)
\;=\;
\frac{1}{N} \sum_{i=1}^{N} \widehat{P}_i(a).
$$

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

- the fixed-effect coefficients $\boldsymbol{\beta}$, and  
- the site-level random-intercept variance $\sigma^2$.

---

### 9.1 Fixed-effect uncertainty

Fixed-effect coefficients are drawn from a multivariate normal distribution:

$$
\boldsymbol{\beta}^{(b)} \sim
\mathcal{N}\!\left( \boldsymbol{\beta}, \; \mathbf{V}_{\boldsymbol{\beta}} \right),
$$

where $\mathbf{V}_{\boldsymbol{\beta}}$ is the estimated variance–covariance matrix of $\boldsymbol{\beta}$.

---

### 9.2 Random-effect uncertainty

When available, uncertainty in the site-level standard deviation $\sigma$ is approximated by drawing from a log-normal distribution calibrated to the profile-likelihood confidence interval for $\sigma$. When a valid profile interval cannot be obtained, $\sigma$ is held fixed at its maximum-likelihood estimate.

---

### 9.3 Joint simulation

For each simulation draw $b = 1, \dots, B$:

1. Draw $\boldsymbol{\beta}^{(b)}$ and $\sigma^{(b)}$.  
2. Recompute predicted probabilities $\widehat{P}^{(b)}(a)$ using the same Monte Carlo marginalization over random effects.  

Pointwise 95% confidence intervals are then obtained as:

$$
\text{CI}_{0.95}(a)
\;=\;
\left[
\text{quantile}_{0.025}\!\left\{ \widehat{P}^{(b)}(a) \right\},
\;
\text{quantile}_{0.975}\!\left\{ \widehat{P}^{(b)}(a) \right\}
\right].
$$

---

## 10. Summary

Predicted probabilities of home birth are:

- population-averaged with respect to outreach-site effects,  
- standardized over the empirical distribution of all other covariates,  
- generated across the observed range of ANC visits, and  
- accompanied by simulation-based confidence intervals that account for uncertainty in both fixed and random effects.

This approach yields interpretable risk curves that reflect both individual-level predictors and residual between-site heterogeneity.
