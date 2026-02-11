# Predictive Modeling and Bayesian Inference in R

Two applied statistical modeling exercises implemented from scratch in R using a reproducible RMarkdown workflow. The focus is on likelihood-based modeling, predictive evaluation, and simulation-based Bayesian inference.

---

## Part 1 — Predicting 3D Printed Object Weight

Modeled the relationship between CAD-estimated weight $x$ and actual weight $y$ using heteroscedastic Gaussian regression models.

### Mathematical Framework

Each observation is modeled as:

$$
Y_i \sim \mathcal{N}(\mu_i, \sigma_i^2)
$$

with mean structure

$$
\mu_i = \beta_1 + \beta_2 x_i
$$

Two alternative variance structures were implemented.

**Model A**

$$
\sigma_i^2 = \exp(\beta_3 + \beta_4 x_i)
$$

**Model B**

$$
\sigma_i^2 = \exp(\beta_3) + \exp(\beta_4) x_i^2
$$

Parameters were estimated via maximum likelihood by minimizing the negative log-likelihood.

### Predictive Evaluation

For each observation, predictive distributions were constructed and evaluated using:

- Leave-one-out cross-validation  
- Squared Error (point accuracy)  
- Dawid–Sebastiani score (distributional accuracy)  

A Monte Carlo test was used to assess whether one model provided significantly better predictive performance.

---

## Part 2 — Bayesian Estimation of Burial Counts

Estimated total population size $N$ and recovery probability $\phi$ from observed femur counts.

### Model

$$
Y_1, Y_2 \sim \text{Binomial}(N, \phi)
$$

with priors

$$
N \sim \text{Geometric}(\xi)
$$

$$
\phi \sim \text{Beta}(a, b)
$$

Posterior expectations were approximated using Monte Carlo integration:

$$
\hat{E}[N \mid y] =
\frac{\sum_k N^{(k)} p(y \mid N^{(k)}, \phi^{(k)})}
{\sum_k p(y \mid N^{(k)}, \phi^{(k)})}
$$

All likelihood calculations were implemented using log-Gamma functions for numerical stability.

---

## Implementation Notes

- Custom log-likelihood functions  
- Direct maximum likelihood estimation  
- Predictive distribution construction  
- Manual leave-one-out cross-validation  
- Simulation-based Bayesian computation  
- Fully reproducible RMarkdown workflow  

---
