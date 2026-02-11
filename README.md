# Predictive Modeling and Bayesian Inference in R

Two applied statistical modeling exercises implemented from scratch in R using a reproducible RMarkdown workflow. The focus is on likelihood-based modeling, predictive evaluation, and simulation-based Bayesian inference.

---

## Part 1 — Predicting 3D Printed Object Weight

Modeled the relationship between CAD-estimated weight `x` and actual weight `y` using heteroscedastic Gaussian regression models.

### Mathematical Framework

Each observation is modeled as:

$begin:math:display$
Y\_i \\sim \\mathcal\{N\}\(\\mu\_i\, \\sigma\_i\^2\)
$end:math:display$

with:

$begin:math:display$
\\mu\_i \= \\beta\_1 \+ \\beta\_2 x\_i
$end:math:display$

Two alternative variance structures were implemented:

**Model A**

$begin:math:display$
\\sigma\_i\^2 \= \\exp\(\\beta\_3 \+ \\beta\_4 x\_i\)
$end:math:display$

**Model B**

$begin:math:display$
\\sigma\_i\^2 \= \\exp\(\\beta\_3\) \+ \\exp\(\\beta\_4\) x\_i\^2
$end:math:display$

Parameters were estimated via maximum likelihood by minimizing the negative log-likelihood.

### Predictive Evaluation

For each observation, predictive distributions were constructed and evaluated using:

- Leave-one-out cross-validation  
- Squared Error (point accuracy)  
- Dawid–Sebastiani score (distributional accuracy)  

A Monte Carlo test was used to assess whether one model provided significantly better predictive performance.

---

## Part 2 — Bayesian Estimation of Burial Counts

Estimated total population size `N` and recovery probability `φ` from observed femur counts.

### Model

$begin:math:display$
Y\_1\, Y\_2 \\sim \\text\{Binomial\}\(N\, \\phi\)
$end:math:display$

with priors:

$begin:math:display$
N \\sim \\text\{Geometric\}\(\\xi\)
$end:math:display$

$begin:math:display$
\\phi \\sim \\text\{Beta\}\(a\, b\)
$end:math:display$

Posterior expectations were approximated using Monte Carlo integration:

$begin:math:display$
\\hat\{E\}\[N \\mid y\] \=
\\frac\{\\sum\_k N\^\{\(k\)\} p\(y \\mid N\^\{\(k\)\}\, \\phi\^\{\(k\)\}\)\}
\{\\sum\_k p\(y \\mid N\^\{\(k\)\}\, \\phi\^\{\(k\)\}\)\}
$end:math:display$

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

## Structure

```
├── report.Rmd
├── code.R
├── report.html
└── README.md
```
