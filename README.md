# Semiparametric Regression for Misclassified Competing Risks Data

R code accompanying the paper:

*"Semiparametric Regression for Misclassified Competing Risks Data"*

This repository contains all functions and reproducible example code needed to:

1. Simulate competing risks data with time- and covariate-dependent misclassification.
2. Estimate misclassification probabilities using an external validation sample.
3. Fit the proposed semiparametric regression model using a B-spline-based pseudo-likelihood approach.
4. Conduct sensitivity analyses to assess robustness against transportability assumption violations.
5. Perform bootstrap inference for regression parameters while accounting for uncertainty in the estimated misclassification probabilities.
6. Compute and plot cumulative incidence functions (CIFs) under each sensitivity-analysis setting.

---

## **Repository contents**

### **Core estimation functions**

* **`bssmle.R`**
  Implements the proposed semiparametric regression estimator using B-splines and externally estimated misclassification probabilities. Returns baseline hazard spline coefficients and cause-specific regression effects.

* **`pseudo_likelihood_estimation_Mpofu.R`**
  Estimates time- and covariate-dependent misclassification probabilities using the pseudo-likelihood approach of Mpofu et al. (2020).

### **Simulation and example**

* **`simulate_data.R`**
  Generates competing-risks data with uni- or bidirectional misclassification, supporting both external validation and main-study simulation.

* **`example_analysis.Rmd`**
  A fully reproducible end-to-end example demonstrating:

  * Simulation of external validation data
  * Estimation of misclassification probabilities
  * Generation of the main analysis dataset
  * Application of the proposed semiparametric estimator
  * Sensitivity analysis using an η-grid
  * Bootstrap inference that accounts for uncertainty in the estimated misclassification probabilities
  * Reconstruction and plotting of CIFs

* **`example_analysis.html`**
  Rendered output of the complete example for easy viewing.

---

## Requirements

R (≥ 4.0) and the following packages:

```r
install.packages(c("alabama", "splines", "survival", "MASS", "boot", "sandwich",
                   "numDeriv", "Hmisc", "TeachingDemos", "lmtest"))
```

## Quick start

Open `example_analysis.Rmd` in RStudio and knit it. The example is
self-contained: it simulates an external validation dataset, fits the
misclassification model, generates the main analysis dataset, applies
the proposed semiparametric estimator, runs the η-grid sensitivity
analysis, performs bootstrap inference, and plots CIFs.

---

## **Method overview**

### **1. External validation and misclassification modeling**

Misclassification probabilities are estimated as:
$p_{21} = P(C^* = 2 \mid C = 1, T, Z)$ and
$p_{12} = P(C^* = 1 \mid C = 2, T, Z)$ using a double-sampling pseudo-likelihood approach (Mpofu et al., 2020).

The model uses:

* logistic regression for predictive values
* logistic regression for misclassification
* time- and covariate-dependent misclassification
* externally observed true causes

The estimated coefficients are later used to compute predicted misclassification probabilities in the main analysis dataset.

---

### **2. Semiparametric regression model**

The proposed estimator:

* uses B-splines to model baseline cumulative hazards
* incorporates external misclassification probabilities ($p_{12}, p_{21}$)
* fits a semiparametric proportional cause-specific hazards model
* simultaneously models both causes through a unified pseudo-likelihood
* produces bias-corrected cause-specific regression estimates

This method corrects for outcome misclassification using externally estimated misclassification probabilities.

---

### **3. Sensitivity analysis (η-shift)**

To examine violations of the transportability assumption, misclassification probabilities in the main dataset are adjusted via

$logit(p_{jh}(η)) = logit(p_{jh}) + η$, with $η \in$ {-0.5, -0.25, 0, 0.25, 0.5}.

Each modified misclassification scenario yields a new set of regression estimates, allowing assessment of robustness against transportability assumption violations.

---

### **4. Bootstrap inference**

Bootstrap samples are used to compute:

* bootstrap standard errors
* Wald statistics
* confidence intervals

The bootstrap procedure accounts for uncertainty in both the main analysis data and the estimated misclassification probabilities. In each bootstrap sample, the misclassification model parameters are regenerated using their estimated variance-covariance matrix, and the corresponding misclassification probabilities are recalculated before fitting the proposed model.

This provides finite-sample inference for regression effects.

---

### **5. CIF computation**

Cumulative incidence functions are reconstructed using:

* the spline-based baseline cumulative hazards
* the estimated cause-specific regression coefficients
* the resulting cause-specific hazards and survival function

---

## **Citation**

If you use this code in your research, please cite:

Balanos T, Yiannoutsos CT, Pabon-Rodriguez FM, Nan H, Bakoyannis G. *Semiparametric Regression for Misclassified Competing Risks Data*. arXiv preprint arXiv:2605.16652. 2026. https://doi.org/10.48550/arXiv.2605.16652

Code repository:
https://github.com/tbalanos/misclassified-competing-risks-bssmle

---

## **Contact**

**Theofanis Balanos, Ph.D.**  
Department of Biostatistics and Health Data Science
Richard M. Fairbanks School of Public Health
Indiana University Indianapolis 

Email: **tbalanos@iu.edu**

---
