# **Semiparametric Regression for Misclassified Competing Risks Data Using External Validation**

R code accompanying my dissertation paper:
**“Semiparametric Regression Analysis of Misclassified Competing Risks Data Using External Validation.”**

This repository contains all functions and reproducible example code needed to:

1. Simulate competing-risks data with time- and covariate-dependent misclassification
2. Estimate time- and covariate-dependent misclassification probabilities using an external validation sample
3. Fit the proposed semiparametric regression model using a B-spline pseudo-likelihood
4. Conduct sensitivity analysis for violations of the transportability assumption
5. Perform independent-subject bootstrap inference for regression parameters
6. Compute and plot cumulative incidence functions (CIFs) under each sensitivity setting

---

## **Repository Contents**

### **Core Estimation Functions**

* **`bssmle.R`**
  Implements the proposed semiparametric regression estimator using B-splines and externally estimated misclassification probabilities. Returns baseline hazard spline coefficients and cause-specific regression effects.

* **`pseudo_likelihood_estimation_Mpofu.R`**
  Pseudo-likelihood estimation of time- and covariate-dependent misclassification probabilities following Mpofu et al. (2020).

### **Simulation and Example**

* **`simulate_data.R`**
  Generates competing-risks data with uni- or bidirectional misclassification, supporting external validation and main-study simulation.

* **`example_analysis.Rmd`**
  A fully reproducible end-to-end example demonstrating:

  * Simulation of external validation data
  * Estimation of misclassification probabilities
  * Main analysis dataset generation
  * Application of the proposed semiparametric estimator
  * Sensitivity analysis using an $η$-grid
  * Independent-subject bootstrap inference
  * Reconstruction and plotting of CIFs

* **`example_analysis.html`**
  Rendered output of the complete example for easy viewing access.

---

## **Method Overview**

### **1. External Validation & Misclassification Modeling**

Misclassification probabilities are estimated as:
$p_{21} = P(C^* = 2 \mid C = 1, T, Z)$, and $p_{12} = P(C^* = 1 \mid C = 2, T, Z)$ using a double-sampling pseudo-likelihood approach (Mpofu et al., 2020).

The model uses:
* logistic regression for predictive values
* logistic regression for misclassification
* time- and covariate-dependent misclassification
* externally observed true causes

The estimated coefficients are later used to compute predicted misclassification probabilities in the main analysis dataset.

---

### **2. Semiparametric Regression Model**

The proposed estimator:

* uses B-splines to model baseline cumulative hazards
* incorporates external misclassification probabilities ($p_{12}, p_{21}$)
* fits a semiparametric proportional cause-specific hazards model
* simultaneously models both causes through a unified pseudo-likelihood
* produces bias-corrected cause-specific regression estimates

This method corrects for outcome misclassification using externally estimated misclasification probabilities.

---

### **3. Sensitivity Analysis (η-Shift)**

To examine violations of the transportability assumption, misclassification probabilities in the main dataset are adjusted via:

$logit(p_{jh}(η)) = logit(p_{jh}) + η$ with $η$ ∈ {-0.5, -0.25, 0, 0.25, 0.5}.

Each modified misclassification scenario yields a new set of regression estimates, allowing assessment of robustness to misclassification misspecification.

---

### **4. Independent-Subject Bootstrap**

Independent subjects are resampled with replacement to compute:

* bootstrap standard errors
* Wald statistics
* confidence intervals

This provides finite-sample inference for regression effects.

---

### **5. CIF Computation**

Cumulative incidence functions are reconstructed using:

* the spline-based baseline cumulative hazards
* the estimated cause-specific regression coefficients
* the resulting cause-specific hazards and survival function

---
Here is the **correct, fully adjusted, Paper 1 version** — simple, accurate, and matching your repo structure and filenames.

---

## **How to Run the Example**

### **1. Clone or download the repository**

```bash
git clone https://github.com/tbalanos/misclassified-competing-risks-bssmle.git
```

### **2. Open R or RStudio and source the scripts**

```r
source("simulate_data.R")
source("pseudo_likelihood_estimation_Mpofu.R")
source("bssmle.R")
```

### **3. Run the full analysis**

Open:

```
example_analysis.Rmd
```

and click **Knit**, or view the pre-rendered:

```
example_analysis.html
```

---

## **Citation**

If you use this code in your research, please cite:

**Balanos, T. (2026).**
*Semiparametric Regression Analysis of Misclassified Competing Risks Data Using External Validation.*
Indiana University Indianapolis.

(Full journal citation will be added once published.)

---

## **Contact**

**Theofanis Balanos**
Department of Biostatistics & Health Data Science
Indiana University Indianapolis

Email: **[tbalanos@iu.edu](mailto:tbalanos@iu.edu)**

---




