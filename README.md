# Cup-to-Disc Ratio as an Effect Modifier of GLP-1-Associated NAION Risk
### A Simulation Study

**Author:** Holli Groff, OD, MPH  
**Date:** October 2026  
**Language:** R

---

## Background

Non-arteritic anterior ischemic optic neuropathy (NAION) is a leading cause of sudden, permanent vision loss in adults over 50, caused by infarction of the optic nerve head. Recent observational studies and meta-analyses have identified a statistically significant association between GLP-1 receptor agonists (particularly semaglutide) and increased NAION risk, with pooled hazard ratios of approximately 2.0–2.2 compared to other antidiabetic medications (PLOS Medicine, 2025).

However, a critical and unexplored question remains: **does pre-existing optic disc anatomy modify this risk?**

The leading structural hypothesis for NAION is the "disc at risk" — a small, crowded optic disc with a low cup-to-disc ratio (CDR) that leaves little anatomical reserve when perfusion pressure drops. If GLP-1 agonists reduce optic nerve perfusion (through blood pressure lowering, altered CSF dynamics, or direct vascular effects), patients with small CDR may be disproportionately vulnerable compared to those with normal disc anatomy.

This is a **novel research question** that no published study has addressed. As an optometrist with 26 years of clinical experience, the author identified this gap from direct clinical observation of NAION patients and knowledge of optic disc anatomy.

---

## Research Question

> Does pre-existing cup-to-disc ratio (CDR) modify the risk of NAION in patients taking a GLP-1 receptor agonist?

This is an **effect modification** (interaction) question — CDR as a biological moderator of the GLP-1/NAION relationship.

---

## Study Design

**Type:** Simulation study to establish methodology and estimate sample size requirements for a proposed real-world study

**Simulated population:** 50,000 adults with type 2 diabetes or obesity  
**Follow-up:** 3 years  
**Outcome:** NAION diagnosis  
**Exposure:** GLP-1 receptor agonist use  
**Effect modifier:** Cup-to-disc ratio (CDR), categorised as:
- Small / disc at risk: CDR ≤ 0.30
- Normal: CDR 0.31–0.60
- Large: CDR > 0.60

### Simulation parameters (grounded in published literature)

| Parameter | Value | Source |
|---|---|---|
| Background NAION incidence | 18.9 per 100,000/year | Liu et al., *Ophthalmology* 2025 |
| GLP-1 main effect HR | 2.17 | PLOS Medicine meta-analysis, 2025 |
| CDR distribution | Mean 0.40, SD 0.12 | Jonas et al., ophthalmology literature |
| Small CDR prevalence | ~15–20% | Population estimates |
| Interaction OR (hypothesized) | 3.0 | **Novel — the primary estimand** |

---

## Methods

1. **Population simulation** using biologically plausible distributions for age, sex, GLP-1 exposure, CDR, hypertension, and smoking
2. **Outcome simulation** via logistic model incorporating main effects and a GLP-1 × CDR interaction term
3. **Multivariable logistic regression** with interaction term (`glm`, binomial family)
4. **Stratified analysis** — crude odds ratios for GLP-1 within each CDR stratum
5. **Sample size calculation** for 80% power to detect interaction OR ≥ 3.0 (α = 0.05)

### Key model term
```r
glm(naion ~ glp1 + cdr_small + cdr_large +
      I(glp1 * cdr_small) +   # interaction term
      age_c + htn + smoker,
    family = binomial(link = "logit"))
```

---

## Key Results (Simulated)

| CDR Category | GLP-1 NAION rate | Non-GLP-1 rate | Crude OR |
|---|---|---|---|
| Small (≤0.30) | ~1,573 per 100,000 | ~144 per 100,000 | **11.10** (5.34–23.08) |
| Normal (0.31–0.60) | ~194 per 100,000 | ~92 per 100,000 | **2.10** (1.20–3.70) |

**Interaction term (GLP-1 × small CDR): OR 5.50 (95% CI 2.35–12.90), p < 0.001**

The diverging NAION rates between CDR strata — visible in the interaction plot — are the hallmark of effect modification.

**Estimated sample size for a real-world study:** ~31,000 patients, expecting ~160 NAION events. Feasible with large ophthalmic registries (e.g., IRIS Registry) or linked EHR + claims data.

---

## Outputs

Running `naion_cdr_simulation.R` produces three figures:

| File | Description |
|---|---|
| `naion_cdr_stratified_forest.png` | Forest plot of GLP-1 OR within each CDR stratum |
| `naion_cdr_interaction_plot.png` | Interaction plot showing diverging NAION rates |
| `naion_cdr_regression.png` | Multivariable logistic regression OR plot |

---

## Requirements

```r
install.packages(c("dplyr", "ggplot2", "patchwork", "broom"))
```

R ≥ 4.0. All packages installed automatically on first run if missing.

---

## Clinical Implications

If confirmed in real-world data, these findings would suggest:
- Patients with small CDR ("disc at risk") may warrant ophthalmology consultation before initiating GLP-1 therapy
- CDR screening could identify a high-risk subgroup for closer monitoring
- The GLP-1/NAION association may be largely driven by anatomically susceptible patients, rather than representing uniform risk

---

## Proposed Real-World Study

The ideal dataset to test this hypothesis would link:
- **Pharmacy/claims data** (GLP-1 exposure, diabetes/obesity diagnoses, NAION diagnosis codes)
- **Ophthalmic imaging records** (CDR from fundus photography or OCT)

Suitable data sources: **IRIS Registry** (American Academy of Ophthalmology), **VA Healthcare System**, or a linked EHR + claims dataset with ophthalmic data.

---

## References

- Liu et al. (2025). *Ophthalmology*. NAION incidence in semaglutide users vs comparators.
- PLOS Medicine meta-analysis (2025). Semaglutide-associated NAION risk: systematic review of 5 observational studies. Pooled HR 2.17 (1.73–2.74).
- Hathaway et al. (2024). *JAMA Ophthalmology*. Semaglutide and NAION risk in T2D and obesity cohorts.
- Amini et al. (2025). *Acta Ophthalmologica*. Review of semaglutide/NAION association.

---

## About the Author

Holli Groff is an optometrist (OD) with an MPH in Epidemiology and Medical Analytics from King's College London. This project sits at the intersection of her 26 years of ophthalmic clinical experience and her training in real-world evidence, pharmacovigilance, and health data science. She is transitioning into life sciences industry roles in RWE, HEOR, and medical affairs.

🔗 [LinkedIn](https://www.linkedin.com/in/holligroff)
