# microsumulation
Microsumulation: A repo for calculating summary outcomes with microsimulation data and similar individual-level lifetime longitudinal data

Summary outcomes:
- Prevalence (age-specific)
- Incidence (age-specific)
- Age-conditional and lifetime risk
- Stage distribution of cancer at diagnosis
- Distribution of number of concurrent precancerous lesions (lesion multiplicity)

File directory:
- R: Folder for R functions
    - epi_functions.R: Functions for calculating summary outcomes
    - simulate_data.R: Functions for simulating data for unit tests and benchmarking
- vignettes.R: Vignettes for using the functions
- benchmarking: Folder for function benchmarking analyses
- tests: Folder for unit tests

### Citation

If you use this code or find it helpful, please cite our preprint using:  
> **Selina Pi, Jeremy D. Goldhaber-Fiebert, Fernando Alarid-Escudero.** *Calculating epidemiological outcomes from simulated longitudinal data.* medRxiv [Preprint]. May 2, 2025. [https://doi.org/10.1101/2025.04.30.25326766v1](https://doi.org/10.1101/2025.04.30.25326766v1)

or

```bibtex
@misc{pi_calculating_2025,
	title = {Calculating epidemiological outcomes from simulated longitudinal data},
	author = {Pi, Selina and Goldhaber-Fiebert, Jeremy D. and Alarid-Escudero, Fernando},
	year = {2025},
	publisher = {medRxiv [Preprint]},
	doi = {10.1101/2025.04.30.25326766},
	url = {https://www.medrxiv.org/content/10.1101/2025.04.30.25326766v1}
}
```
