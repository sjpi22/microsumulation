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
- benchmarking: Folder for benchmarking functions in epi_functions.R
- tests: Folder for unit tests

Please cite the upcoming pre-print using:  
@misc{pi_calculating_2025,  
	title = {Calculating epidemiological outcomes from simulated longitudinal data},  
	author = {Pi, Selina and Goldhaber-Fiebert, Jeremy D and Alarid-Escudero, Fernando},  
	year = {2025},  
}