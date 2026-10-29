
# CORSIKA Atmospheric Muon Flux Validation

Computational reconstruction and validation of atmospheric muon energy spectra and angular distributions using **CORSIKA 7.8050**, with comparisons against experimental measurements and established physical expectations.

This repository contains Python analysis tools developed as part of an undergraduate Physics thesis at **Universidad de los Andes (Bogotá, Colombia)**, investigating the thermal and electromagnetic modulation of atmospheric muons.

## 1. Overview

The objective of this project is to validate the physical behavior of atmospheric muons generated with CORSIKA before studying the influence of atmospheric electric fields on the muon flux at ground level.

Three complementary aspects of the simulation framework are considered in the broader research:

- **Energy-spectrum validation:** reconstruction of the atmospheric muon momentum spectrum and comparison with L3+C experimental measurements.
- **Angular-distribution validation:** verification of the expected zenith-angle dependence and azimuthal uniformity of atmospheric muons at sea level.
- **Electric-field implementation:** independent validation of the CORSIKA EFIELD configuration, performed as part of the broader simulation project.

This repository focuses on the energy-spectrum and angular-distribution analyses.

## 2. Simulation Framework

The simulations were performed using CORSIKA 7.8050 with EPOS and UrQMD as the high- and low-energy hadronic interaction models, respectively.

| Parameter | Configuration |
|---|---|
| Simulation software | CORSIKA 7.8050 |
| High-energy interaction model | EPOS |
| Low-energy interaction model | UrQMD |
| Atmospheric model | ATMOD 7 |
| Primary particles | Protons and helium, depending on the validation |
| Analysis language | Python 3 |
| Simulation output | CORSIKA binary particle files |

The energy-spectrum and angular-distribution validations use separate simulation campaigns and analysis configurations.

## 3. Energy-Spectrum Validation

### 3.1. Objective

The first validation investigates whether the simulated atmospheric muon spectrum reproduces the momentum dependence measured experimentally by the **L3+C experiment**.

L3+C provides a reference measurement of the atmospheric muon momentum spectrum, allowing the reconstructed simulation results to be evaluated against experimental observations.

### 3.2. Methodology

The analysis uses CORSIKA simulations of proton and helium primaries.

Primary-event information and secondary muon properties are extracted from the CORSIKA binary files. The resulting muon contributions are weighted according to the primary cosmic-ray flux parametrization of **Wiebel-Sooth**.

Contributions from heavier primary nuclei are approximated using a proton-superposition treatment.

The reconstructed muon spectrum is then compared with the L3+C experimental measurements.

The computational workflow is:

1. Extract primary cosmic-ray information from CORSIKA.
2. Extract secondary muons at the observation levels.
3. Apply the primary cosmic-ray spectral weights.
4. Reconstruct the atmospheric muon momentum spectrum.
5. Compare the simulated spectrum with L3+C measurements.

### 3.3. Results

The reconstructed spectrum reproduces the overall momentum dependence measured by L3+C.

Of the 13 experimental comparison points considered, **10 agree with the reconstructed spectrum within ±15%**.

These results support the ability of the simulation and reconstruction framework to reproduce the atmospheric muon momentum spectrum within the conditions studied.

### 3.4. Figures

**Reconstructed atmospheric muon spectrum and L3+C comparison**

[View energy-spectrum validation (PDF)](L3C_WS_FINAL_spectrum.pdf)

**Ratio between reconstructed and experimental spectra**

[View spectrum-ratio validation (PDF)](L3C_WS_FINAL_ratio.pdf)

### 3.5. Analysis Scripts

| Script | Description |
|---|---|
| `extract_corsika_primaries.py` | Extracts primary cosmic-ray properties from CORSIKA binary output |
| `extract_corsika_muons.py` | Extracts secondary muon properties at the observation levels |
| `compute_wiebel_sooth_muon_spectrum.py` | Reconstructs the muon spectrum using the Wiebel-Sooth primary flux parametrization |
| `compare_muon_spectrum_with_l3c.py` | Compares the reconstructed spectrum with L3+C experimental measurements |

## 4. Angular-Distribution Validation

### 4.1. Objective

The second validation investigates whether the simulated muons reproduce the expected angular characteristics of the atmospheric muon flux at sea level.

Two properties are examined:

- The zenith-angle dependence of atmospheric muons.
- The approximate uniformity of the azimuthal distribution.

### 4.2. Simulation Configuration

A dedicated simulation campaign was performed at sea level using proton primaries.

| Parameter | Value |
|---|---|
| Primary particle | Proton |
| Primary energy range | 10–100 GeV |
| Primary spectral index | −2.7 |
| Primary zenith-angle range | 0°–45° |
| Observation altitude | Sea level |
| Atmospheric model | ATMOD 7 |
| Number of simulated showers | 105,000 |
| Total extracted muons | 22,811 |
| Muons selected for angular analysis | 7,978 |

The selected muons satisfy:

$$
1 \leq p_\mu \leq 5\ \mathrm{GeV}/c,
$$

$$
0^\circ \leq \theta < 30^\circ.
$$

The complete dataset combines 21 independent simulation runs, each containing 5,000 primary showers.

### 4.3. Zenith-Angle Distribution

At sea level, the atmospheric muon intensity is commonly approximated by a dependence proportional to

$$
I(\theta)\propto\cos^2\theta
$$

within an appropriate momentum and angular range.

Including the angular measure, the corresponding zenith-angle distribution is

$$
\frac{dN}{d\theta}
\propto
\cos^2\theta\sin\theta.
$$

The reconstructed distribution was compared with this reference using six zenith-angle bins between 0° and 30°.

Statistical uncertainties were estimated from the contributions of individual simulated showers.

The comparison yielded:

$$
\chi^2 = 7.495,
\qquad
\nu = 5,
$$

$$
\chi^2_{\mathrm{red}} = 1.499,
\qquad
p = 0.1864.
$$

The simulated distribution is statistically compatible with the expected zenith-angle dependence.

**Zenith-angle validation**

[View zenith-angle validation (PDF)](zenith_validation_cos2.pdf)

### 4.4. Azimuthal Distribution

The azimuthal distribution was reconstructed using twelve equal-width bins covering the interval

$$
0^\circ \leq \phi < 360^\circ.
$$

The resulting distribution was compared with a constant reference, corresponding to azimuthal uniformity.

The goodness-of-fit test yielded:

$$
\chi^2 = 4.845,
\qquad
\nu = 11,
$$

$$
\chi^2_{\mathrm{red}} = 0.440,
\qquad
p = 0.9384.
$$

The simulated azimuthal distribution is statistically compatible with a constant.

An initial pilot simulation of 5,000 showers exhibited larger fluctuations and an apparent disagreement with uniformity. After increasing the sample to 105,000 showers, this discrepancy was no longer observed.

**Azimuthal validation**

[View azimuthal validation (PDF)](azimut_validation.pdf)

### 4.5. Statistical Summary

| Distribution | Reference | Reduced chi-square | p-value |
|---|---|---:|---:|
| Zenith | cos²θ sinθ | 1.499 | 0.1864 |
| Azimuth | Constant | 0.440 | 0.9384 |

Both distributions are statistically compatible with their respective reference models.

These results support the ability of the simulation framework to reproduce the expected angular behavior of atmospheric muons within the momentum and angular ranges considered.

### 4.6. Analysis Scripts

| Script | Description |
|---|---|
| `validate_corsika_angular_distribution.py` | Combines simulation runs, reconstructs angular distributions, estimates shower-level statistical uncertainties, and performs chi-square goodness-of-fit tests |
| `plot_muon_zenith_distribution.py` | Generates the final zenith-angle validation figure using the cos²θ sinθ reference |

## 5. Repository Structure

```text
corsika-muon-spectrum-analysis/
│
├── README.md
├── LICENSE
├── .gitignore
│
├── extract_corsika_primaries.py
├── extract_corsika_muons.py
├── compute_wiebel_sooth_muon_spectrum.py
├── compare_muon_spectrum_with_l3c.py
│
├── validate_corsika_angular_distribution.py
├── plot_muon_zenith_distribution.py
│
├── L3C_WS_FINAL_spectrum.pdf
├── L3C_WS_FINAL_ratio.pdf
├── zenith_validation_cos2.pdf
└── azimut_validation.pdf
```

## 6. Software Requirements

The analysis was developed in Python 3 using the following scientific libraries:

- NumPy
- Pandas
- SciPy
- Matplotlib

Install the dependencies with:

```bash
python3 -m pip install numpy pandas scipy matplotlib
```

The extraction scripts require CORSIKA binary simulation output.

The energy-spectrum comparison additionally requires the corresponding experimental reference dataset.

## 7. Computational Workflow

### Energy-spectrum reconstruction

```text
CORSIKA simulations
        |
        v
Primary and muon extraction
        |
        v
Primary cosmic-ray flux weighting
        |
        v
Muon momentum-spectrum reconstruction
        |
        v
Comparison with L3+C measurements
        |
        v
Energy-spectrum validation
```

### Angular-distribution reconstruction

```text
CORSIKA simulations
        |
        v
Primary and muon extraction
        |
        v
Combination of simulation runs
        |
        v
Muon selection and angular reconstruction
        |
        v
Statistical uncertainty estimation
        |
        v
Zenith and azimuth distributions
        |
        v
Chi-square goodness-of-fit tests
        |
        v
Angular-distribution validation
```

For the angular analysis, shower contributions are identified by both the simulation run and shower number, preventing different showers from being incorrectly combined when their local identifiers coincide.

## 8. Scope and Limitations

The results are specific to the primary-particle configurations, atmospheric model, observation altitudes, momentum intervals, and angular selections used in the simulations.

Agreement with L3+C measurements and with the expected angular distributions supports the physical consistency of the reconstructed muon flux.

However, these validations do not independently establish the accuracy of every aspect of the hadronic interaction models or atmospheric parametrization.

The repository contains analysis scripts and validation figures. CORSIKA executables, large binary simulation outputs, and externally sourced experimental datasets are not redistributed.

## 9. Research Context

This repository supports an undergraduate Physics thesis at **Universidad de los Andes**, focused on the thermal and electromagnetic modulation of the atmospheric muon flux.

The broader research investigates how atmospheric conditions, particularly electric fields associated with thunderstorms, may modify the flux of muons reaching ground level.

The energy-spectrum and angular-distribution validations establish a baseline for subsequent numerical investigations involving atmospheric electric fields and extended simulation campaigns.

## 10. References

1. Heck, D., Knapp, J., Capdevielle, J. N., Schatz, G., and Thouw, T. (1998). *CORSIKA: A Monte Carlo Code to Simulate Extensive Air Showers*. Forschungszentrum Karlsruhe, FZKA 6019.

2. Achard, P., et al. (L3 Collaboration). (2004). *Measurement of the atmospheric muon spectrum from 20 to 3000 GeV*. Physics Letters B, 598, 15–32.

3. Gaisser, T. K., Engel, R., and Resconi, E. (2016). *Cosmic Rays and Particle Physics*. 2nd edition. Cambridge University Press.

## Author

**Jacobo Trejos**

Undergraduate Physics  
Universidad de los Andes  
Bogotá, Colombia

Developed as part of undergraduate research in atmospheric muon physics.

