# Ambient Noise Seismology: From Theory to Practice

A comprehensive, semester-long Jupyter notebook course on **ambient noise seismology**, covering theory, mathematics, and hands-on Python coding.

## 📖 [View the Course Online](https://miichenry.github.io/ambient-noise-seismology-course/)

## Course Content

| Module | Topic |
|--------|-------|
| 1 | Introduction to Seismic Ambient Noise |
| 2 | Seismic Data Handling with ObsPy |
| 3 | The Cross-Correlation Theorem |
| 4 | Green's Function Retrieval from Ambient Noise |
| 5 | Ambient Noise Preprocessing |
| 6 | Cross-Correlation Computation |
| 7 | Stacking Methods |
| 8 | Optimal Processing of Noise Correlations |
| 9 | Surface Wave Dispersion Analysis (FTAN) |
| 10 | Ambient Noise Tomography: Inversion |
| 11 | Seismic Interferometry Principles |
| 12 | Seismic Velocity Changes (dv/v) Monitoring |
| 13 | HVSR / MHVSR Method |
| 14 | NoisePy: Large-Scale Processing |
| 15 | Applications and Case Studies |
| 16 | Dense Arrays, Advanced Topics, and Exercises |

## Features

- **Theory + Practice**: Each module pairs mathematical foundations with interactive Python code
- **30+ citations** from the ambient noise seismology literature
- **Working implementations** of key algorithms: temporal normalization, spectral whitening, cross-correlation, stacking (linear, PWS, robust), FTAN, stretching dv/v, HVSR, 2D tomographic inversion
- **Real data examples** using IRIS FDSN web services (ObsPy)
- **Synthetic demonstrations** to build physical intuition
- **Semester exercises** for coursework

## Key References

- Bensen et al. (2007, 2008) — Processing pipeline and US tomography
- Campillo & Roux (2015) — Theoretical foundations
- Fichtner et al. (2017, 2020) — Optimal processing framework
- Jiang & Denolle (2020) — NoisePy
- Cox et al. (2020) — HVSR lognormal statistics
- Cabrera-Pérez et al. (2023) — La Palma ANT
- Ryberg et al. (2022) — LARGE-N mineral exploration
- Stehly et al. (2024) — Noise source dynamics

## Requirements

```bash
pip install obspy numpy scipy matplotlib
```

## Usage

```bash
jupyter notebook Ambient_Noise_Seismology_Course.ipynb
```

## License

This course material is provided for educational purposes.
