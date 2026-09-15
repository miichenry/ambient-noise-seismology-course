# Ambient Noise Seismology — Lecture Notes

## Course Overview

**Format:** 16 modules, designed for a semester (14-16 weeks). Each module pairs theory (markdown cells) with hands-on practice (code cells). The notebook is self-contained — students run cells sequentially in Jupyter.

**Prerequisites:** Basic seismology (wave equation, seismograms), Python fundamentals (NumPy, Matplotlib). No prior experience with ObsPy or ambient noise required.

**Logistics:** Students need a working Python 3 environment with `obspy`, `numpy`, `scipy`, `matplotlib`. Some cells download data from IRIS FDSN — internet access required for those sessions. The synthetic exercises work offline.

---

## Week 1 — Module 1: Introduction to Seismic Ambient Noise
**Cells: 0–8 | Time: ~2 hours (1h lecture + 1h practice)**

### Key Points to Convey

1. **The paradigm shift.** Traditional seismology waits for earthquakes. Ambient noise seismology uses what used to be thrown away — the "noise" between events. This is a revolution because:
   - No need for earthquakes (works in aseismic regions)
   - Continuous monitoring possible (not event-dependent)
   - Works at all scales: 10 m to 1000 km

2. **What generates the noise.** Draw a frequency axis on the board:
   - **> 30 s**: Earth's hum — ocean infragravity waves, atmospheric coupling
   - **10-20 s**: Primary microseism — direct ocean swell hitting the coast (same frequency as ocean waves)
   - **3-10 s**: Secondary microseism — **the main one** — wave-wave interaction of opposing ocean swells. Explain Longuet-Higgins (1950): two opposing waves create a standing wave with pressure oscillations at double the ocean frequency. This is the dominant noise source globally.
   - **< 1 s**: Cultural noise — traffic, industry, wind. Dominant in cities.

3. **Stationarity matters.** Introduce the idea early: the cross-correlation theory assumes isotropic, stationary noise. In reality:
   - Noise sources are seasonal (North Atlantic storms in winter)
   - Mediterranean stations see non-stationary events from the Adriatic/Aegean (Stehly et al., 2024)
   - This is why we need **months** of data — averaging over seasons approximates isotropy

### Practice Session

- **Practice 1.1**: Download real data, compute PSD, identify the microseismic peaks. Let students compare stations at different locations — coastal vs. inland, quiet vs. noisy.
- Point out the Peterson NLNM/NHNM models. Ask: "Is this station noisy or quiet?"
- **Tip**: If IRIS is slow, have a pre-downloaded dataset as backup.

### Homework Suggestion
Download 24h of data from a station near the coast and one far inland. Compare PSDs. Where are the microseismic peaks stronger? Why?

---

## Week 2 — Module 2: Seismic Data Handling with ObsPy
**Cells: 9–12 | Time: ~2 hours**

### Key Points

1. **ObsPy is the lingua franca.** Every ambient noise code (NoisePy, MSNoise, SeisNoise) builds on ObsPy. Students must be fluent in:
   - `UTCDateTime` — all times are UTC, arithmetic with seconds
   - `Trace` — single channel: `.data` (numpy array), `.stats` (metadata)
   - `Stream` — collection of Traces: `st.select()`, `st.filter()`, `st.merge()`
   - `Client("IRIS")` — FDSN data access: `get_waveforms()`, `get_stations()`

2. **SEED naming convention**: Network.Station.Location.Channel (e.g., `IU.ANMO.00.BHZ`)
   - Channel codes: band (B=broad, H=high), instrument (H=high-gain), orientation (Z/N/E)
   - Location code `00` vs `10` — different sensors at same station

3. **Station pair for CC**: We use ANMO-CCM (~1400 km). Explain why this distance:
   - Far enough for surface waves to develop
   - Close enough that travel time (~400 s at 3.5 km/s) fits in a reasonable lag window
   - Both are GSN stations with reliable, continuous data

### Practice Session

- **Practice 2.1**: Create Traces by hand, manipulate stats, basic operations
- **Practice 2.2**: Download ANMO-CCM data, plot waveforms
- Let students try different stations — suggest they pick stations from their home region

### Common Pitfall
Students may get confused by different sampling rates. Emphasize that both traces must be at the same sample rate before cross-correlation. ObsPy's `resample()` or `decimate()` handles this.

---

## Week 3 — Module 3: The Cross-Correlation Theorem
**Cells: 13–17 | Time: ~3 hours (theory-heavy)**

### This is the Core of the Course

1. **Start simple (Practice 3.1).** Two Ricker wavelets with a known delay. Cross-correlate. The peak appears at the delay time. This is just signal processing — nothing specific to seismology yet.

   **Write on the board:**
   $$C_{AB}(\tau) = \int_0^T A(t) \cdot B(t+\tau) \, dt$$
   
   And the frequency-domain equivalent (the theorem):
   $$C_{AB}(f) = \hat{U}_A^*(f) \cdot \hat{U}_B(f)$$
   
   Emphasize: this is why it's fast — FFT is O(N log N) vs. O(N^2) in time domain.

2. **The magic (Practice 3.2).** Now the paradigm-changing experiment: cross-correlate *noise* from two receivers surrounded by random sources. **Peaks emerge at the inter-station travel time.** This should produce a genuine "wow" moment.

   **Explain the physics step by step:**
   - Each source emits noise that reaches both stations with different delays
   - The cross-correlation measures the *differential* delay for each source
   - Sources along the station-pair axis (the "end-fire lobes") produce a consistent differential delay = inter-station travel time
   - Sources perpendicular produce zero differential delay
   - Summing over many sources: end-fire contributions add coherently, everything else averages to zero

3. **End-fire lobes (Practice 3.3).** Three scenarios:
   - Isotropic sources → symmetric CC (both causal and acausal peaks)
   - One-sided sources → asymmetric CC (only one peak)
   - Perpendicular sources → no peak at the expected time
   
   **This directly illustrates why stacking works**: different days sample different source directions. Over months, you approximate isotropy.

### Pedagogical Tip
The synthetic simulations use frequency-domain phase shifts (`exp(-2j*pi*f*tau)`) for accurate delays. If students ask why not just shift in time domain — explain that integer-sample shifts truncate the signal and lose the coherent delay structure. Phase shifts are exact.

---

## Week 4 — Module 4: Green's Function Retrieval
**Cells: 18 | Time: ~1.5 hours (pure theory)**

### The Mathematical Foundation

1. **Three approaches to the same result** (give all three — students have different intuitions):

   - **Time-reversal analogy** (heuristic): Correlation is like a time-reversal experiment. Recording noise at A, time-reversing it, and re-emitting it — the wavefield focuses at B. The correlation of the two recordings captures this focusing.

   - **Stationary phase** (geometric): Only noise sources in the end-fire lobes contribute. The stationary phase condition selects paths along the inter-station great circle. The CC is the sum of Green's functions for those paths.

   - **Representation theorem** (rigorous): Wapenaar (2004) showed formally that for an isotropic noise field:
     $$\langle u_1(\omega) u_2^*(\omega) \rangle \propto \text{Im}[G(\vec{r}_1, \vec{r}_2; \omega)]$$
     
     The CC of the noise field equals the imaginary part of the Green's function.

2. **The Bessel function connection** (Aki, 1957): For a 2D isotropic field of plane waves:
   $$C_{12}(\omega) \propto J_0(k|\vec{r}_1 - \vec{r}_2|)$$
   
   This was known since 1957 but nobody connected it to Green's function retrieval until Lobkis & Weaver (2001) and Campillo & Paul (2003).

3. **Historical timeline** — draw on the board:
   - 1957: Aki — spatial autocorrelation = Bessel function
   - 1968: Claerbout — daylight imaging concept
   - 2001: Lobkis & Weaver — proof in ultrasonics
   - 2003: Campillo & Paul — first seismic demonstration (coda waves)
   - 2004: Shapiro & Campillo — broadband Rayleigh waves from noise
   - 2005: Shapiro et al. — first ANT (California)

### Key Assumptions to Emphasize
- **Diffuse wavefield** (isotropy OR equipartition)
- Green's function retrieval is an *asymptotic* result — requires infinite sources or infinite time
- In practice, we get an *approximation* — and the quality depends on source distribution and stacking duration

---

## Week 5 — Module 5: Ambient Noise Preprocessing
**Cells: 19–22 | Time: ~3 hours**

### The Bensen et al. (2007) Recipe

This paper defined the field. Walk through the four phases:

**Phase 1 — Single-station preparation:**
1. Remove instrument response
2. Demean, detrend
3. Bandpass filter
4. Cut to 1-day segments
5. **Temporal normalization** (the key step)
6. **Spectral whitening**

**Phase 2:** Cross-correlate and stack (Modules 6-7)
**Phase 3:** Dispersion measurement (Module 9)
**Phase 4:** Quality control

### Temporal Normalization (Practice 5.1)

Why? Earthquakes and transients dominate the amplitude of raw recordings. We need to suppress them while preserving the phase information in the ambient noise.

**Five methods** (show all, explain trade-offs):

| Method | How it works | Pro | Con |
|--------|-------------|-----|-----|
| One-bit | `sign(data)` | Simple, robust | Destroys all amplitude info |
| Clipping | Threshold at RMS | Preserves some amplitude | Arbitrary threshold |
| Event removal | Zero out detected EQs | Clean | May miss small events |
| RAM | Divide by running mean of |data| | **Best balance** — recommended | Needs tuning of window width |
| Water-level | Iterative down-weighting | Gentle | Slower |

**RAM equation** (Bensen et al., 2007, Eq. 1):
$$w_n = \frac{1}{2N+1} \sum_{j=n-N}^{n+N} |d_j|, \qquad \tilde{d}_n = d_n / w_n$$

**Practical tip**: The RAM window width (2N+1) should be about half the minimum period of interest. For microseismic band (7-20 s), use N ~ 50-100 samples at 1 Hz.

### Spectral Whitening (Practice 5.2)

Why? The noise spectrum is peaked at the microseismic frequencies. Without whitening, the CC is dominated by 7 s energy and you get poor resolution at other periods.

How: Divide the spectrum by its smoothed amplitude, preserving only the phase.

**Show the "before/after"**: the peaked spectrum becomes flat, and the CC becomes broadband.

### Practice 5.3: Complete Pipeline

Students implement the full Bensen et al. (2007) Phase 1 pipeline. This is the function they'll reuse throughout the rest of the course.

---

## Week 6 — Module 6: Cross-Correlation Computation
**Cells: 23–26 | Time: ~2.5 hours**

### From Preprocessed Data to CC

1. **Windowed CC**: Don't correlate the entire day at once. Cut into windows (e.g., 30 min), correlate each, then stack. This allows:
   - Quality control (reject bad windows)
   - Sub-stacking for temporal monitoring
   - Memory management for large datasets

2. **CC methods** (Practice 6.2, NoisePy Table 1):

   | Method | Formula | When to use |
   |--------|---------|-------------|
   | Pure xcorr | $\hat{U}_A^* \hat{U}_B$ | Standard, most common |
   | Coherency | $\hat{U}_A^* \hat{U}_B / |\hat{U}_A||\hat{U}_B|$ | Equalized, good for heterogeneous noise |
   | Deconvolution | $\hat{U}_A^* \hat{U}_B / |\hat{U}_A|^2$ | Source-normalized |
   | Phase CC | $e^{i\phi_A^*} e^{i\phi_B}$ | Pure phase, most aggressive |

3. **Limitations of our computation.** Be upfront:
   - 1 day of data is **not enough** — we're illustrating the workflow, not producing publishable results
   - Bensen et al. (2007): SNR ~ T^{1/n}, n ≈ 2.5-3. Need months for convergence.
   - The filtered stacked CC may not show a clear surface wave arrival — explain this is expected and why (insufficient stacking, non-isotropic noise on a single day, long station distance)
   - For FTAN and subsequent analysis, we switch to synthetic data

---

## Week 7 — Module 7: Stacking Methods
**Cells: 27–29 | Time: ~2 hours**

### Why Stack?

**Connect back to the theory from Week 3:**
- The Green's function retrieval theorem has an ensemble average ⟨...⟩
- At any instant, noise is NOT isotropic
- Stacking = practical implementation of the ensemble average
- Ergodicity assumption: time average ≈ ensemble average

**Two assumptions of stacking:**
1. **Ergodicity** — over time, different source configurations sample enough azimuths
2. **Stationarity of the medium** — the Earth doesn't change during the stacking window

These are in tension for monitoring: long stacks → better SNR but assumes frozen medium.

### Three Stacking Methods (Practice 7.1)

1. **Linear stack** — simple mean. Most common. `stack = mean(CC_i)`

2. **Phase-weighted stack** (Schimmel & Paulssen, 1997):
   - Compute instantaneous phase of each CC via Hilbert transform
   - Weight by coherence of the phase across windows
   - Down-weights incoherent noise more aggressively than linear stack
   - Equation: `pws = linear_stack * |mean(exp(i*phi_k))|^v`

3. **Robust stack** (Pavlis & Vernon, 2010):
   - Iteratively down-weight outliers
   - Median-based: resistant to a few bad windows
   - More computationally expensive but handles transient contamination

### Convergence (Practice 7.2)

Show SNR vs. number of stacked windows. Students see that SNR grows — but slowly. This motivates why real studies use months to years of data.

---

## Week 8 — Module 8: Optimal Processing (Optional/Advanced)
**Cells: 30 | Time: ~1 hour (lecture only, or skip)**

### When to Cover This

This module (Fichtner et al., 2020) is more advanced. Options:
- **Skip entirely** for an introductory course
- **Cover as a lecture** without the practice for an intermediate course
- **Full coverage** for a graduate seminar

### The Core Insight

Standard nonlinear processing (one-bit, whitening) breaks linear wave physics. The transfer coefficient framework quantifies this:

$$T_{ik}^{(n)} = f^{(n)} \cdot g_{ik} + e_{ik}^{(n)}$$

The residual $e$ creates an "unphysical wavefield" — different station pairs see different effective source distributions. This means:
- One-bit: ~10 dB amplitude errors, ~3% traveltime bias
- Whitening: smaller errors, mostly amplitude

**The punchline**: For tomography, these biases may be small enough to ignore. For full waveform inversion or high-precision monitoring, they matter. The "optimal processing" (OPUS) framework eliminates the unphysical component.

**Message for students**: Don't blindly equate CC = Green's function. Understand what your processing does.

---

## Week 9 — Module 9: Surface Wave Dispersion (FTAN)
**Cells: 31–32 | Time: ~3 hours**

### Why Synthetic Data?

**Be transparent with students:** Our real-data CC from Module 6 doesn't show a clear surface wave signal because 1 day of data is insufficient. Rather than pretend, we use synthetic data with known dispersion to teach the method properly. This is pedagogically better — students can verify their FTAN picks against the true model.

### The FTAN Method (Bensen et al., 2007, Phase 3)

1. Apply a narrow **Gaussian bandpass filter** centered at frequency f_0:
   $$G(\omega, \omega_0) = \exp\left[-\alpha\left(\frac{\omega - \omega_0}{\omega_0}\right)^2\right]$$
   
   The width parameter alpha controls the trade-off: narrow filter = good period resolution but poor spatial resolution.

2. Compute the **envelope** via Hilbert transform

3. Map the envelope to **velocity space**: v = distance / t_arrival

4. **Pick the maximum** at each period → group velocity dispersion curve

### The Synthetic CC

We simulate a 15 km station pair (realistic for a 20 km radius network) with dispersion:
- v_g = 0.8 + 0.55 * T (km/s)
- 0.5 s period → 1.1 km/s (shallow sediments)
- 5.0 s period → 3.6 km/s (upper crust)

The spectrum is built by **integrating the group delay** at each frequency — this is the correct way to construct a dispersive waveform.

### Depth Sensitivity

Surface waves at different periods sense different depths. Rule of thumb: peak sensitivity at depth ~ T * v / 3. For a 20 km radius network (periods 0.5-5 s), you image the upper ~5 km.

### Data Selection Criteria (Bensen et al., 2007)

Drill these into students — they will need them for their projects:
- **3-wavelength rule**: distance > 3 * lambda (avoid near-field effects)
- **SNR > 10** (can relax to 7 for sparse networks)
- **Temporal repeatability**: seasonal subsets should agree within 100 m/s

---

## Week 10 — Module 10: Tomographic Inversion
**Cells: 33–34 | Time: ~3 hours**

### From Dispersion Curves to Velocity Maps

**Two-step inversion:**

1. **2D tomography at each period**: Travel times → velocity map
   - Forward problem: ray-theoretic path integrals
   - Inverse problem: penalized least squares (Bensen et al., 2008)
   - Regularization: damping (identity matrix) + smoothing (Laplacian)
   
   Write the objective function on the board and explain each term.

2. **Depth inversion at each grid point**: Dispersion curve → 1D Vs(z) profile
   - Use sensitivity kernels: ∂c/∂Vs
   - Or transdimensional Bayesian (Ryberg et al., 2022; Cabrera-Pérez et al., 2023) — the data choose the number of layers

### Practice 10.1: Synthetic 2D Tomography

Students implement a simplified tomography with:
- Random station locations
- Synthetic velocity model with a slow basin and fast intrusion
- Straight-ray forward operator
- Gauss-Newton inversion with regularization

**Discussion points:**
- Where is resolution good/poor? → depends on ray coverage
- What happens if you change damping? → trade-off between fit and smoothness
- How many data per model parameter is enough? → ideally > 2

### Nonlinear Methods

Mention MANgOSTA (Cabrera-Pérez et al., 2021): for complex settings (islands, volcanoes), linear inversion fails because ray paths depend on the velocity model. The multiscale approach iterates, recomputing rays at each step.

---

## Week 11 — Module 11: Seismic Interferometry
**Cells: 35–36 | Time: ~2 hours**

### Broader Context

Seismic interferometry is the umbrella — ambient noise is just one application. Cover:

1. **Virtual sources**: The CC between stations A and B creates a "virtual seismogram" as if there were a source at A and a receiver at B (or vice versa). No real source needed.

2. **Beyond ambient noise**:
   - Earthquake coda (Campillo & Paul, 2003)
   - Active sources (drill-bit noise, traffic)
   - Autocorrelation (single-station body wave reflections — Claerbout, 1968)
   - C3 method: correlations of correlations

3. **Historical motivation**: Claerbout (1968) proposed "daylight imaging" — turn transmitted noise into reflection seismograms. It took 35 years for the seismology community to realize this works.

---

## Week 12 — Module 12: Velocity Change Monitoring (dv/v)
**Cells: 37–39 | Time: ~3 hours**

### The Stretching Technique

1. **Principle**: If the medium velocity changes by dv/v, all arrival times in the coda shift by dt/t = -dv/v. So we stretch/compress a reference CC and find the stretching factor that maximizes the correlation coefficient.

   $$\text{CC}(\epsilon) = \frac{\int u_{\text{ref}}(t(1-\epsilon)) \cdot u_{\text{cur}}(t) \, dt}{\sqrt{\int u_{\text{ref}}^2 \, dt \int u_{\text{cur}}^2 \, dt}}$$
   
   Maximum at $\epsilon_0 = dv/v$.

2. **Precision**: Can detect velocity changes of 10^{-5} (0.001%). This is extraordinary — no other method achieves this precision for in-situ velocity monitoring.

3. **Applications** (Brenguier et al., 2016):
   - **Volcanoes**: dv/v decreases before eruptions (pressurization opens cracks)
   - **Fault zones**: dv/v drops after large earthquakes (damage), then recovers (healing)
   - **Seasonal cycles**: temperature and groundwater effects (~0.1% annual variation)
   - **Geothermal reservoirs**: fluid injection monitoring

### Practice 12.1-12.2

- Students implement the stretching technique
- Then simulate 365 days of dv/v with seasonal variation + a volcanic "event"
- Practice detecting the event above the seasonal background

### Key Discussion

The tension: short stacks (1 day) give temporal resolution but are noisy. Long stacks (1 month) are precise but blur fast transients. For volcanic monitoring, 5-10 day stacks are typical.

---

## Week 13 — Module 13: HVSR Method
**Cells: 40–41 | Time: ~2.5 hours**

### Different from Everything Else

HVSR is a **single-station** method — no cross-correlation between stations. It uses the spectral ratio of horizontal to vertical components of ambient noise at one location.

1. **What it measures**: The fundamental site resonance frequency f_0 = Vs / (4H)
   - H = thickness of soft layer
   - Vs = shear wave velocity of soft layer
   - f_0 = frequency where horizontal motion is amplified relative to vertical

2. **When to use it**: Quick site characterization, hazard assessment, mapping bedrock depth across a large area with a single portable sensor.

3. **Lognormal statistics** (Cox et al., 2020): f_0 follows a lognormal distribution, NOT normal. Report lognormal median and sigma_ln, not mean ± std. The window-rejection algorithm operates in log-space.

4. **SESAME guidelines**: Minimum recording duration T > 200/f_0, window length > 10/f_0, peak amplitude > 2 for a clear resonance.

### Practice

- Practice 13.1: Compute HVSR from 3-component synthetic data with known resonance
- Students see how the peak frequency maps to layer thickness

---

## Week 14 — Module 14-15: NoisePy Workflow & Case Studies
**Cells: 42–45 | Time: ~3 hours**

### NoisePy Architecture (Jiang & Denolle, 2020)

Walk through the three scripts:
- **S0**: Data download and preparation → HDF5/ASDF format
- **S1**: Cross-correlation for all station pairs (FFT + conj multiply)
- **S2**: Stacking (linear, PWS, robust)

Show how the `ConfigParameters` dataclass maps to everything they've learned:
- `cc_len`, `step` → windowing (Module 6)
- `time_norm`, `freq_norm` → preprocessing (Module 5)
- `maxlag` → lag window
- `stack_method` → stacking (Module 7)

### Case Studies (Module 15)

Cover 3-4 examples from the literature:

1. **Continental scale** — Bensen et al. (2008): US tomography, 203 stations, 24 months
2. **Dense array** — Chmiel et al. (2019): Groningen, 400 stations, 350 m spacing
3. **Volcanic island** — Cabrera-Pérez et al. (2023): La Palma, pre-eruptive velocity change
4. **Mineral exploration** — Ryberg et al. (2022): 400 nodes, detecting tin-tungsten deposits

**Key message**: The same method works at all scales — the physics doesn't change, only the station density and period range.

---

## Week 15 — Module 16: Dense Arrays & Exercises
**Cells: 46–48 | Time: ~2 hours + exercise assignment**

### LARGE-N Arrays

- Low-cost nodal seismometers have enabled hundreds to thousands of sensors
- Short deployments (days to weeks) sufficient at local scale
- All N(N-1)/2 station pairs are computed
- Higher-mode surface waves extractable from dense arrays

### Semester Exercises

Assign 2-3 of the following based on student level:

1. **Effect of preprocessing** (Week 3-4 due): Vary normalization parameters, compare CCs
2. **Convergence analysis** (Week 5-6 due): Stack increasing durations, measure SNR
3. **Mini-tomography** (Week 7-9 due): 10-20 stations, all pairs, 2D velocity map
4. **dv/v monitoring** (Week 10-12 due): 6-12 months, stretching technique, seasonal correction
5. **HVSR survey** (Week 13-14 due): Field deployment, 3-5 sites, lognormal analysis

---

## General Teaching Tips

### Things That Confuse Students

1. **Sign conventions**: `np.correlate(A, B)` vs. `conj(FA)*FB` give different sign conventions for the lag axis. Be consistent and explain your convention clearly.

2. **Causal vs. acausal**: The CC has two sides. Positive lag = wave traveling A→B, negative = B→A. Asymmetry tells you about noise source direction, not about the Earth.

3. **Why doesn't my CC show a signal?** Almost certainly insufficient stacking. 1 day is not enough for most station pairs. This is not a bug — it's the physics.

4. **Group vs. phase velocity**: Group velocity is what the FTAN picks (envelope arrival). Phase velocity is faster and requires inter-station phase measurement. Both contain information about Earth structure but sample depth differently.

5. **The Green's function is approximate**: Students may think CC = exact Green's function. Emphasize it's an approximation that depends on source distribution, stacking duration, and processing choices.

### Board Work vs. Notebook

- **Derive on the board**: CC definition, CC theorem, Green's function relation, Bessel function, RAM equation, stretching technique formula
- **Run in the notebook**: All visualizations, parameter exploration, "what-if" experiments

### Assessment Ideas

- **Quiz**: Theory questions (what is the secondary microseism? what does temporal normalization do? why stack?)
- **Practical exam**: Give students a dataset and ask them to produce a dispersion curve
- **Term paper**: Pick a published ANT study, reproduce (simplified) with the notebook tools, discuss limitations
- **Group project**: Small network deployment (if possible), process with the notebook pipeline

---

## Quick Reference: Key Equations

| # | Equation | Reference |
|---|----------|-----------|
| 1 | $C_{AB}(f) = \hat{U}_A^*(f) \cdot \hat{U}_B(f)$ | CC Theorem |
| 2 | $\langle u_1 u_2^* \rangle \propto \text{Im}[G]$ | Campillo & Roux (2015), Eq. 7 |
| 3 | $C(\omega) = J_0(k\|\vec{r}_1 - \vec{r}_2\|)$ | Aki (1957) |
| 4 | $w_n = \frac{1}{2N+1}\sum|d_j|$ | Bensen et al. (2007), Eq. 1 |
| 5 | $G(\omega) = \exp[-\alpha((\omega-\omega_0)/\omega_0)^2]$ | FTAN filter, Bensen et al. (2007), Eq. 6 |
| 6 | $dv/v = -dt/t = \epsilon_0$ | Stretching technique |
| 7 | $f_0 = V_s / (4H)$ | Site resonance |
| 8 | $\text{SNR} \propto T^{1/n}$, $n \approx 2.5$ | Convergence, Bensen et al. (2007) |

---

## Recommended Reading Order for Students

**Essential (assign these):**
1. Bensen et al. (2007) — The processing bible. Read Sections 1-3.
2. Campillo & Roux (2015) — Theory chapter. Read Sections 1-2 for equations.
3. Jiang & Denolle (2020) — NoisePy paper. Read for workflow understanding.

**Recommended (for motivated students):**
4. Bensen et al. (2008) — US tomography, shows the full workflow applied
5. Brenguier et al. (2016) — Volcano monitoring with dv/v
6. Cox et al. (2020) — HVSR statistics
7. Fichtner et al. (2020) — Optimal processing (advanced)

**Background reading:**
8. Shearer (2019), Chapter 11 — Earth noise in a seismology textbook
9. Nakata, Gualtieri & Fichtner (2019) — The definitive ambient noise book
