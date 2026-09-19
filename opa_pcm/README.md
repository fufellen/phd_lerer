# PCM phase-shifter errors in a one-dimensional optical phased array

This subproject propagates phase-shifter uncertainty from the EIM/planar PCM
LR-DLSPP comparison to array-level metrics.

## Scientific scope

The array is a uniform one-dimensional array with isotropic elements,
wavelength 1.55 um and pitch lambda/2. Arrays with 16, 32 and 64 channels are
evaluated. The model separates:

1. systematic phase-scale error caused by an error in predicted L_pi;
2. phase-dependent insertion loss of a segmented 0..2pi PCM phase shifter;
3. independent random phase error.

The phase-shifter inputs are the EIM and strict planar metrics for GSST, Sb2S3
and Sb2Se3. The planar result is not presented as full-vector truth. The gap
between the EIM prediction and planar response is propagated as model
uncertainty.

For an EIM-designed device evaluated with the planar model:

```text
phase scale = L_pi,EIM / L_pi,planar
```

The 0..2pi device has total length 2 L_pi. A channel with commanded wrapped
phase q uses crystalline fraction f=q/(2pi):

```text
loss_dB = 2 * length_scale * ((1-f) IL_a + f IL_c)
```

## Run

```powershell
python opa_pcm/scripts/simulate_pcm_opa.py
```

Only NumPy and Matplotlib are required.

## Outputs

- `results/deterministic_metrics.csv`: N=16/32/64 at a 30 degree target;
- `results/steering_sweep_n32.csv`: target angles from -60 to 60 degrees;
- `results/random_phase_tolerance.csv`: Monte Carlo tolerance for 0..30 degree
  RMS phase errors;
- `results/patterns_n32_theta30.png`;
- `results/metrics_n32_theta30.png`;
- `results/random_phase_tolerance.png`;
- `results/steering_sweep_n32.png`.

## Validation

The script checks:

- the ideal broadside peak;
- the analytical half-power beamwidth approximation;
- the -13.15 dB sidelobe level of a uniformly excited linear array;
- unity ideal peak and coherent efficiency;
- equivalence of wrapped and unwrapped phase commands;
- invariance to a common phase;
- agreement of Monte Carlo target efficiency with the analytical expectation
  for independent Gaussian phase errors.

## Full-vector characteristic and PCM levels (`fem_and_levels.py`)

```powershell
python opa_pcm/scripts/fem_and_levels.py
```

The same array model with the full-vector (FEM, COMSOL mode analysis) characteristic of the
shifter: L_pi and attenuation of both states from the APEDE-2026 check (110 nm of PCM under the
gold in the FEM model). The EIM lengths are taken at the same geometry, with the second EIM stage
in TM polarization (as in the article and in VIBR_T) and in TE polarization (the polarization of
the vertical field of the mode). A device sized by EIM is evaluated with the FEM response.

Second part: a FEM-calibrated device whose command phase is rounded to M PCM levels
(M = 2...64), swept over target angles -60...60 deg for N = 32.

Checks: a FEM-calibrated unquantized beam points exactly at every tested angle; for lossless
elements the mean target efficiency over angles follows 1/N + (1 - 1/N) sinc^2(pi/M) within 0.02.

Outputs: `results/fem_metrics.csv`, `results/fem_steering_sweep_n32.csv`,
`results/pcm_levels_n32.csv`, `results/fem_levels_n32.png`.

Result (N = 32): with TE on the second stage the EIM-sized device steers within 0.05 deg of the
target over -60...60 deg and keeps the target efficiency within 0.5 % of the calibrated device;
with TM, Sb2S3 (phase scale 1.76) forms false global maxima (target efficiency down to 0.06) and
Sb2Se3 (scale 1.38) raises the worst sidelobe to -3.6 dB. Phase quantization of the low-loss
shifters costs 0.87 / 0.21 / 0.05 dB on average at M = 4 / 8 / 16 and limits the worst sidelobe
to -6.0 / -9.4 / -11.2 dB against -13.2 dB without quantization. GSST is loss-limited (target
efficiency 0.19-0.99 depending on the angle) and is left out of the quantization panels.

## Limits

The result is an array-factor study. It does not include an emitter element
factor, mutual coupling, splitter imbalance, thermal crosstalk, wavelength
dispersion or full-vector PCM phase-shifter fields. The phase-dependent-loss
law assumes a segmented device with linearly mixed amorphous and crystalline
lengths.
