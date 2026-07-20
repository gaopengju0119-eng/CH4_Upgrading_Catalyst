# Electrocatalyst Screening — CH₄ → CH₃OH Conversion

High-throughput computational screening for selective methane-to-methanol electrocatalyst
design using Materials Project data and pymatgen. The workflow pins down the stable phase at
the operating condition, aligns its DFT density of states (DOS) to the operational Fermi level,
and computes a spin-aware selectivity index.

---

## ⚡ Quick Start

### 1. Create the environment
```bash
conda env create -f environment.yml
conda activate catalyst-screening
```
Or update an existing env in place:
```bash
conda env update -f environment.yml --prune
```

### 2. Get a Materials Project API key
Visit <https://materialsproject.org/api> and copy your free API key.

### 3. Set your API key
Edit the `MP_API_KEY` variable near the top of `main.py` (line 41):
```python
MP_API_KEY = "your_key_here"
```
> ⚠️ The repository currently contains a hardcoded key. Replace it with your own and avoid
> committing real keys — prefer an environment variable or Colab Secret in shared settings.

### 4. Run
```bash
python main.py
```

---

## 📖 What This Script Does

`main.py` runs a three-step diagnostic for a single target metal:

### Step 1 — Pourbaix Phase Stability
- Fetches Pourbaix entries from the Materials Project (`get_pourbaix_entries`).
- Converts the applied potential from the RHE to the SHE scale via the Nernst relation.
- Finds the thermodynamically stable phase at `(pH, E_SHE)` and resolves its `mp-id`.
- Raises a clear error if the stable domain is ionic or spans multiple solids (no unique bulk DOS).

### Step 2 — Electronic Structure & Magnetic Analysis
- Retrieves the summary document and DFT DOS for the stable `mp-id`.
- Rigidly aligns the DOS to the operational Fermi level on the vacuum scale.
- Detects spin polarization from the actual `Spin.up` / `Spin.down` channels and diagnoses
  magnetism from summary metadata, total moment, and DOS spin-asymmetry.
- Plots the spin-resolved alignment diagnostic.

### Step 3 — Selectivity Index
- Integrates the DOS over the CH₄ activation window and the CH₃OH over-oxidation window.
- Computes the selectivity index `S = I_total_CH4 / I_total_CH3OH` (spin-summed).

### Outputs
- `screening_total_dos_<metal>.csv` — full numerical diagnostic table (e.g. `screening_total_dos_Ni.csv`).
- Two matplotlib figures shown on screen: the alignment diagnostic and a dual-panel DOS window plot.
  (Figures are displayed via `plt.show()`, not saved to disk.)

---

## 🎛️ Configuration

All user parameters live in the **`USER INPUT PARAMETERS`** block (`main.py`, lines 41–53).
Defaults as shipped:

| Parameter       | Default   | Meaning                                   |
| --------------- | --------- | ----------------------------------------- |
| `metal_symbol`  | `"Ni"`    | Target transition metal                   |
| `E_RHE`         | `2.0`     | Applied potential vs. RHE (V)             |
| `pH`            | `13.0`    | Electrolyte pH                            |
| `temp_k`        | `298.15`  | Temperature (K)                           |
| `ion_conc`      | `1e-6`    | Ion concentration (mol/L)                 |
| `E_SHE_abs`     | `-4.44`   | Absolute SHE reference (eV, vacuum scale) |
| `E_HOMO_CH4`    | `-12.80`  | CH₄ HOMO level (eV)                        |
| `E_HOMO_CH3OH`  | `-10.00`  | CH₃OH HOMO level (eV)                      |

### Change metal & conditions
```python
metal_symbol = "Co"   # e.g. cobalt
E_RHE = 1.2           # lower potential
pH = 14.0             # more alkaline
```

### Screen multiple metals
A commented scaffold sits at the end of `main.py` (Section 9, ~lines 336–344):
```python
metals_to_screen = ["Pt", "Ir", "Co", "Ni", "Cu", "Mn"]
```
Wrap Steps 1–3 in a function to run the loop in production.

---

## 📋 Requirements

Managed via `environment.yml`:

- **Python 3.11**
- **numpy** (`<2.1`), **scipy**, **pandas** (`==2.2.2`), **matplotlib**
- **mp-api** — Materials Project API client
- **mpcontribs-client** — **required** by `get_pourbaix_entries` for ion reference data
- **pymatgen** — materials science toolkit
- **pyarrow**, **boltons**, **requests** (`==2.32.4`)
- **Internet connection** — for Materials Project API queries
- **Materials Project API key** — free from materialsproject.org

---

## 🔍 Troubleshooting

### `AttributeError: 'NoneType' object has no attribute 'query_contributions'`
`get_pourbaix_entries` needs the MPContribs client, which only initializes if
`mpcontribs-client` is installed. Install it (already pinned in `environment.yml`):
```bash
pip install mpcontribs-client
```

### Pourbaix / API errors
Check your internet connection and verify the API key at <https://materialsproject.org/api>.

### "stable Pourbaix domain is ionic or contains multiple solids"
The stable phase at your `(pH, E_RHE)` has no unique bulk `mp-id`. Adjust the potential/pH,
or add an explicit phase-selection rule.

### Magnetic material with single-channel DOS
The script refuses to fabricate spin-up/spin-down curves. Use a spin-resolved MP calculation.

---

## 📊 Understanding the CSV

Key columns in `screening_total_dos_<metal>.csv`:
- **Stable Phase / mp-id** — the resolved stable bulk phase.
- **Magnetic / Magnetic Ordering / Total Magnetization / DOS Spin Asymmetry** — magnetic diagnosis.
- **DFT Fermi / Energy Shift / Vacuum Fermi** — DOS alignment on the vacuum scale.
- **Center_elec** — occupied-state electronic center (eV, vacuum).
- **I_total_CH4 / I_total_CH3OH** — hybridization integrals over each window.
- **S_selectivity** — selectivity index (`>1` favors CH₄ activation).

---

## 📝 Citation

If using this script in research, please cite:
- Materials Project — <https://materialsproject.org>
- pymatgen — <https://pymatgen.org>

---

**Last Updated:** 2026-07-20
