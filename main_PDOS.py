# -*- coding: utf-8 -*-
"""High-throughput electrocatalyst screening (Colab-compatible Python script).

Selective methane-to-methanol upgrading via electrochemical oxidation.
Materials Project + pymatgen.
"""

# ============================================================
# 0. Install Dependencies
# ============================================================
# In Google Colab, run this in a separate installation cell:
# !pip install --upgrade pyarrow mp-api mpcontribs-client "pymatgen>=2026.1" boltons "numpy<2.1" -q
# !pip install pandas==2.2.2 "requests>=2.32.5" --no-deps --force-reinstall -q
#
# Note: mp-api pulls in emmet-core>=0.87, which imports StructureGraph/MoleculeGraph
# from pymatgen.core.graphs (added in pymatgen 2026.1). An older pymatgen raises
# "ModuleNotFoundError: No module named 'pymatgen.core.graphs'".

# ============================================================
# 1. Imports
# ============================================================
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mp_api.client import MPRester
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter
from pymatgen.core import Element
from pymatgen.electronic_structure.core import OrbitalType, Spin

warnings.filterwarnings("ignore")


def _trapz(y, x):
    """Universal trapezoidal integration helper (NumPy 1.x and 2.x)."""
    if hasattr(np, "trapezoid"):
        return np.trapezoid(y, x)
    return np.trapz(y, x)


# ============================================================
# 2. USER INPUT PARAMETERS
# ============================================================
MP_API_KEY = "VhG42T1WDDotMqdqcbId4MJGJVblvtgP"  # Set this in Colab Secrets/environment.
metal_symbol = "Ru"       # Target transition metal
E_RHE = 2.0                # Applied potential vs. RHE (V)
pH = 13.0                  # Electrolyte pH
temp_k = 298.15            # Temperature (K)
ion_conc = 1e-6            # Ion concentration (mol/L)

# Fixed reference energy levels (vacuum scale, eV)
E_SHE_abs = -4.44
E_HOMO_CH4 = -12.80
E_HOMO_CH3OH = -10.80
MAG_DOS_REL_TOL = 1e-3     # Relative integrated spin-asymmetry threshold
MAG_MOMENT_TOL = 1e-3      # µB per formula unit

# Output artifacts (plots + CSV) are written to an "Output" folder next to this
# script. Falls back to the current working directory when __file__ is undefined
# (e.g. in a Colab/Jupyter cell).
try:
    _BASE_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _BASE_DIR = os.getcwd()
OUTPUT_DIR = os.path.join(_BASE_DIR, "Output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"[Output] Artifacts will be saved to: {OUTPUT_DIR}")


# ============================================================
# 3. UNIT CONVERSION: E_RHE -> E_SHE
# E_SHE = E_RHE - (kT ln(10)/e) * pH
# Note: Pourbaix diagrams use V vs. SHE.
# ============================================================
nernst_slope = 8.617333262e-5 * temp_k * np.log(10.0)  # V per pH unit
E_SHE = E_RHE - nernst_slope * pH
print(f"[Unit Conversion] E_RHE = {E_RHE:.3f} V -> E_SHE = {E_SHE:.3f} V (pH={pH:g})")

# ============================================================
# 4. Operando Phase Stability — Pourbaix Analysis
# ============================================================
print("\n" + "=" * 60)
print("STEP 1: Pourbaix Phase Stability Diagnosis")
print("=" * 60)

# Initialize MPRester - get_pourbaix_entries requires ion reference data from contributions
mpr = MPRester(MP_API_KEY)
pbx_entries = mpr.get_pourbaix_entries([metal_symbol])

pbx = PourbaixDiagram(
    pbx_entries,
    conc_dict={metal_symbol: ion_conc},
    filter_solids=True,
)
stable_entry = pbx.find_stable_entry(pH=pH, V=E_SHE)
stable_formula = stable_entry.name

# ---- Pourbaix diagram plot with the operating point marked ----
plotter = PourbaixPlotter(pbx)
ax_pbx = plotter.get_pourbaix_plot(
    title=f"Pourbaix diagram: {metal_symbol} (ion conc = {ion_conc:g} M)",
)
ax_pbx.scatter(
    [pH], [E_SHE],
    marker="*", s=350, color="red", edgecolor="black", zorder=5,
    label=f"Operating point (pH={pH:g}, {E_SHE:.2f} V vs. SHE)",
)
ax_pbx.set_ylabel("E (V vs. SHE)")
ax_pbx.legend(loc="best", fontsize=10)
pourbaix_png = os.path.join(OUTPUT_DIR, f"pourbaix_{metal_symbol}.png")
ax_pbx.figure.savefig(pourbaix_png, dpi=300, bbox_inches="tight")
print(f"Saved: {pourbaix_png}")
plt.show()


def extract_solid_mp_id(entry):
    """Return a unique solid mp-id from a PourbaixEntry/MultiEntry, else None."""
    candidates = []
    components = getattr(entry, "entry_list", None) or [entry]
    for component in components:
        inner = getattr(component, "entry", component)
        entry_id = getattr(inner, "entry_id", None) or getattr(component, "entry_id", None)
        phase_type = str(getattr(component, "phase_type", getattr(entry, "phase_type", ""))).lower()
        if entry_id and str(entry_id).startswith("mp-") and phase_type != "ion":
            candidates.append(str(entry_id))
    unique_ids = sorted(set(candidates))
    return unique_ids[0] if len(unique_ids) == 1 else None


stable_entry_id = extract_solid_mp_id(stable_entry)
# Extract base mp-id if it contains calculation type suffix (e.g., "mp-123-GGA+U" -> "mp-123")
if stable_entry_id and "-" in stable_entry_id:
    # Split and take only the first part (base mp-id)
    base_mp_id = stable_entry_id.split("-")[0] + "-" + stable_entry_id.split("-")[1]
    if base_mp_id.startswith("mp-"):
        stable_entry_id = base_mp_id
print(f"  [Stable phase]  : {stable_formula}")
print(f"  [Stable mp-id]  : {stable_entry_id or 'none (ion or multi-solid domain)'}")
if stable_entry_id is None:
    raise RuntimeError(
        "The stable Pourbaix domain is ionic or contains multiple solids, so a unique bulk DOS "
        "cannot be assigned without an explicit phase-selection rule."
    )

# ============================================================
# 5. Electronic Structure & Magnetic Analysis
# ============================================================
print("\n" + "=" * 60)
print("STEP 2: Electronic Structure & Magnetic Analysis")
print("=" * 60)

# ---- 2a. Active Fermi level (vacuum scale) ----
E_F_active = E_SHE_abs - E_SHE

# ---- 2b. DOS retrieval & magnetic metadata ----
mpr = MPRester(MP_API_KEY)
summary_docs = mpr.materials.summary.search(
    material_ids=[stable_entry_id],
    fields=["material_id", "formula_pretty", "is_magnetic", "ordering", "total_magnetization"],
)
if not summary_docs:
    raise RuntimeError(f"No Materials Project summary document for {stable_entry_id}.")
summary = summary_docs[0]
dos_data = mpr.get_dos_by_material_id(stable_entry_id)

if dos_data is None:
    raise RuntimeError(f"Could not retrieve DOS for {stable_entry_id}.")

formula_pretty = getattr(summary, "formula_pretty", stable_formula)
summary_is_magnetic = bool(getattr(summary, "is_magnetic", False))
total_mag_raw = getattr(summary, "total_magnetization", 0.0)
total_mag = float(total_mag_raw) if total_mag_raw is not None else 0.0
ordering = str(getattr(summary, "ordering", "Unknown"))

efermi_dft = float(dos_data.efermi)
energies_raw = np.asarray(dos_data.energies, dtype=float)
energy_shift = E_F_active - efermi_dft
energies_abs = energies_raw + energy_shift

# ---- 2c. Handle spin polarization and diagnose magnetism ----
# Use the actual channel mapping. get_densities() with no spin argument returns
# a spin-summed array and therefore cannot be used to detect spin polarization.
spin_densities = dos_data.densities
has_up = Spin.up in spin_densities
dos_up = np.asarray(spin_densities[Spin.up], dtype=float) if has_up else None
has_down = Spin.down in spin_densities
dos_down = np.asarray(spin_densities[Spin.down], dtype=float) if has_down else None
is_spin_polarized = has_up and has_down

if not has_up:
    raise RuntimeError("The retrieved DOS contains no Spin.up channel.")

if is_spin_polarized:
    spin_difference = _trapz(np.abs(dos_up - dos_down), energies_raw)
    spin_weight = _trapz(np.abs(dos_up) + np.abs(dos_down), energies_raw)
    spin_asymmetry = spin_difference / max(spin_weight, 1e-30)
else:
    spin_asymmetry = 0.0

# Summary metadata is the primary magnetic label; DOS asymmetry and moment provide
# transparent fallbacks for older/incomplete records.
is_magnetic = (
    summary_is_magnetic
    or abs(total_mag) > MAG_MOMENT_TOL
    or spin_asymmetry > MAG_DOS_REL_TOL
)

if is_magnetic and not is_spin_polarized:
    raise RuntimeError(
        f"{stable_entry_id} is reported magnetic, but its MP DOS has only one spin channel. "
        "Refusing to fabricate spin-up/spin-down curves; use a spin-resolved calculation."
    )

if is_spin_polarized:
    total_dos = np.abs(dos_up) + np.abs(dos_down)
else:
    total_dos = np.abs(dos_up)

print(f"  [Material ID]   : {stable_entry_id}")
print(f"  [Formula]       : {formula_pretty}")
print(f"  [MP ordering]   : {ordering}")
print(f"  [Total moment]  : {total_mag:.6f} uB")
print(f"  [Spin channels] : {'up + down' if is_spin_polarized else 'single channel'}")
print(f"  [DOS asymmetry] : {spin_asymmetry:.6e}")
print(f"  [Is magnetic]   : {is_magnetic}")

plt.figure(figsize=(9, 6))
if is_magnetic:
    plt.plot(dos_up, energies_abs, color="tab:blue", lw=1.8, label="Spin up")
    plt.plot(-dos_down, energies_abs, color="tab:red", lw=1.8, label="Spin down")
    plt.fill_betweenx(energies_abs, 0, dos_up, color="tab:blue", alpha=0.15)
    plt.fill_betweenx(energies_abs, 0, -dos_down, color="tab:red", alpha=0.15)
else:
    plt.plot(total_dos, energies_abs, color="black", lw=1.8, label="Total DOS")
    plt.fill_betweenx(energies_abs, 0, total_dos, color="gray", alpha=0.2)

plt.axvline(0, color="black", lw=0.8, alpha=0.5)
plt.axhline(E_F_active, color="green", ls="--", lw=2, label=f"Fermi level ({E_F_active:.2f} eV)")
plt.ylabel("Energy (eV vs. vacuum)")
plt.xlabel("DOS (states/eV)")
plt.title(f"Spin-resolved alignment diagnostic: {formula_pretty} ({stable_entry_id})")
plt.legend(loc="best")
plt.grid(alpha=0.3)
plt.tight_layout()
alignment_png = os.path.join(OUTPUT_DIR, f"dos_alignment_{metal_symbol}.png")
plt.savefig(alignment_png, dpi=300, bbox_inches="tight")
print(f"Saved: {alignment_png}")
plt.show()

# ---- 2d. Center calculation (spin-summed) ----
occupied_mask = energies_abs <= E_F_active
occ_E = energies_abs[occupied_mask]
occ_dos = total_dos[occupied_mask]
occ_norm = _trapz(occ_dos, occ_E) if occ_E.size >= 2 else 0.0
electronic_center = _trapz(occ_E * occ_dos, occ_E) / occ_norm if occ_norm > 1e-12 else np.nan

print("  [Alignment Summary]")
print(f"      - Shift applied: {energy_shift:.4f} eV")
print(f"      - Vacuum Fermi (eV): {E_F_active:.4f} eV")

# ============================================================
# 5b. Orbital-Projected DOS (metal d-band, O 2p-band)
# ============================================================
# get_dos_by_material_id returns a CompleteDos, which carries the element- and
# orbital-projected DOS (the same data the Materials Project site plots). Not every
# MP entry stores projections, and O only exists for oxide/hydroxide phases, so every
# projection lookup is guarded and the script degrades to total-DOS-only cleanly.
print("\n" + "=" * 60)
print("STEP 2b: Orbital-Projected DOS (PDOS)")
print("=" * 60)


def _sum_spins(dos_obj):
    """Spin-summed density array for a pymatgen Dos, on its native energy grid."""
    return sum(np.asarray(d, dtype=float) for d in dos_obj.densities.values())


def occupied_band_center(density):
    """DOS-weighted centroid of the occupied states (vacuum scale, eV)."""
    density = np.asarray(density, dtype=float)
    mask = energies_abs <= E_F_active
    e_occ = energies_abs[mask]
    g_occ = density[mask]
    norm = _trapz(g_occ, e_occ) if e_occ.size >= 2 else 0.0
    return _trapz(e_occ * g_occ, e_occ) / norm if norm > 1e-12 else np.nan


metal_d_center = np.nan
o_p_center = np.nan
projected_curves = []  # list of (label, density, color)

if getattr(dos_data, "pdos", None):
    # Metal d-band (always attempt the target metal).
    try:
        metal_d = dos_data.get_element_spd_dos(Element(metal_symbol))[OrbitalType.d]
        metal_d_density = _sum_spins(metal_d)
        metal_d_center = occupied_band_center(metal_d_density)
        projected_curves.append((f"{metal_symbol} d", metal_d_density, "tab:blue"))
    except Exception as exc:  # noqa: BLE001 - projection may be absent for this entry
        print(f"  [PDOS] no {metal_symbol}-d projection: {exc}")

    # Oxygen 2p-band (only present for oxide/hydroxide phases).
    try:
        o_p = dos_data.get_element_spd_dos(Element("O"))[OrbitalType.p]
        o_p_density = _sum_spins(o_p)
        o_p_center = occupied_band_center(o_p_density)
        projected_curves.append(("O p", o_p_density, "tab:red"))
    except Exception:  # noqa: BLE001 - no oxygen in this phase
        print("  [PDOS] no oxygen in this phase; skipping O-p")
else:
    print("  [PDOS] No projected DOS stored for this material; total DOS only.")

print(f"  [Metal d-band center] : {metal_d_center:.4f} eV (vacuum)")
print(f"  [O 2p-band center]    : {o_p_center:.4f} eV (vacuum)")

if projected_curves:
    plt.figure(figsize=(9, 6))
    for label, density, color in projected_curves:
        plt.plot(density, energies_abs, color=color, lw=1.8, label=label)
        plt.fill_betweenx(energies_abs, 0, density, color=color, alpha=0.12)
    plt.axvline(0, color="black", lw=0.8, alpha=0.5)
    plt.axhline(E_F_active, color="green", ls="--", lw=2, label=f"Fermi level ({E_F_active:.2f} eV)")
    plt.axhline(E_HOMO_CH4, color="purple", ls=":", lw=1.2, label=f"CH4 HOMO ({E_HOMO_CH4:.2f} eV)")
    plt.axhline(E_HOMO_CH3OH, color="brown", ls=":", lw=1.2, label=f"CH3OH HOMO ({E_HOMO_CH3OH:.2f} eV)")
    plt.ylabel("Energy (eV vs. vacuum)")
    plt.xlabel("Projected DOS (states/eV)")
    plt.title(f"Orbital-projected DOS: {formula_pretty} ({stable_entry_id})")
    plt.legend(loc="best", fontsize=9)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    pdos_png = os.path.join(OUTPUT_DIR, f"pdos_{metal_symbol}.png")
    plt.savefig(pdos_png, dpi=300, bbox_inches="tight")
    print(f"Saved: {pdos_png}")
    plt.show()

# ============================================================
# 6. Quantum Descriptors & Selectivity Metrics (Spin-Summed)
# ============================================================
print("\n" + "=" * 60)
print("STEP 3: Selectivity Index (Spin-Polarized Aware)")
print("=" * 60)


def integrate_spin_dos(energies, up, down, E_low, E_high):
    """Integrate each available spin channel, then return up, down, and total."""
    energies = np.asarray(energies, dtype=float)
    mask = (energies >= E_low) & (energies <= E_high)
    if mask.sum() < 2:
        return 0.0, 0.0, 0.0
    i_up = float(_trapz(np.abs(up[mask]), energies[mask]))
    i_down = float(_trapz(np.abs(down[mask]), energies[mask])) if down is not None else 0.0
    return i_up, i_down, i_up + i_down


I_up_CH4, I_dn_CH4, I_tot_CH4 = integrate_spin_dos(
    energies_abs, dos_up, dos_down, E_HOMO_CH4, E_F_active
)
I_up_CH3OH, I_dn_CH3OH, I_tot_CH3OH = integrate_spin_dos(
    energies_abs, dos_up, dos_down, E_HOMO_CH3OH, E_F_active
)
S_selectivity = I_tot_CH4 / (I_tot_CH3OH + 1e-20)

print("[1] Spin-resolved integration results:")
print(f"    Is magnetic: {is_magnetic}")
print(f"    I_CH4   up/down/total: {I_up_CH4:.6f} / {I_dn_CH4:.6f} / {I_tot_CH4:.6f}")
print(f"    I_CH3OH up/down/total: {I_up_CH3OH:.6f} / {I_dn_CH3OH:.6f} / {I_tot_CH3OH:.6f}")
print(f"    Selectivity index: {S_selectivity:.6f}")

# ============================================================
# 7. OUTPUT: Simplified Diagnostic Report
# ============================================================
print("\n" + "=" * 60)
print("DIAGNOSTIC REPORT SUMMARY (TOTAL DOS)")
print("=" * 60)

results_dict = {
    "Metal Symbol": metal_symbol,
    "E_RHE (V)": E_RHE,
    "pH": pH,
    "Stable Phase": stable_formula,
    "mp-id": stable_entry_id,
    "Magnetic": is_magnetic,
    "Magnetic Ordering": ordering,
    "Total Magnetization (uB)": round(total_mag, 6),
    "DOS Spin Asymmetry": round(spin_asymmetry, 8),
    "DFT Fermi (eV)": round(efermi_dft, 4),
    "Energy Shift (eV)": round(energy_shift, 4),
    "Vacuum Fermi (eV)": round(E_F_active, 4),
    "Center_elec (eV)": round(electronic_center, 4),
    "Metal_d_center (eV)": round(metal_d_center, 4),
    "O_p_center (eV)": round(o_p_center, 4),
    "I_up_CH4": round(I_up_CH4, 6),
    "I_down_CH4": round(I_dn_CH4, 6),
    "I_total_CH4": round(I_tot_CH4, 6),
    "I_up_CH3OH": round(I_up_CH3OH, 6),
    "I_down_CH3OH": round(I_dn_CH3OH, 6),
    "I_total_CH3OH": round(I_tot_CH3OH, 6),
    "S_selectivity": round(S_selectivity, 6),
}

df_results = pd.DataFrame([results_dict]).T
df_results.columns = ["Value"]
print(df_results.to_string())
output_csv = os.path.join(OUTPUT_DIR, f"screening_total_dos_{metal_symbol}.csv")
df_results.to_csv(output_csv)
print(f"\nSaved: {output_csv}")

# ============================================================
# 8. VISUALIZATION: Dual-Scale DOS Analysis (Spin-Resolved)
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 8), sharey=True)
y_min_vac = min(E_HOMO_CH4, E_HOMO_CH3OH) - 3.0
y_max_vac = E_F_active + 5.0

for ax, window_lo, window_label, color_shade in zip(
    axes,
    [E_HOMO_CH4, E_HOMO_CH3OH],
    ["CH4 Activation Window", "CH3OH Over-oxidation Window"],
    ["#AED6F1", "#ABEBC6"],
):
    if is_magnetic:
        ax.plot(dos_up, energies_abs, color="tab:blue", lw=1.5, label="Spin up")
        ax.plot(-dos_down, energies_abs, color="tab:red", lw=1.5, label="Spin down")
    else:
        ax.plot(total_dos, energies_abs, color="black", lw=1.5, label="Total DOS")

    ax.axhspan(window_lo, E_F_active, color=color_shade, alpha=0.35, label=window_label)
    ax.axhline(E_F_active, color="green", ls="--", lw=1.5, label="Fermi level")
    ax.axhline(window_lo, color="purple", ls=":", lw=1.2)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_ylim(y_min_vac, y_max_vac)
    ax.set_xlabel("DOS (states/eV; spin down mirrored)")
    ax.set_title(window_label)
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=9)

axes[0].set_ylabel("Energy (eV vs. vacuum)")
fig.suptitle(
    f"Selectivity Screening: {formula_pretty} ({stable_entry_id})\n"
    f"{'Spin-resolved magnetic DOS' if is_magnetic else 'Nonmagnetic DOS'}",
    fontsize=14,
)
plt.tight_layout()
selectivity_png = os.path.join(OUTPUT_DIR, f"selectivity_dos_{metal_symbol}.png")
fig.savefig(selectivity_png, dpi=300, bbox_inches="tight")
print(f"Saved: {selectivity_png}")
plt.show()