# Electrocatalyst Screening - CH4 to CH3OH Conversion

High-throughput computational screening for methane-to-methanol electrocatalyst design using Materials Project data and pymatgen.

---

## 📚 Documentation Overview

Start here based on your needs:

### **🚀 Just Want to Run It?**
→ **[QUICK_START.md](QUICK_START.md)** - 3-step guide to get running in 5 minutes

### **🔑 Need API Key Help?**
→ **[API_KEY_SETUP.md](API_KEY_SETUP.md)** - Complete Materials Project API key guide

### **🛠️ Want to Understand Parameters?**
→ **[REQUIRED_USER_INPUTS.md](REQUIRED_USER_INPUTS.md)** - Detailed parameter explanations

### **🐛 Debugging & Fixes**
- **[BUG_FIX_SUMMARY.md](BUG_FIX_SUMMARY.md)** - AttributeError fix details
- **[DEBUG_SUMMARY.md](DEBUG_SUMMARY.md)** - IDE warnings and fixes
- **[IDE_WARNINGS_FIXED.md](IDE_WARNINGS_FIXED.md)** - All 9 PyCharm warnings resolved
- **[FIXES_APPLIED.md](FIXES_APPLIED.md)** - Initial Jupyter notebook conversion fixes

---

## ⚡ Quick Start (TL;DR)

### 1. Install Dependencies
```bash
pip install -U mp-api pymatgen numpy scipy pandas matplotlib
```

### 2. Get API Key
Visit https://materialsproject.org/api → Copy your free API key

### 3. Set API Key
```bash
# Windows PowerShell
$env:MP_API_KEY = "your_key_here"

# Mac/Linux Bash
export MP_API_KEY="your_key_here"
```

### 4. Run
```bash
python main.py
```

That's it! 🎉

---

## 📖 What This Script Does

The script performs a complete computational screening workflow:

### **Step 1: Pourbaix Phase Stability**
- Queries Materials Project Pourbaix diagram data
- Identifies thermodynamically stable metal oxide at operating conditions
- Visualizes phase diagram with operating point

### **Step 2: Electronic Structure Analysis**
- Downloads DFT-calculated density of states (DOS)
- Aligns DOS to operational Fermi level (vacuum scale)
- Extracts d-band center and valence band maximum

### **Step 3: Quantum Descriptors**
- Calculates orbital overlap integrals for CH4 vs CH3OH
- Computes selectivity indices
- Estimates reaction energetics

### **Outputs:**
- `pourbaix_diagram.png` - Phase stability visualization
- `dos_windows_*.png` - Density of states analysis
- `screening_results_*.csv` - Numerical data table

---

## 🎯 Example Usage

### Run for Platinum at 1.5V vs RHE, pH 10
```bash
python main.py
# Uses defaults: Pt, 1.5V, pH=10
```

### Change Metal & Conditions
Edit lines 29-33:
```python
metal_symbol = "Co"      # Change to cobalt
E_RHE = 1.2              # Lower potential
pH = 14.0                # Higher pH (alkaline)
```

### Screen Multiple Metals
Uncomment lines at end of script (~365) and add your metals:
```python
metals_to_screen = ["Pt", "Ir", "Co", "Ni", "Cu"]
```

---

## 📋 Requirements

- **Python 3.8+**
- **pymatgen** - Materials science toolkit
- **mp-api** - Materials Project API client
- **numpy, scipy** - Scientific computing
- **pandas** - Data handling
- **matplotlib** - Plotting
- **Internet connection** - For MP API queries
- **Materials Project API key** - Free from materialsproject.org

---

## ⚠️ Important Notes

### API Key Security
- **Never commit** your API key to GitHub
- **Use environment variables** in production
- **Keep private** - don't share publicly
- **Regenerate** if exposed

### Computational Considerations
- Requires internet (MP API queries)
- ~30-60 seconds per metal (API + DOS download)
- DOS data ~100KB per material

### Physics Approximations
- Rigid DOS alignment (not work function calc)
- Vacuum scale energy reference (Trasatti convention)
- DFT functionals from Materials Project

---

## 🔍 Troubleshooting

### "API KEY AUTHENTICATION FAILED"
→ See [API_KEY_SETUP.md](API_KEY_SETUP.md)

### "POURBAIX ANALYSIS FAILED"
→ Check internet connection, verify API key validity

### "DOS data unavailable"
→ Script will use fallback descriptors (d-band center estimate)

### IDE warnings/errors
→ See [IDE_WARNINGS_FIXED.md](IDE_WARNINGS_FIXED.md)

---

## 📊 Understanding Outputs

### Pourbaix Diagram
Shows thermodynamic stability of different phases (solids, ions) across pH and potential ranges. Red dot = your operating point.

### DOS Windows
Left: CH4 activation window (HOMO = -12.80 eV)  
Right: CH3OH over-oxidation window (HOMO = -10.00 eV)  
Shaded regions = integration windows for hybridization

### CSV Results
Key metrics:
- **E_F_active** - Operational Fermi level (vacuum scale)
- **d_band_center** - Descriptor for catalytic activity
- **I_tot_CH4 / I_tot_CH3OH** - Hybridization integrals
- **S_kinetic** - Selectivity index (>1 = selective for CH4)

---

## 🚀 Next Steps

1. **Run the script** with your API key
2. **Examine the outputs** in PNG and CSV files
3. **Modify parameters** to explore different conditions
4. **Screen multiple metals** by editing the loop

---

## 📝 Citation

If using this script in research, please cite:
- Materials Project: https://materialsproject.org
- pymatgen: https://pymatgen.org
- Pourbaix analysis methodology in your work

---

## 📧 Support

For issues:
1. Check the relevant documentation file above
2. Verify API key is valid (test at https://materialsproject.org/api)
3. Ensure all dependencies installed: `pip install -U mp-api pymatgen`

---

## ✨ Features

✅ Complete automation - single script, one command  
✅ Robust error handling - clear error messages  
✅ Production-ready - type hints, validation  
✅ Well documented - 6 guides + inline comments  
✅ Secure - environment variable support  
✅ Cross-platform - Windows/Mac/Linux  

---

**Last Updated:** 2026-07-20  
**Status:** Production Ready ✅
