# W51 au scripts (A2: adopted distance 5.41 kpc)

Python-script versions of the three notebooks that produce every distance-dependent (au) value in
the paper. The code is copied from the notebooks in `W51-Masers-and-Starforming-Regions/notebooks/`.
Each `# %% cell N` marker is one notebook cell.

Apart from comments, only these changes were made to the notebook code:

- **Figures save, never open.** Every script starts with `matplotlib.use('Agg')`.
  - `SED_modeling.py`: in `fit_2dgaussian`, `plt.show()` became `plt.savefig("gaussian2D_fit_<source>_<freq>GHz.png")` + `plt.close()`.
  - `new_publication_plots.py`: the final `plt.show()` and `pl.show()` were removed.
- **Results print.** Each cell-ending expression the notebook displayed is wrapped in `print()`, with a
  label for the au results. `SED_modeling.py` also prints `ca3.physprops()` and `d2.physprops()`
  (size in au, density, mass, Nlyc, EM).
- **Figure 1 saves the real plot.** `fig.savefig('Overview_W51IRS2_6', dpi=800)` became
  `ax.figure.savefig(...)`. `fig` is an empty figure; `w51north_plot()` draws on a new one.

| Script | Source notebook | Produces (main.tex) |
|---|---|---|
| `ca3_near_masers.py` | `notebooks/ca3_near_masers.ipynb` | Maser separations in arcsec/au and W51 North P.A. (lines 161, 166, 169) |
| `SED_modeling.py` | `notebooks/SED_modeling.ipynb` | CA3 radius 48 au, density, Nlyc (lines 263, 279, 288); SED panels for Figure 4 |
| `new_publication_plots.py` | `notebooks/new_publication_plots.ipynb` | Figure 1 and its 250 au / 1400 au / 0.02 pc scale bars |

## How to run

Put the scripts anywhere and run them from a terminal in the same Python environment as the notebooks:

```bash
python ca3_near_masers.py
python SED_modeling.py
python new_publication_plots.py
```

Results print to the terminal. To also keep them in a file: `python ca3_near_masers.py | tee ca3_output.txt`.

Files each script saves (in the folder you run it from):
- `ca3_near_masers.py`: none (results are printed).
- `SED_modeling.py`: `CA3_HII_region_fit.png`, `d2_HII_region_fit.png`, and 10 Gaussian-fit check plots `gaussian2D_fit_<source>_<freq>GHz.png`.
- `new_publication_plots.py`: `Overview_W51IRS2_6.png` (plus `w51north_overplot_apr25.png` if you uncomment [OPTION]).

The file paths are the ones in your notebooks (`/Volumes/GalagaResearch/...`), so the external drive
must be mounted. Figures are saved to the folder you run the script from.

Requirements (already in your notebook environment): numpy, pandas, astropy, matplotlib,
spectral-cube, radio-beam, pyspeckit, and `dust_emissivity`, which `HII_model.py` imports. PyPI's
`dust_emissivity` package is broken; install it from GitHub if needed:
`pip install git+https://github.com/keflavich/dust_emissivity`.

## Suggestions (commented out)

Each suggestion is a tagged comment block. Uncomment the code lines in a block to apply it.
Descriptive lines stay as comments.

- **[A2]** uses the adopted distance, 5.41 kpc (Sato et al. 2010). au = arcsec × pc.
  In `new_publication_plots.py`, uncomment the helper in cell 2 first; the cell 11 bar lines use it.
- **[FIX]** corrects a bug found while converting:
  - `ca3_near_masers.py` cell 15: the position angle is reversed by 180°.
  - `ca3_near_masers.py` cell 13: the d2 loop covers only 7 of the 9 d2 masers. The maximum is unchanged.
- **[STYLE]** writes "au" in lowercase on the Figure 1 labels.
- **[OPTION]** also saves Figure 1 as `w51north_overplot_apr25.png`, the file name main.tex uses.
- **[NOTE]** is information only.

## Expected results with the [A2] suggestions applied

| Quantity | Old | New (5.41 kpc) |
|---|---|---|
| CA3 max separation (0.020649″) | 106 au | 112 au |
| d2 max separation (0.062995″) | 323 au | 341 au |
| W51 North, 61.45 km/s maser (0.225420″) | 1156 au, P.A. 129° | 1220 au, P.A. 309° (with [FIX]) |
| W51 North, 60.94 km/s maser (0.407092″) | 2088 au, P.A. 282° | 2202 au, P.A. 102° (with [FIX]) |
| CA3 radius (SED fit) | 48 au | ≈51 au |
| d2 radius (SED fit) | 1.8e2 au | ≈1.9e2 au (186) |
| Figure 1 bars at 0.007″/px | 7 / 39.2 / 140 px | 6.6 / 37.0 / 108.9 px |

The old separations used **5128 pc**, hardcoded in `ca3_near_masers.ipynb`. The SED fit and the
Figure 1 bars used **5.1 kpc**.

## What was tested

Each script was run without `MPLBACKEND`, confirming that the scripts themselves keep windows from opening.

- **`ca3_near_masers.py`:** run on your real CSVs from the repo, as-is and with every suggestion
  uncommented. It reproduces the paper's current numbers at 5128 pc and gives the new values above.
- **`SED_modeling.py`:** run as-is and with every suggestion uncommented. The CSV inputs were real, but
  your five continuum FITS images are on your drive, so stand-in images with the Table 2 peak fluxes
  were used. These reproduce 48 au at 5.1 kpc and give 51 au at 5.41 kpc.
- **`new_publication_plots.py`:** run as-is and with every suggestion uncommented. The maser CSV was
  real. The CASA image and the outflow FITS files were stand-ins, so the picture itself was not checked;
  the scale bars, labels and save were.

Run the scripts on your machine with the real images to get the final values.
