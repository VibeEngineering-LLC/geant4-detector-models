Write ONE Python 3 script, no prose, no markdown fences, only code. UTF-8, `sys.stdout.reconfigure(encoding="utf-8")`. Standard library only.

TASK: turn the coincidence table of (K X-ray, gamma) pairs into sum-peak rows for a gamma-spectrum library, mapping simulated energies onto the library's own line energies.

INPUTS (paths relative to the script folder, use os.path.join(os.path.dirname(os.path.abspath(__file__)), ...)):
1. `../audit/gs49_pairs.csv` — header `nuclide,e_x_keV,e_gamma_keV,e_sum_keV,n_events,p_raw,p_raw_err,p_corr`; p_* are probabilities per decay.
2. `data/gs49_library_lines.csv` — lines starting with '#' ignored; header `nuclide,e_keV,i_pct,kind`; kind is `x` (X-ray line) or `g` (gamma line).
OUTPUT: `configs/gs49_xg_sums.json` — JSON (ensure_ascii=False, indent=1) list of objects
`{"nuclide": str, "e_x_keV": float, "e_g_keV": float, "e_sum_keV": float, "p_pair": float, "p_err": float, "note": str}`.

ALGORITHM:
For each pairs row with p_corr >= 1e-4:
 - map e_x_keV to the nearest library line with same nuclide and kind 'x' provided |difference| <= 1.0 keV; otherwise drop the row;
 - map e_gamma_keV to the nearest library line with same nuclide and kind 'g' provided |difference| <= 0.3 keV; otherwise drop the row;
 - skip the row if the mapped X-ray energy is below 50 keV.
Merge rows that map to the same (nuclide, X line, gamma line): p_pair = sum of p_corr, p_err = sqrt(sum of p_raw_err^2) (propagate the p_corr/p_raw ratio is NOT needed).
e_sum = e_x_lib + e_g_lib rounded to 3 decimals. Drop merged rows with e_sum > 2900.0 or p_pair < 2e-4.
note = "#GS-49 X+gamma: " + nuclide + " X " + f"{e_x:.3f}" + " + g " + f"{e_g:.3f}" + "; p=" + f"{p_pair:.3e}" (p per decay of the nuclide, decay-tally run, gs49_pairs.py).
Sort by descending p_pair. Print in Russian: number of input rows used, number dropped for each reason (X unmapped, gamma unmapped, X below 50 keV, below p threshold, above 2900 keV),
the number of output rows and the ten largest with nuclide, e_x, e_g, e_sum, p_pair. Create the `configs` folder if missing.
Wrap in main(); `if __name__ == "__main__": main()`.
