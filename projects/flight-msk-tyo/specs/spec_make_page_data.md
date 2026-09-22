Write ONE Python 3 script `make_page_data.py` (numpy allowed; standard library `subprocess`, `csv`, `json`, `os`, `sys`, `tempfile`). Output ONLY the code (no markdown fences, no prose). Comments in Russian, short. Start with `sys.stdout.reconfigure(encoding="utf-8")`.

## Purpose
Prepare the data file `data.js` for a static web page that shows a simulated gamma spectrum with its components.

## Command line
`make_page_data.py <merged_prefix> <lines_prefix> <out_data_js> [--emin 20] [--emax 10000]`
where `<merged_prefix>` is the prefix of merged simulation files:
* `<merged_prefix>_total_all.csv` and `<merged_prefix>_total_<comp>.csv` for comp in `neutron, proton, mup, mum, em, ep, gamma` (format: line 1 = `# ...`, line 2 = header `bin_keV,sumw,sumw2,light_sumw,light_sumw2`, then N+1 rows (last = overflow); `sumw` = rate in s^-1 per 1 keV bin, `sumw2` = squared uncertainty in s^-2);
* `<merged_prefix>_cat_all.csv`: header `pg,proc,vol,pgname,procname,volname,rate_total,bin_0,...,bin_N`, one row per category, values = rate in s^-1 per 1 keV bin.
`<lines_prefix>` is the prefix of `<lines_prefix>_lines.csv` (header `E_keV,reaction,place_class,net_cph,sigma_cph,signif,found,origin`) and `<lines_prefix>_unlisted.csv` (header `E_keV,net_cph,sigma_cph,signif,origin`).

## Steps
1. Smearing is done by an existing script `spec_smear.py` located in the same directory as this script: `python spec_smear.py <in_total.csv> <out.csv> --emin <emin> --emax <emax>` where `<in_total.csv>` has the format above; its output CSV has the header `E_keV,counts_per_h_per_keV,sigma_counts_per_h_per_keV`. Call it with `subprocess.run([sys.executable, path_to_spec_smear, ...], capture_output=True)` with the environment variable `PYTHONIOENCODING=utf-8`; check `returncode == 0`, else print stderr and `sys.exit(2)`.
2. Series: (a) `total` from `_total_all.csv` (label "Полный спектр", group "total", keep the sigma column); (b) the seven components from their `_total_<comp>.csv` with labels `neutron`->"нейтроны", `proton`->"протоны", `mup`->"мюоны +", `mum`->"мюоны −", `em`->"электроны", `ep`->"позитроны", `gamma`->"гамма-кванты", group "component"; (c) the 12 categories of `_cat_all.csv` with the largest sum of the bins 0..N-1 restricted to `emin <= bin < emax`; for each write a temporary total-format file (`sumw` = the category bins, `sumw2` = 0, other columns 0, overflow row = last bin) into a `tempfile.TemporaryDirectory()`, smear it the same way; label = `"<pgname> / <procname> / <volname>"`, group "origin"; (d) one more series "прочие происхождения" = total smeared minus the sum of the 12 origin series (clip at 0).
3. Round the y values to 4 significant digits (`float("%.4g" % v)`); the energy grid is `E0 = emin + 0.5`, step 1 keV (store `E0` and `dE`, not the whole array).
4. Lines: read `_lines.csv` with the `csv` module (fields may contain quotes; the last field `origin` may be empty): keep rows with `found == 1`, at most 80 with the largest `signif`, each as `{"E": float, "reaction": str, "net": float, "sig": float, "signif": float, "origin": str}`; from `_unlisted.csv` keep at most 30 rows with the largest `signif` as `{"E":..., "net":..., "sig":..., "signif":..., "origin":...}` in a second list.
5. Write `out_data_js` as UTF-8 text `window.SPECTRA = <json>;` where json = `{"E0": ..., "dE": 1, "unit": "отсчётов в час на кэВ", "series": [{"id": ..., "label": ..., "group": ..., "y": [...], "sig": [...] (only for total)}], "lines": [...], "unlisted": [...]}` written with `json.dump(..., ensure_ascii=False)`; print the file size in KB and the number of series.
