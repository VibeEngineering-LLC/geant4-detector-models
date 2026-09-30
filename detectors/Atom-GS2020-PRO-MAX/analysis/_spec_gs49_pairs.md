Write ONE Python 3 script, no prose, no markdown fences, only code. UTF-8, `sys.stdout.reconfigure(encoding="utf-8")`. Standard library only.

TASK: parse decay tally files and build a table of coincidence probabilities "K X-ray + gamma in the same decay".

INPUT files: `C:/g4work/gs2020/gs49/tally_<NUC>.txt` for NUC in Ra228, Ac228, Th228, Ra224, Rn220, Pb212, Bi212, Tl208 (skip missing).
Format: lines starting with '#' ignored; lines `ion_Z,<Z>`, `ion_A,<A>`, `seed,<s>`, `n_events,<N>`, `n_signatures,<k>` (key,value);
the line `count;signature`; then data lines `<count>;P:<e1> <e2> ...;E:<e1> <e2> ...` where P-list are photon energies in keV
(may be empty), E-list conversion/Auger electron kinetic energies in keV (may be empty). Split data lines on ';' into exactly three parts
(count, "P:...", "E:..."). Every data line = `count` events with exactly this photon/electron content; N = sum of counts = n_events.

PARAMETERS:
- Daughter element Z of the decay (atom where the conversion happens) DZ = {Ra228:89, Ac228:90, Th228:88, Ra224:86, Rn220:84, Pb212:83, Bi212:84 (beta branch; alpha branch gives Tl-208, Z=81, treat K-X by energy window, see below), Tl208:82}.
- K X-ray windows in keV (photon energy in the window => K X-ray of that element): element Z -> (lo, hi):
  81 Tl: 70.0-90.0; 82 Pb: 72.0-90.0; 83 Bi: 74.0-93.0; 84 Po: 76.0-96.0; 86 Rn: 81.0-101.0; 88 Ra: 85.0-107.0; 89 Ac: 88.0-110.0; 90 Th: 89.0-112.0.
  For nuclides with two daughter elements (Bi212: Po (beta) or Tl (alpha)) classify a photon as X-ray of Z if it lies in that Z window; use windows of both.
  A photon in a K-X window counts as X-ray; every other photon with energy >= 100 keV counts as gamma. (Photons <100 keV outside windows are ignored.)
- K binding energies B_K keV: 81:85.530, 82:88.0045, 83:90.526, 84:93.105, 86:98.404, 88:103.922, 89:106.755, 90:109.651.
- ICC correction table: file os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "gs49_icc_ratio_lnhb_over_g4.csv")
  (columns daughter_Z,e_tr_keV,ratio; lines starting with '#' ignored; first non-comment line is the header).

ALGORITHM per nuclide, per data line (count c):
1. Split photons into X-rays (K windows) and gammas (others >= 100 keV). If no X-ray or no gamma: nothing for pairs.
2. Correction weight w for the event: find the K-conversion electron: for each electron e in the E-list and each Z in the nuclide's daughter elements,
   E_tr = e + B_K[Z]; if a row of the ICC table with same daughter_Z has |e_tr_keV - E_tr| <= 0.6 then w = its ratio (take the closest row); otherwise w = 1.0.
   If several matches, use the one with the largest |ratio-1|. If the event has X-rays but no electrons, w = 1.0.
3. For every distinct pair (x, g) with x in the event's X-ray energies and g in the event's gamma energies (each distinct energy value once per event, ignoring duplicates):
   accumulate raw[(x,g)] += c, corr[(x,g)] += c*w.
OUTPUT: write `../audit/gs49_pairs.csv` relative to the script folder, columns:
`nuclide,e_x_keV,e_gamma_keV,e_sum_keV,n_events,p_raw,p_raw_err,p_corr` where p_raw = raw/N, p_raw_err = sqrt(raw)/N, p_corr = corr/N, e_sum = x+g rounded to 3 decimals,
energies printed with 3 decimals, sorted by nuclide then descending p_corr. Keep only rows with raw >= 10.
Also print to stdout in Russian: per nuclide N and the number of rows; the ten largest p_corr rows per nuclide; and for each nuclide the total fraction of events with at least one K-X-ray.
Wrap parsing in a function; `if __name__ == "__main__": main()`.
