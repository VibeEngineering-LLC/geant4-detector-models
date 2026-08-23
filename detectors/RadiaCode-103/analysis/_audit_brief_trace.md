You are a strict adversarial reviewer. You did NOT write this code and have not seen
any prior review. Do NOT search for generic patterns — instead MENTALLY EXECUTE the
function print_window_residuals with this concrete input and check every intermediate
number by hand:

e_meas = np.array([605.0, 607.0, 609.0, 611.0, 613.0])  # keV, 2 keV spacing
cps_meas = np.array([0.01, 0.02, 0.05, 0.02, 0.01])
models = [("A", np.array([0.02, 0.02, 0.02, 0.02, 0.02]))]
fwhm_fn = lambda e0: 0.09 * e0  # 9% FWHM, typical CsI(Tl) at 609 keV -> FWHM ~ 54.8 keV

For the window "Bi214 609" (e0=609.3): compute sigma, compute the window mask bounds
(e0 - 2.5*sigma, e0 + 2.5*sigma) numerically, state which of the 5 channels fall
inside, compute ym (sum of cps_meas in that mask) and rm (sum of model in that mask)
by hand, then compute the printed percent value by hand. Then check: does the sign and
magnitude match what the code in coverage_windows.py would actually print? Report any
discrepancy between your by-hand calculation and what the code computes — this is the
most reliable way to catch an off-by-one or a sign flip that pattern-matching misses.

Also check: for the label "mu 2700-2828", is the boundary condition (>= 2700, < 2828)
consistent with how e_meas would realistically be spaced (e.g. if a channel lands
exactly on 2828.0, is it included or excluded, and does that match the comment/label
"2700-2828" which implies an inclusive range)?

Report ONLY concrete numeric discrepancies you found by hand-tracing, not generic
advice.
