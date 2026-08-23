import sys, os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import numpy as np
import fit_coverage as fc
import fit_nuclides as fn

def main():
    names, A, e_meas, cps_meas, live = fc.assemble()
    
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        amps_dict, a_mu = fn.main()
    
    amp = np.array([amps_dict.get(n, 0.0) for n in names])
    model = A @ amp
    diff = cps_meas - model

    print("=== DEFICIT SPEKTRA (izmerenie minus model iz linij) ===")
    
    bands = ((20, 100), (100, 300), (300, 700), (700, 1500), (1500, 2000),
             (2000, 2400), (2400, 2830))
    print("%5s-%-6s %10s %10s %10s %8s" % ("band", "keV", "meas", "model", "deficit", "%def"))
    
    for e_low, e_high in bands:
        mask = (e_meas >= e_low) & (e_meas < e_high)
        meas_sum = np.sum(cps_meas[mask])
        model_sum = np.sum(model[mask])
        diff_sum = np.sum(diff[mask])
        pct = 100 * diff_sum / meas_sum if meas_sum > 0 else float('nan')
        print("%5d-%-6d %10.5f %10.5f %10.5f %8.1f%%" % (e_low, e_high, meas_sum, model_sum, diff_sum, pct))

    print("\n=== ZNACHIMYE PIKI V DEFICITE (net > 3 sigma) ===")
    found = []
    for nuc, e0 in fn.DIAG:
        net, sd = fn.net_window(diff, e_meas, e0, live)
        if sd > 0 and net > 3.0 * sd:
            found.append((nuc, e0, net, sd))
    
    if found:
        for nuc, e0, net, sd in found:
            print("%-8s %8.1f keV   net=%10.3e cps   n_sigma=%7.1f" % (nuc, e0, net, net/sd))
    else:
        print("net znachimyh pikov - deficit gladkij (kontinuum/rasseyanie/fon), ne otdelnye linii")

    # Save full diff to CSV
    out_path = os.path.join(_HERE, "..", "results", "deficit_spectrum_20260822.csv")
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write("energy_keV,diff_cps\n")
        for e, d in zip(e_meas, diff):
            f.write(f"{e:.1f},{d}\n")

    return 0

if __name__ == "__main__":
    sys.exit(main())
