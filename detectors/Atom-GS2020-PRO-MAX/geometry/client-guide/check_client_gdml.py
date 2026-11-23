import sys, math
import xml.etree.ElementTree as ET
sys.stdout.reconfigure(encoding="utf-8")
fails = 0
def ok(t): print(f"OK   {t}")
def fail(t):
    global fails; fails += 1; print(f"FAIL {t}")
def ev(e, ns):
    try: return float(eval(e, {"__builtins__": {}}, ns))
    except: return None
def get_num(el, a, ns):
    v = el.get(a)
    if v is None: return 0.0
    r = ev(v, ns)
    return 0.0 if r is None else r
def main():
    if len(sys.argv) < 2: sys.exit(1)
    try: root = ET.parse(sys.argv[1]).getroot()
    except: fail("XML parse error"); print(f"ИТОГ: FAIL={fails}"); return
    req = ["define", "materials", "solids", "structure", "setup"]
    if not all(root.find(r) is not None for r in req): fail("Missing root sections")
    ns = {"sqrt": math.sqrt, "sin": math.sin, "cos": math.cos, "pi": math.pi}
    for _ in range(10):
        for c in root.find("define").findall("constant"):
            if c.get("name") not in ns:
                v = ev(c.get("value"), ns)
                if v is not None: ns[c.get("name")] = v
    mats = {m.get("name"): m for m in root.find("materials").findall("material")}
    solids = {s.get("name"): s for s in root.find("solids")}
    vols = {v.get("name"): v for v in root.find("structure").findall("volume")}
    for v in vols.values():
        mr, sr = v.find("materialref"), v.find("solidref")
        if mr is None or sr is None: fail(f"Volume {v.get('name')} missing refs")
        elif sr.get("ref") not in solids: fail(f"Solid {sr.get('ref')} not found")
        elif mr.get("ref") not in mats and not mr.get("ref").startswith("G4_"): fail(f"Material {mr.get('ref')} not found")
    for pv in root.find("structure").iter("physvol"):
        vr = pv.find("volumeref")
        if vr is None or vr.get("ref") not in vols: fail(f"Physvol {pv.get('name')} bad ref")
    w = root.find("setup").find("world")
    if w is None or w.get("ref") not in vols: fail("World ref invalid")
    req_vols = ["LV_NaI_crystal", "LV_Reflector", "LV_Al_housing_1", "LV_Al_housing_2", "LV_Silica_window", "LV_Epoxy_seal_mid", "LV_Epoxy_seal_top", "LV_SourceMatrix", "LV_Headspace", "LV_VesselWall"]
    for rv in req_vols:
        if rv not in vols: fail(f"Missing volume {rv}")
    if "LV_NaI_crystal" in vols:
        v = vols["LV_NaI_crystal"]
        s = solids.get(v.find("solidref").get("ref"))
        if s is None or s.tag != "genericPolycone": fail("Crystal solid not genericPolycone")
        else:
            pts = [p for p in s.findall("rzpoint")]
            rs = [get_num(p, "r", ns) for p in pts]
            zs = [get_num(p, "z", ns) for p in pts]
            ok(f"Crystal: max_r={max(rs):.2f}mm, height={max(zs)-min(zs):.2f}mm")
            if v.find("materialref").get("ref") != "G4_SODIUM_IODIDE": fail("Crystal material wrong")
    if len(sys.argv) > 2:
        try: rroot = ET.parse(sys.argv[2]).getroot()
        except: fail("Ref parse error")
        else:
            rsolids = {s.get("name"): s for s in rroot.find("solids")}
            rvols = {v.get("name"): v for v in rroot.find("structure").findall("volume")}
            diffs = 0
            for name, rs in rsolids.items():
                if name == "WorldBox": continue
                if name not in solids: diffs += 1; fail(f"Missing solid {name}")
                elif rs.tag != "tessellated" and ET.tostring(rs) != ET.tostring(solids[name]): diffs += 1; fail(f"Solid {name} mismatch")
                elif rs.tag == "tessellated" and len(rs.findall("triangular")) != len(solids[name].findall("triangular")): diffs += 1; fail(f"Solid {name} tri count")
            for name, rv in rvols.items():
                if name == "World": continue
                if name not in vols: diffs += 1; fail(f"Missing vol {name}")
                else:
                    cv = vols[name]
                    if (cv.find("materialref").get("ref"), cv.find("solidref").get("ref")) != (rv.find("materialref").get("ref"), rv.find("solidref").get("ref")):
                        diffs += 1; fail(f"Vol {name} refs")
            if diffs == 0: ok(f"Reference check passed ({len(rsolids)-1} solids, {len(rvols)-1} vols)")
            else: fail(f"Reference check: {diffs} diffs")
    else: print("SKIP reference check")
    if "LV_SourceMatrix" in vols:
        v = vols["LV_SourceMatrix"]
        s = solids.get(v.find("solidref").get("ref"))
        if s is None or s.tag != "polycone": fail("SourceMatrix not polycone")
        else:
            zps = [p for p in s.findall("zplane")]
            V = 0
            for i in range(len(zps)-1):
                z1, z2 = get_num(zps[i], "z", ns), get_num(zps[i+1], "z", ns)
                r1, r2 = get_num(zps[i], "rmax", ns), get_num(zps[i+1], "rmax", ns)
                m1, m2 = get_num(zps[i], "rmin", ns), get_num(zps[i+1], "rmin", ns)
                V += math.pi * abs(z2-z1) * ((r1**2-m1**2) + (r2**2-m2**2)) / 2
            ok(f"sample volume = {V/1000:.4f} mL")
            mref = v.find("materialref").get("ref")
            if mref in mats:
                d = ev(mats[mref].find("D").get("value"), ns)
                fs = sum(float(f.get("n")) for f in mats[mref].findall("fraction"))
                if abs(fs - 1) > 1e-6: fail("Fractions sum != 1")
                else: ok(f"density={d} g/cm3, mass={d*V/1000:.4f} g")
            elif mref.startswith("G4_"): ok("density from NIST, mass not computed")
    if "LV_SourceMatrix" in vols and "LV_VesselWall" in vols and "LV_Headspace" in vols:
        sm = solids[vols["LV_SourceMatrix"].find("solidref").get("ref")]
        vw = solids[vols["LV_VesselWall"].find("solidref").get("ref")]
        hs = solids[vols["LV_Headspace"].find("solidref").get("ref")]
        sm_zs = [get_num(p, "z", ns) for p in sm.findall("zplane")]
        vw_zs = [get_num(p, "z", ns) for p in vw.findall("zplane")]
        sm_rs = [get_num(p, "rmax", ns) for p in sm.findall("zplane")]
        vw_rs = [get_num(p, "rmax", ns) for p in vw.findall("zplane")]
        hs_zs = [get_num(p, "z", ns) for p in hs.findall("zplane")]
        if min(sm_zs) < min(vw_zs) - 1e-6 or max(sm_zs) > max(vw_zs) + 1e-6: fail("SourceMatrix z out of VesselWall")
        if max(sm_rs) >= max(vw_rs): fail("SourceMatrix r >= VesselWall r")
        if abs(hs_zs[0] - sm_zs[-1]) > 1e-6: fail("Headspace z mismatch")
    print(f"ИТОГ: FAIL={fails}")
    print("Это проверка структуры. Окончательная проверка - Geant4: /geometry/test/run (накладки) и пробный расчёт.")
if __name__ == "__main__":
    main()
    sys.exit(0 if fails == 0 else 1)
