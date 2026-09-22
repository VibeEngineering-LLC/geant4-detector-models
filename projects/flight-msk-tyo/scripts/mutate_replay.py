import sys, subprocess
sys.stdout.reconfigure(encoding="utf-8")

def read_file(path):
    with open(path, "r", encoding="utf-8") as f:
        return f.read()

def write_file(path, content):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(content)

def run_build():
    result = subprocess.run(["pwsh", "-NoProfile", "-File", "scripts/build_flight.ps1"], capture_output=True)
    stdout = result.stdout.decode("utf-8", errors="replace").replace("\x00", "")
    if "build rc=0" not in stdout:
        return False
    return True

def run_test():
    result = subprocess.run(["pwsh", "-NoProfile", "-Command", '. "C:\\g4work\\g4setup-11.4.2.ps1"; Set-Location results\\selftest; & "C:\\g4work\\build\\flight-msk-tyo\\replay_selftest.exe"'], capture_output=True)
    stdout = result.stdout.decode("utf-8", errors="replace").replace("\x00", "")
    lines = [line for line in stdout.splitlines() if "SELFTEST" in line or " FAIL" in line]
    return lines

def mutate_and_test(name, old_text, new_text, file_path):
    original_content = read_file(file_path)
    try:
        if old_text not in original_content:
            print(f"{name} -> ERROR: fragment not found")
            return False
        content = original_content.replace(old_text, new_text, 1)
        write_file(file_path, content)
        if not run_build():
            print(f"{name} -> BUILD FAILED")
            return False
        lines = run_test()
        if any("SELFTEST FAIL" in line for line in lines):
            print(f"{name} -> RED")
            for line in lines:
                print(line)
            return True
        else:
            print(f"{name} -> STILL GREEN (mutation not detected)")
            for line in lines:
                print(line)
            return False
    finally:
        write_file(file_path, original_content)
        run_build()

def mutate_fill(name, file_path):
    original_content = read_file(file_path)
    try:
        old_text = "if (open_ && reader_.Next(buf_)) {"
        new_text = "if (open_ && reader_.Next(buf_)) { buf_.evt = 7;"
        if old_text not in original_content:
            print(f"{name} -> ERROR: fragment not found")
            return False
        content = original_content.replace(old_text, new_text, 1)
        write_file(file_path, content)
        if not run_build():
            print(f"{name} -> BUILD FAILED")
            return False
        lines = run_test()
        if any("SELFTEST FAIL" in line for line in lines):
            print(f"{name} -> RED")
            for line in lines:
                print(line)
            return True
        else:
            print(f"{name} -> STILL GREEN (mutation not detected)")
            for line in lines:
                print(line)
            return False
    finally:
        write_file(file_path, original_content)
        run_build()

mutations = [
    ("M1 weight without n", "w = (mk > 0) ? n * ell_cm / ((L_cm - 2 * h) * mk) : 0.0;", "w = (mk > 0) ? ell_cm / ((L_cm - 2 * h) * mk) : 0.0;"),
    ("M2 window half as wide", "std::abs(ev[j].y_cm - zd_cm) <= r", "std::abs(ev[j].y_cm - zd_cm) <= 0.5 * r"),
    ("M3 no cut at the tube ends", "if (std::abs(zd_cm) > L_cm / 2 - h) return false;", ""),
    ("M4 events merged across files", "", "")
]

file_path = "src/replay.hh"
count = 0

try:
    for name, old_text, new_text in mutations[:-1]:
        if not mutate_and_test(name, old_text, new_text, file_path):
            count += 1
    if not mutate_fill("M4 events merged across files", file_path):
        count += 1
finally:
    pass

print(f"mutants not detected: {count}")
sys.exit(1 if count != 0 else 0)
