Write ONE Python 3 script `mutate_replay.py`. Output ONLY the code (no markdown fences, no prose). Comments in Russian, short. Start with `import sys, subprocess` and `sys.stdout.reconfigure(encoding="utf-8")`.

## Purpose
Mutation check of a C++ self-test: for each mutation, patch the header `src/replay.hh` (relative to the current directory), rebuild, run the self-test and print whether it turned RED. ALWAYS restore the original header at the end (use try/finally) and rebuild once more after restoring.

## Mutations (list of tuples name, old_text, new_text; replace the FIRST occurrence; if `old_text` is not in the file print `<name> -> ERROR: fragment not found` and continue)
1. ("M1 weight without n", "w = n * ell_cm / (L_cm * m);", "w = ell_cm / (L_cm * m);")
2. ("M2 window half as wide", "std::abs(ev[j].y_cm - zd_cm) <= r", "std::abs(ev[j].y_cm - zd_cm) <= 0.5 * r")
3. ("M3 no cut at the tube ends", "if (std::abs(zd_cm) > L_cm / 2) return false;", "")
4. ("M4 events merged across files", "have_ = true;\n        break;", "have_ = true;\n        break;") — replace this one by a mutation of the function `Fill`: old_text `if (open_ && reader_.Next(buf_)) {`, new_text `if (open_ && reader_.Next(buf_)) { buf_.evt = 7;` (all records get the same event number).

## Build and run commands (exact)
* Build: `subprocess.run(["pwsh", "-NoProfile", "-File", "scripts/build_flight.ps1"], capture_output=True)`; a build failure (`b'build rc=0'` not present in the decoded stdout) must be reported as `<name> -> BUILD FAILED` and counted as NOT red.
* Run: `subprocess.run(["pwsh", "-NoProfile", "-Command", '. "C:\\g4work\\g4setup-11.4.2.ps1"; Set-Location results\\selftest; & "C:\\g4work\\build\\flight-msk-tyo\\replay_selftest.exe"'], capture_output=True)`; decode stdout with `"utf-8", errors="replace"`, remove NUL characters, keep the lines that contain `SELFTEST` or ` FAIL`.
* A mutation is RED when the output contains `SELFTEST FAIL`; print `<name> -> RED` or `<name> -> STILL GREEN (mutation not detected)` followed by the kept lines.
* At the end print `mutants not detected: <count>` and exit with code 1 if the count is not 0, else 0.
* Read the original header once with `encoding="utf-8"`, write with `encoding="utf-8", newline="\n"`.
