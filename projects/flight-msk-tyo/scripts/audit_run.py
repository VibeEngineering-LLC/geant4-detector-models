"""Прогон команд-доказательств самоаудита (#SA-11): вывод записывает скрипт, не автор.
Использование: python audit_run.py <claims.txt> <out.md>; строка claims.txt: ID|описание утверждения|команда (bash), '#' — комментарий."""
import sys, subprocess, datetime
sys.stdout.reconfigure(encoding="utf-8")
rows = [l.rstrip("\n").split("|", 2) for l in open(sys.argv[1], encoding="utf-8") if l.strip() and not l.startswith("#")]
out = ["# Прогон команд-доказательств " + datetime.datetime.now().isoformat(timespec="seconds"), ""]
bad = 0
for cid, claim, cmd in rows:
    r = subprocess.run([__import__("os").environ.get("BASH_EXE", "bash"), "-c", cmd], capture_output=True, timeout=600)
    txt = (r.stdout + r.stderr).decode("utf-8", "replace").replace("\x00", "").strip().splitlines()
    bad += r.returncode != 0
    out += ["## %s — %s" % (cid, claim), "команда: `%s`" % cmd, "код возврата: %d" % r.returncode, "```"] + txt[-25:] + ["```", ""]
out.append("команд с ненулевым кодом: %d из %d" % (bad, len(rows)))
open(sys.argv[2], "w", encoding="utf-8", newline="\n").write("\n".join(out) + "\n")
print(out[-1])
