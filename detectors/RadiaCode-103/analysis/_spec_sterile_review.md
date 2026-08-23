OLLAMA TASK SPEC — generate ONE Python module, output ONLY raw Python code, no markdown fences.

# Target file: ollama_sterile_review.py

## Purpose
Reusable driver: send one or more source files plus a review brief to qwen3-coder:30b
via guarded_generate, asking it to find bugs (NOT to praise the code), and print the
JSON verdict to stdout. Used for IRON MODE sterile audit cycles (#FIT-1, 2026-08-22).

## Behavior — CLI script (argparse-free, sys.argv based)

Usage: python ollama_sterile_review.py <brief.md> <file1.py> [file2.py ...]

1. sys.stdout.reconfigure(encoding="utf-8"); sys.stdin.reconfigure(encoding="utf-8", errors="replace")
2. Import guarded_generate the same way gen_code.py does:
   sys.path.insert(0, path to "C:/Users/1/.claude/skills/workflow/scripts")
   from vram_guard_reference import guarded_generate
3. Read brief.md (utf-8) — this is the review instructions, written by the caller,
   containing NO conclusions, only the review criteria and acceptance rules.
4. For each source file argument: read it (utf-8), and build a combined prompt:
   the brief text, then for each file a section "=== FILE: <path> ===" followed by
   the raw file content.
5. Append fixed instructions to the prompt: respond ONLY with JSON, schema:
   {"found_issues": bool, "issues": [{"file": str, "severity": "critical"|"major"|"minor",
   "description": str, "suggested_fix": str}]}. If nothing wrong, found_issues=false and
   issues=[].
6. Call guarded_generate(model="qwen3-coder:30b", prompt=combined_prompt, fmt="json",
   want_gpu=True, priority=50, max_wait_s=600, temperature=0.4, num_ctx=16384,
   extra_options={"num_predict": 3000})
7. Parse resp["response"] as JSON (strip markdown fences if present). If parse fails,
   print to stderr "PARSE FAILED: <raw text first 500 chars>" and print
   '{"found_issues": false, "issues": [], "parse_error": true}' to stdout, exit 0.
8. On success, print the parsed JSON (json.dumps, ensure_ascii=False, indent=2) to stdout.

## Constraints
- Python 3.10, stdlib + the guarded_generate import only.
- No classes, one linear top-to-bottom script under `if __name__ == "__main__":`.
- All comments ASCII (transliterate Russian) — this file may run on Windows console.
- Must not raise on missing files — print "FILE NOT FOUND: <path>" to stderr and skip it.
