import sys
import os
import json

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(1)

    sys.stdout.reconfigure(encoding="utf-8")
    sys.stdin.reconfigure(encoding="utf-8", errors="replace")

    # Add the path to the workflow scripts directory
    workflow_path = "C:/Users/1/.claude/skills/workflow/scripts"
    if os.path.exists(workflow_path):
        sys.path.insert(0, workflow_path)
        from vram_guard_reference import guarded_generate
    else:
        print("ERROR: Workflow scripts path not found", file=sys.stderr)
        sys.exit(1)

    brief_file = sys.argv[1]
    source_files = sys.argv[2:]

    try:
        with open(brief_file, "r", encoding="utf-8") as f:
            brief_text = f.read()
    except Exception as e:
        print(f"ERROR: Failed to read brief file {brief_file}: {e}", file=sys.stderr)
        sys.exit(1)

    combined_prompt = brief_text + "\n\n"

    for file_path in source_files:
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                content = f.read()
            combined_prompt += f"=== FILE: {file_path} ===\n{content}\n\n"
        except Exception as e:
            print(f"FILE NOT FOUND: {file_path}", file=sys.stderr)
            continue

    fixed_instructions = (
        "Respond ONLY with JSON, schema:\n"
        '{"found_issues": bool, "issues": [{"file": str, "severity": "critical"|"major"|"minor", '
        '"description": str, "suggested_fix": str}]}. If nothing wrong, found_issues=false and issues=[].'
    )

    combined_prompt += fixed_instructions

    try:
        resp = guarded_generate(
            model="qwen3-coder:30b",
            prompt=combined_prompt,
            fmt="json",
            want_gpu=True,
            priority=50,
            max_wait_s=600,
            temperature=0.4,
            num_ctx=16384,
            extra_options={"num_predict": 3000}
        )
    except Exception as e:
        print(f"ERROR: guarded_generate failed: {e}", file=sys.stderr)
        sys.exit(1)

    response_text = resp["response"].strip()
    
    # Remove markdown code fences if present
    if response_text.startswith("```json"):
        response_text = response_text[7:]
    if response_text.endswith("```"):
        response_text = response_text[:-3]
    response_text = response_text.strip()

    try:
        parsed = json.loads(response_text)
    except Exception as e:
        print(f"PARSE FAILED: {response_text[:500]}", file=sys.stderr)
        result = {"found_issues": False, "issues": [], "parse_error": True}
        print(json.dumps(result, ensure_ascii=False, indent=2))
        sys.exit(0)

    print(json.dumps(parsed, ensure_ascii=False, indent=2))
