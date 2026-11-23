import sys
import hashlib

def main():
    if len(sys.argv) != 4:
        print("Usage: python make_client_gdml.py <template.gdml> <device.gdml> <out.gdml>", file=sys.stderr)
        sys.exit(1)

    template_path = sys.argv[1]
    device_path = sys.argv[2]
    out_path = sys.argv[3]

    with open(device_path, "r", encoding="utf-8", newline="") as f:
        device = f.read()

    # 1. DEFINE
    start = device.find("<define>")
    end = device.find("</define>")
    if start == -1 or end == -1:
        print("Error: DEFINE not found", file=sys.stderr)
        sys.exit(1)
    define = device[start + len("<define>"):end]

    # 2. MATERIALS
    start = device.find("<materials>")
    end = device.find("</materials>")
    if start == -1 or end == -1:
        print("Error: MATERIALS not found", file=sys.stderr)
        sys.exit(1)
    materials = device[start + len("<materials>"):end]

    # 3. SOLIDS
    # Find the line containing name="WorldBox"
    worldbox_line_start = device.find('name="WorldBox"')
    if worldbox_line_start == -1:
        print("Error: WorldBox not found", file=sys.stderr)
        sys.exit(1)
    # Find the start of this line
    line_start = device.rfind("\n", 0, worldbox_line_start)
    if line_start == -1:
        line_start = 0
    else:
        line_start += 1
    # Find the end of this line
    line_end = device.find("\n", worldbox_line_start)
    if line_end == -1:
        line_end = len(device)
    # SOLIDS is from line_end to </solids>
    solids_end = device.find("</solids>")
    if solids_end == -1:
        print("Error: </solids> not found", file=sys.stderr)
        sys.exit(1)
    solids = device[line_end:solids_end]

    # 4. VOLUMES
    start = device.find("<structure>")
    if start == -1:
        print("Error: <structure> not found", file=sys.stderr)
        sys.exit(1)
    start += len("<structure>")
    end = device.find('<volume name="World">')
    if end == -1:
        print("Error: World volume not found", file=sys.stderr)
        sys.exit(1)
    volumes = device[start:end]

    # 5. PHYSVOLS
    start = device.find('<solidref ref="WorldBox"/>')
    if start == -1:
        print("Error: WorldBox solidref not found", file=sys.stderr)
        sys.exit(1)
    start += len('<solidref ref="WorldBox"/>')
    # Find the LAST </volume> before </structure>
    structure_end = device.find("</structure>")
    if structure_end == -1:
        print("Error: </structure> not found", file=sys.stderr)
        sys.exit(1)
    # Find the last </volume> before structure_end
    last_volume_end = device.rfind("</volume>", 0, structure_end)
    if last_volume_end == -1:
        print("Error: </volume> not found", file=sys.stderr)
        sys.exit(1)
    physvols = device[start:last_volume_end]

    # Read template
    with open(template_path, "r", encoding="utf-8", newline="") as f:
        template = f.read()

    # Check markers
    markers = ["@DEVICE_DEFINE@", "@DEVICE_MATERIALS@", "@DEVICE_SOLIDS@", "@DEVICE_VOLUMES@", "@DEVICE_PHYSVOLS@"]
    for marker in markers:
        count = template.count(marker)
        if count != 1:
            print(f"Error: marker {marker} count is {count}, expected 1", file=sys.stderr)
            sys.exit(1)

    # Replace markers
    result = template.replace("@DEVICE_DEFINE@", define)
    result = result.replace("@DEVICE_MATERIALS@", materials)
    result = result.replace("@DEVICE_SOLIDS@", solids)
    result = result.replace("@DEVICE_VOLUMES@", volumes)
    result = result.replace("@DEVICE_PHYSVOLS@", physvols)

    # Write output
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        f.write(result)

    # Verify
    pieces = [define, materials, solids, volumes, physvols]
    for piece in pieces:
        if piece not in result:
            print("Error: verification failed", file=sys.stderr)
            sys.exit(1)

    # Get file size
    import os
    size = os.path.getsize(out_path)

    # Get sha256 of device file
    with open(device_path, "rb") as f:
        sha256 = hashlib.sha256(f.read()).hexdigest()

    # Print results
    print(f"output: {out_path}")
    print(f"size: {size} bytes")
    print(f"define length: {len(define)}")
    print(f"materials length: {len(materials)}")
    print(f"solids length: {len(solids)}")
    print(f"volumes length: {len(volumes)}")
    print(f"physvols length: {len(physvols)}")
    print(f"device sha256: {sha256}")

if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
