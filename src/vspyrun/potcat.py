import os
import sys

def generate_potcar(inpos, potcar_dir, outpot='POTCAR'):
    # Read the POSCAR file to get the element line
    if os.path.isfile(inpos):
        with open(inpos, 'r') as f:
            lines = f.readlines()

        # Find the element line, typically the 6th line; fallback to first all-alpha token line
        try:
            elements = lines[5].split()
        except IndexError:
            elements = []
        if not elements:
            for ln in lines:
                toks = ln.split()
                if toks and all(tok.isalpha() for tok in toks):
                    elements = toks
                    break
        if not elements:
            raise ValueError("Could not determine element line from POSCAR.")
        
        # Initialize POTCAR content
        potcar_content = ""
        
        # Append the POTCAR for each element
        for element in elements:
            element_potcar_path = os.path.join(potcar_dir, element, 'POTCAR')
            if not os.path.isfile(element_potcar_path):
                raise FileNotFoundError(f"POTCAR file for element '{element}' not found in directory: {element_potcar_path}")
            
            with open(element_potcar_path, 'r') as potcar_file:
                potcar_content += potcar_file.read()

        # Write the combined POTCAR file
        with open(outpot, 'w') as output_file:
            output_file.write(potcar_content)

        print(f"POTCAR file successfully generated as '{outpot}'.")
    else:
        print(f"{inpos} file missing. Bye!")
        sys.exit(0)
