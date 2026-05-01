## Visualise tRNA structures using rnaplot with .fa input files from MiTFi (python script)
# usage python3 draw_trna.py <species.struct.fa>

import os
import sys
import re
import subprocess
import tempfile
import math

def parse_struct_fa(filepath):
    trnas = []
    with open(filepath) as f:
        lines = [l.rstrip() for l in f if l.strip()]
    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            header = lines[i]
            seq = lines[i+1] if i+1 < len(lines) else ''
            db  = lines[i+2] if i+2 < len(lines) else ''
            db  = re.sub(r'[^.()\[\]]', '', db)
            minlen = min(len(seq), len(db))
            trnas.append((header, seq[:minlen], db[:minlen]))
            i += 3
        else:
            i += 1
    return trnas

def get_tRNA_name(header):
    parts = header.split('|')
    aa_map = {
        'A':'Alanine','C':'Cysteine','D':'Aspartic acid','E':'Glutamic acid',
        'F':'Phenylalanine','G':'Glycine','H':'Histidine','I':'Isoleucine',
        'K':'Lysine','L1':'Leucine-1','L2':'Leucine-2','M':'Methionine',
        'N':'Asparagine','P':'Proline','Q':'Glutamine','R':'Arginine',
        'S1':'Serine-1','S2':'Serine-2','T':'Threonine','V':'Valine',
        'W':'Tryptophan','Y':'Tyrosine'
    }
    if len(parts) >= 7:
        return aa_map.get(parts[6], parts[6])
    return header

def get_species_name(filepath):
    base = os.path.splitext(os.path.basename(filepath))[0]
    name = re.sub(r'[-_](EL|RR|PAL)\w+.*$', '', base)
    return name.replace('-', ' ')

def rnaplot_svg(seq, db, name, workdir, idx=0):
    safe_name = re.sub(r'[^a-zA-Z0-9_-]', '_', name) + f'_{idx}'
    infile = os.path.join(workdir, f'{safe_name}.fa')
    with open(infile, 'w') as f:
        f.write(f'>{safe_name}\n{seq}\n{db}\n')
    subprocess.run(
        ['RNAplot', '-f', 'svg', '-t', '1', '-i', infile],
        cwd=workdir,
        capture_output=True
    )
    outfile = os.path.join(workdir, f'{safe_name}_ss.svg')
    if os.path.exists(outfile):
        return outfile
    return None

def post_process_svg(svgfile):
    with open(svgfile, 'r', encoding='utf-8') as f:
        content = f.read()

    # remove xml declaration and script blocks
    content = re.sub(r'<\?xml[^>]+\?>', '', content).strip()
    content = re.sub(r'<script[^>]*>.*?</script>', '', content, flags=re.DOTALL)

    # make backbone lighter grey
    content = re.sub(r'(stroke:\s*)grey', r'\1#d0d0d0', content)

    # replace red with soft faded blue
    content = re.sub(r'(stroke:\s*)red', r'\1#7aade0', content)

    # extract inner content
    inner = re.sub(r'^<svg[^>]*>', '', content).strip()
    inner = re.sub(r'</svg>\s*$', '', inner).strip()

    return inner

def make_grid(svg_files, labels, output_path, ncols=8, species_name=''):
    src_size  = 452
    cell_size = 200
    scale     = cell_size / src_size
    label_h   = 20
    title_h   = 50
    pad       = 15

    n = len(svg_files)
    nrows = math.ceil(n / ncols)
    total_w = pad + ncols * (cell_size + pad)
    total_h = title_h + nrows * (cell_size + label_h + pad) + pad

    lines = []
    lines.append('<?xml version="1.0" encoding="UTF-8"?>')
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" '
                 f'width="{total_w}" height="{total_h}" version="1.1">')
    lines.append(f'  <rect width="{total_w}" height="{total_h}" fill="white"/>')
    lines.append(f'  <text x="{total_w//2}" y="35" text-anchor="middle" '
                 f'font-size="20" font-style="italic" font-family="Arial" '
                 f'fill="black">{species_name}</text>')

    for idx, (svgf, label) in enumerate(zip(svg_files, labels)):
        row, col = divmod(idx, ncols)
        x = pad + col * (cell_size + pad)
        y = title_h + row * (cell_size + label_h + pad)

        if svgf and os.path.exists(svgf):
            inner = post_process_svg(svgf)
            lines.append(f'  <g transform="translate({x},{y}) scale({scale:.6f})">')
            lines.append(inner)
            lines.append('  </g>')

        lines.append(f'  <text x="{x + cell_size//2}" y="{y + cell_size + 15}" '
                     f'text-anchor="middle" font-size="11" font-family="Arial" '
                     f'fill="black">{label}</text>')

    lines.append('</svg>')

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f'Saved: {output_path}')

def process_species(filepath):
    trnas = parse_struct_fa(filepath)
    if not trnas:
        print(f'No tRNAs found in {filepath}')
        return

    species = get_species_name(filepath)
    base = os.path.splitext(filepath)[0]
    output_svg = base + '_tRNA_grid.svg'

    with tempfile.TemporaryDirectory() as workdir:
        svg_files = []
        labels = []
        for header, seq, db in trnas:
            name = get_tRNA_name(header)
            print(f'  Plotting {name}...')
            svgf = rnaplot_svg(seq, db, name, workdir, idx=len(svg_files))
            svg_files.append(svgf)
            labels.append(name)
        make_grid(svg_files, labels, output_svg, ncols=8, species_name=species)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python3 draw_trna_v3.py <species.struct.fa>')
        sys.exit(1)
    process_species(sys.argv[1])
