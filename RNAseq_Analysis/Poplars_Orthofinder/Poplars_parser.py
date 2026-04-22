#!/usr/bin/env python3
import os, re, sys
from pathlib import Path

# Cambia si tu carpeta está en otra ruta
FOLDER = "Single_Copy_Orthologue_Sequences"

out_lines = ["file\tHAP2_id\tALBA_id\tlen_HAP2\tlen_ALBA"]
fa_re = re.compile(r'\.(fa|fasta|faa)(\.gz)?$', re.IGNORECASE)

def clean_seq(s):
    # quita espacios, * final (stop), convierte a mayúsculas
    return re.sub(r'\*+$','', re.sub(r'\s+','', s)).upper()

def header_to_species_id(h):
    # detecta especie por patrón del ID
    h = h.strip().lstrip('>')
    if re.match(r'^(PtXaAlbH|PtremulaxPopulusalba|Potra|Potri)\.', h):
        return ('HAP2', h)
    if re.match(r'^(XP|NP)_[0-9]+\.[0-9]+$', h) or re.match(r'^PALZ_[0-9]+', h):
        return ('ALBA', h)
    # fallback: si no reconoce, devuelve None (lo manejamos luego)
    return (None, h)

base = Path(FOLDER)
if not base.is_dir():
    print(f"ERROR: no encuentro la carpeta {FOLDER}", file=sys.stderr)
    sys.exit(1)

for fa in sorted(base.glob("*")):
    if not fa_re.search(fa.name): 
        continue
    headers, seqs = [], []
    with open(fa, "r") as fh:
        cur = []
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
                headers.append(line.strip())
            else:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur))

    if len(headers) != 2 or len(seqs) != 2:
        # salta archivos que no tienen exactamente 2 secuencias
        continue

    # Limpieza secuencias (quita saltos y * final)
    seqs = [clean_seq(s) for s in seqs]
    pairs = [header_to_species_id(h) for h in headers]

    hap2_id = alba_id = None
    for (sp, hid), s in zip(pairs, seqs):
        if sp == 'HAP2' and hap2_id is None:
            hap2_id = hid; len_hap2 = len(s)
        elif sp == 'ALBA' and alba_id is None:
            alba_id = hid; len_alba = len(s)
        else:
            # si no detectó patrón, intenta asignar por descarte
            pass

    # fallback por orden si no identificó patrones
    if hap2_id is None or alba_id is None:
        # Heurística: PtXaAlbH.* = HAP2; XP_/NP_/PALZ_ = ALBA
        for (sp, hid), s in zip(pairs, seqs):
            if hap2_id is None and re.match(r'^(PtXaAlbH|PtremulaxPopulusalba|Potra|Potri)\.', hid):
                hap2_id = hid; len_hap2 = len(s)
            if alba_id is None and (re.match(r'^(XP|NP)_[0-9]+\.[0-9]+$', hid) or re.match(r'^PALZ_[0-9]+', hid)):
                alba_id = hid; len_alba = len(s)

    # si sigue sin distinguir, asigna por posición (no ideal, pero no bloquea)
    if hap2_id is None:
        hap2_id = headers[0].lstrip('>')
        len_hap2 = len(seqs[0])
    if alba_id is None:
        alba_id = headers[1].lstrip('>')
        len_alba = len(seqs[1])

    out_lines.append(f"{fa.name}\t{hap2_id}\t{alba_id}\t{len_hap2}\t{len_alba}")

# escribe salida
with open("one2one_pairs.tsv", "w") as out:
    out.write("\n".join(out_lines))

print("✅ Generado: one2one_pairs.tsv")

