#!/usr/bin/env bash
# =============================================================================
# SnpEff annotation + F-gene minor-variant extraction
#
# For every sample under MAFS_DIR:
#   1. annotate variant_calling/<TAG>_prot_variants.vcf.gz with SnpEff
#      -> variant_calling/<TAG>_prot_variants.ann.vcf.gz
#   2. extract F-gene variants with 0.05 <= AF <= 0.49
#      -> f_gene_variants/<TAG>_F_minor_variants.tsv
#
# The SnpEff database is built once per subtype from the reference FASTA
# plus nextclade.gff, then cached under SNPEFF_DATA.
#
# Usage:
#   bash snpeff_annotate_F.sh                 # build DB if needed, run all
#   bash snpeff_annotate_F.sh --build-only    # just build the databases
#   bash snpeff_annotate_F.sh --rebuild       # force DB rebuild
#   bash snpeff_annotate_F.sh --subtype A
#   bash snpeff_annotate_F.sh --min-af 0.05 --max-af 0.49
# =============================================================================

set -o pipefail

# -----------------------------------------------------------------------------
# PATHS
# -----------------------------------------------------------------------------
SHARE_2025="path/folder"
RUN_ROOT="path/folder"

MAFS_DIR="${RUN_ROOT}"
OUT_DIR="${RUN_ROOT}/f_gene_variants"

REF_GFF_A="${SHARE_2025}/nextclade.gff"
REF_GFF_B="${SHARE_2025}/nextcladersvb.gff"
REF_RSVA="${SHARE_2025}/minMutFinder/ref_RSVA_EPI_ISL_412866.fasta"
REF_RSVB="${SHARE_2025}/minMutFinder/ref_RSVB_2025_12_09_20.fasta"

# Where the custom SnpEff databases live
SNPEFF_DATA="${RUN_ROOT}/snpeff_data"
SNPEFF_CONFIG="${SNPEFF_DATA}/snpEff.config"

DB_A="rsva_custom"
DB_B="rsvb_custom"

MIN_AF=0.05
MAX_AF=0.49
SUBTYPE="both"
BUILD_ONLY=0
REBUILD=0
VCF_SUFFIX="_prot_variants.vcf.gz"

while [ $# -gt 0 ]; do
  case "$1" in
    --subtype)    SUBTYPE="$2"; shift 2 ;;
    --min-af)     MIN_AF="$2";  shift 2 ;;
    --max-af)     MAX_AF="$2";  shift 2 ;;
    --build-only) BUILD_ONLY=1; shift ;;
    --rebuild)    REBUILD=1;    shift ;;
    --vcf-suffix) VCF_SUFFIX="$2"; shift 2 ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
warn() { printf '[%s] WARNING: %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(date +%H:%M:%S)" "$*" >&2; exit 1; }

# -----------------------------------------------------------------------------
# PREFLIGHT
# -----------------------------------------------------------------------------
command -v snpEff  >/dev/null || die "snpEff not on PATH. conda install -c bioconda snpeff"
command -v bgzip   >/dev/null || die "bgzip not on PATH"
command -v python3 >/dev/null || die "python3 not on PATH"

[ -f "$REF_GFF_A" ] || die "RSV-A GFF not found: $REF_GFF_A"
[ -f "$REF_GFF_B" ] || die "RSV-B GFF not found: $REF_GFF_B"
[ -f "$REF_RSVA" ] || die "RSV-A reference not found: $REF_RSVA"
[ -f "$REF_RSVB" ] || die "RSV-B reference not found: $REF_RSVB"
[ -d "$MAFS_DIR" ] || die "mafs dir not found: $MAFS_DIR"

# -----------------------------------------------------------------------------
# SEQUENCE-NAME SANITY CHECK
# SnpEff silently annotates nothing if the GFF seqid does not match the FASTA
# header. This is the single most common failure mode, so check it up front.
# -----------------------------------------------------------------------------
check_names() {
  local fasta="$1" label="$2" gff="$3"
  local fa_name gff_names
  fa_name="$(grep -m1 '^>' "$fasta" | sed 's/^>//' | awk '{print $1}')"
  gff_names="$(grep -v '^#' "$gff" | awk '{print $1}' | sort -u | tr '\n' ' ')"
  printf '  %-6s FASTA seq name : %s\n' "$label" "$fa_name"
  printf '  %-6s GFF  seq names : %s\n' "$label" "$gff_names"
  if ! echo "$gff_names" | grep -qw "$fa_name"; then
    warn "$label: FASTA name '$fa_name' is NOT among the GFF seqids."
    warn "        SnpEff will annotate nothing. Fix by renaming one to match,"
    warn "        e.g.:  sed 's/^>.*/>${gff_names%% *}/' $fasta > fixed.fasta"
  fi
}

log "Checking sequence-name agreement between references and GFF"
check_names "$REF_RSVA" "RSV-A" "$REF_GFF_A"
check_names "$REF_RSVB" "RSV-B" "$REF_GFF_B"
echo

# -----------------------------------------------------------------------------
# BUILD A SNPEFF DATABASE
# -----------------------------------------------------------------------------
build_db() {
  local db="$1" fasta="$2" gff="$3"
  local dir="${SNPEFF_DATA}/${db}"

  if [ -f "${dir}/snpEffectPredictor.bin" ] && [ "$REBUILD" -eq 0 ]; then
    log "Database '${db}' already built, skipping (use --rebuild to force)"
    return 0
  fi

  log "Building SnpEff database '${db}'"
  mkdir -p "$dir"
  cp "$gff"   "${dir}/genes.gff"
  cp "$fasta" "${dir}/sequences.fa"

  # Minimal config listing our custom genomes
  if [ ! -f "$SNPEFF_CONFIG" ] || [ "$REBUILD" -eq 1 ]; then
    cat > "$SNPEFF_CONFIG" << EOF
data.dir = ${SNPEFF_DATA}
${DB_A}.genome : RSV-A custom (nextclade.gff)
${DB_B}.genome : RSV-B custom (nextcladersvb.gff)
EOF
  fi

  if snpEff build -gff3 -v -noCheckCds -noCheckProtein \
       -config "$SNPEFF_CONFIG" -dataDir "$SNPEFF_DATA" "$db" 2>&1 | tail -20; then
    log "Built '${db}'"
  else
    die "snpEff build failed for '${db}'. Check ${dir}/genes.gff formatting."
  fi
}

if [ "$SUBTYPE" = "both" ] || [ "$SUBTYPE" = "A" ]; then build_db "$DB_A" "$REF_RSVA" "$REF_GFF_A"; fi
if [ "$SUBTYPE" = "both" ] || [ "$SUBTYPE" = "B" ]; then build_db "$DB_B" "$REF_RSVB" "$REF_GFF_B"; fi
echo

[ "$BUILD_ONLY" -eq 1 ] && { log "--build-only given, stopping."; exit 0; }

# -----------------------------------------------------------------------------
# ANNOTATE EACH SAMPLE
# -----------------------------------------------------------------------------
mkdir -p "$OUT_DIR"
n_ok=0; n_skip=0
annotated_list="${OUT_DIR}/.annotated_vcfs.txt"
: > "$annotated_list"

for d in "$MAFS_DIR"/MARMS_RSV*/; do
  [ -d "$d" ] || continue
  tag="$(basename "$d")"

  case "$tag" in
    MARMS_RSVA-*) st=A; db="$DB_A" ;;
    MARMS_RSVB-*) st=B; db="$DB_B" ;;
    *) warn "$tag: cannot determine subtype, skipping"; n_skip=$((n_skip+1)); continue ;;
  esac
  if [ "$SUBTYPE" != "both" ] && [ "$SUBTYPE" != "$st" ]; then continue; fi

  vcf="${d}variant_calling/${tag}${VCF_SUFFIX}"
  if [ ! -f "$vcf" ]; then
    warn "$tag: $(basename "$vcf") not found, skipping"
    n_skip=$((n_skip+1)); continue
  fi

  ann="${d}variant_calling/${tag}_prot_variants.ann.vcf"
  stats="${d}variant_calling/${tag}_snpEff_summary.html"

  log "SnpEff  $tag"
  if snpEff ann -config "$SNPEFF_CONFIG" -dataDir "$SNPEFF_DATA" \
       -stats "$stats" -noLog "$db" "$vcf" > "$ann" 2>"${ann}.err"; then
    bgzip -f "$ann"
    echo "${ann}.gz" >> "$annotated_list"
    n_ok=$((n_ok+1))
    rm -f "${ann}.err"
  else
    warn "$tag: snpEff failed, see ${ann}.err"
    n_skip=$((n_skip+1))
  fi
done

echo
log "Annotated ${n_ok} sample(s), skipped ${n_skip}"
[ "$n_ok" -eq 0 ] && die "Nothing annotated."

# -----------------------------------------------------------------------------
# EXTRACT F-GENE MINOR VARIANTS
# -----------------------------------------------------------------------------
echo
log "Extracting F-gene variants with ${MIN_AF} <= AF <= ${MAX_AF}"

python3 - "$MAFS_DIR" "$OUT_DIR" "$MIN_AF" "$MAX_AF" << 'PYEOF'
import csv, gzip, re, sys
from pathlib import Path

mafs, outdir, min_af, max_af = Path(sys.argv[1]), Path(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])

ANN_FIELDS = ['Allele','Annotation','Annotation_Impact','Gene_Name','Gene_ID',
              'Feature_Type','Feature_ID','Transcript_BioType','Rank','HGVS_c',
              'HGVS_p','cDNA_pos','CDS_pos','AA_pos','Distance','ERRORS']

AA3TO1 = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E',
          'Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F',
          'Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V','Ter':'*'}

COLS = ['sample','subtype','chrom','pos','ref','alt','af','dp','gene','feature_id',
        'annotation','impact','hgvs_c','hgvs_p','aa_change','aa_pos']

def short_aa(h):
    if not h: return '', ''
    m = re.match(r'^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', h)
    if m:
        r,p,a = m.groups()
        return f'{AA3TO1.get(r,r)}{p}{AA3TO1.get(a,a)}', p
    m = re.match(r'^p\.([A-Z][a-z]{2})(\d+)=$', h)
    if m:
        r,p = m.groups()
        return f'{AA3TO1.get(r,r)}{p}{AA3TO1.get(r,r)}', p
    m = re.search(r'(\d+)', h)
    return h, (m.group(1) if m else '')

def info_dict(s):
    d = {}
    for it in s.split(';'):
        if '=' in it:
            k,v = it.split('=',1); d[k]=v
        elif it: d[it]=True
    return d

def get_af(info, fk, sv):
    if 'AF' in info:
        try: return float(str(info['AF']).split(',')[0])
        except ValueError: pass
    if fk and sv and 'AF' in fk:
        try: return float(sv[fk.index('AF')].split(',')[0])
        except (ValueError, IndexError): pass
    if 'AO' in info and 'DP' in info:
        try: return float(str(info['AO']).split(',')[0])/float(info['DP'])
        except (ValueError, ZeroDivisionError): pass
    return None

combined = {'A': [], 'B': []}
gene_names_seen = set()
totals = {'variants':0, 'with_ann':0, 'in_F':0, 'passed':0}

for vcf in sorted(mafs.glob('MARMS_RSV*-*/variant_calling/*_prot_variants.ann.vcf.gz')):
    tag = vcf.parent.parent.name
    st = 'A' if 'RSVA' in tag else 'B'
    sample = tag.split('-', 1)[1] if '-' in tag else tag
    rows = []

    with gzip.open(vcf, 'rt') as fh:
        for line in fh:
            if line.startswith('#'): continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8: continue
            chrom,pos,_,ref,alt,_q,_fl,info_s = f[:8]
            fk = f[8].split(':') if len(f) > 8 else []
            sv = f[9].split(':') if len(f) > 9 else []
            totals['variants'] += 1

            info = info_dict(info_s)
            if 'ANN' not in info: continue
            totals['with_ann'] += 1

            af = get_af(info, fk, sv)
            if af is None: continue
            dp = info.get('DP','')

            for entry in info['ANN'].split(','):
                parts = entry.split('|')
                parts += ['']*(len(ANN_FIELDS)-len(parts))
                ann = dict(zip(ANN_FIELDS, parts))
                gene = ann.get('Gene_Name','')
                if gene: gene_names_seen.add(gene)
                if gene.upper() != 'F': continue
                totals['in_F'] += 1
                if not (min_af <= af <= max_af): break
                totals['passed'] += 1
                aac, aap = short_aa(ann.get('HGVS_p',''))
                rows.append({'sample':sample,'subtype':st,'chrom':chrom,'pos':pos,
                             'ref':ref,'alt':alt,'af':f'{af:.4f}','dp':dp,
                             'gene':gene,'feature_id':ann.get('Feature_ID',''),
                             'annotation':ann.get('Annotation',''),
                             'impact':ann.get('Annotation_Impact',''),
                             'hgvs_c':ann.get('HGVS_c',''),'hgvs_p':ann.get('HGVS_p',''),
                             'aa_change':aac,'aa_pos':aap})
                break

    sub = outdir / f'RSV{st}'
    sub.mkdir(parents=True, exist_ok=True)
    with open(sub / f'{tag}_F_minor_variants.tsv','w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)
    combined[st].extend(rows)
    print(f'  {tag}: {len(rows)} F variants in range')

print()
for st in ('A','B'):
    if not combined[st]: continue
    p = outdir / f'ALL_RSV{st}_F_minor_variants.tsv'
    with open(p,'w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(combined[st])
    n = len({r['sample'] for r in combined[st]})
    print(f'RSV-{st}: {len(combined[st])} variants across {n} samples -> {p}')

print()
print(f"Totals: {totals['variants']} variants, {totals['with_ann']} annotated, "
      f"{totals['in_F']} in F, {totals['passed']} within AF range")
if gene_names_seen:
    print('Gene names present in the annotations:', ', '.join(sorted(gene_names_seen)))
    if not any(g.upper()=='F' for g in gene_names_seen):
        print("WARNING: no gene called 'F'. Check what nextclade.gff calls it and")
        print("         adjust the gene filter accordingly.")
else:
    print("WARNING: no gene names found at all - SnpEff annotated nothing.")
    print("         Almost always a chromosome-name mismatch between the VCF")
    print("         and the GFF. Compare: zcat <vcf> | grep -v '^#' | cut -f1 | head -1")
PYEOF

echo
log "Done. Output in ${OUT_DIR}"