#!/bin/bash
# Makale 3'un kod oncesi uc olcumu (makale/plan_ED_m3_frekans_suzgeci_261026.md §5).
#
# Yalniz OKUR: unitig deposu, cift matrislerinin ozetleri ve assembly'ler. Kosu
# ciktisina dokunmuyor, yalniz results/tables/m3_0* yaziyor; suren kosuyla ayni anda
# gonderilebilir. Henuz yazilmamis bir matris ya da depo hata degil, atlanir.
#
#   sbatch slurm/run_m3_olcumler_env.slurm                        # butun organizmalar
#   AMR_ORGANISM=ecoli sbatch slurm/run_m3_olcumler_env.slurm     # yalniz biri
#
# Sira onemli: 0b'nin sayilari "bir unitig'in her k-mer'i unitig'in varlik vektorunu
# tasiyor" varsayimina dayaniyor, 0a da tam onu olcuyor. Bu yuzden 0a basta.
#
# Internet gerektirmiyor, tek is olarak gidiyor (kuyruk hakki 100 isten yemiyor).
#SBATCH -p barbun -N1 -c40 --mem=300G --time=06:00:00
#SBATCH --exclude=barbun45,barbun46 -o %x-%j.out -e %x-%j.err
#SBATCH --mail-type=END,FAIL --mail-user=eren0demirbas@gmail.com
set -euo pipefail
export APPTAINER_BINDPATH=/arf
W=/arf/scratch/edemirbas/amr; H=/arf/home/edemirbas/ML_AMR_Prediction_v2; SIF=$W/containers/amr.sif
cd $H
RUN="apptainer exec --no-home $SIF python -u"

# AMR_ORGANISM verilmediyse 0a panelin her organizmasini kendi icinde doluyor
ORG="${AMR_ORGANISM:-all}"

echo "== M3-0a ($ORG) =="
$RUN scripts/m3_0a_unitig_kmer_homogeneity.py --organism "$ORG" --threads 40 2>&1

echo "== M3-0b =="
$RUN scripts/m3_0b_kmer_feature_count.py ${AMR_ORGANISM:+--organism "$AMR_ORGANISM"} 2>&1

echo "== M3-0c =="
$RUN scripts/m3_0c_collapse_ratio.py 2>&1

echo "== yazilanlar =="
ls -la results/tables/m3_0* 2>&1 || echo "(dosya yok)"
echo "===== M3 OLCUMLERI BITTI ====="
