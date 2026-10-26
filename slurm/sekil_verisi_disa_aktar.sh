#!/usr/bin/env bash
# =============================================================================================
# Şekillerin ve ek tabloların okuduğu küçük tabloları sonuç ağacından ayırıp paketler.
# =============================================================================================
# Neden: makalenin şekilleri `figures/hepsini_uret.sh <kök> <çıktı>` ile üretiliyor ve kök dizinin
# altında results/ bekliyor. Sonuç ağacı 57 GB ve 500 binden fazla dosya; hepsini yerele çekmek
# anlamsız. Oysa şekiller ve ek tablolar yalnız küçük CSV/JSON tablolarını okuyor (dereceler,
# katmanlar, ölçümler, panel, CABBAGE, dış araç karşılaştırması). Bu betik onları ayırır, böylece
# Sonuçlar yazılırken her sayı için TRUBA'ya gidilmez.
#
# Ağırlığı yapan şey genom başına araç çıktıları (ResFinder BLAST dizinleri genom başına ~21 dosya,
# 2026-10-11'de dosya kotasını doldurmuştu). Onlar dışarıda kalır.
#
#   bash slurm/sekil_verisi_disa_aktar.sh olc                 # ne kadar yer tutacak, ne dışarıda
#   bash slurm/sekil_verisi_disa_aktar.sh paketle <cikti.tgz> # paketle
#
# ⚠️ Paketleme GİRİŞ DÜĞÜMÜNDE ÇALIŞMAZ: tar orada CPU sınırına takılıp SIGXCPU ile ölür. Debug
# düğümünde koş:
#   cd $AMR_WORK && srun -p debug -N1 -c8 --time=00:30:00 --pty bash
# Ölçme (olc) hafiftir, giriş düğümünde koşabilir.
set -euo pipefail

KOMUT=${1:?kullanim: olc | paketle <cikti.tgz>}
KOK=${AMR_WORK:?AMR_WORK tanimli degil}/results
[ -d "$KOK" ] || { echo "sonuc dizini yok: $KOK" >&2; exit 1; }

# Genom başına araç çıktıları: şekillerin hiçbiri okumaz, ama dosya sayısının çoğunu yaparlar
# (ilk ölçümde 65.436 dosyanın 34.000'i; yalnız E. coli'de 2 x 5.763). Şekillerin okuduğu toplu
# tablolar bir üst dizinde durur (<org>/external/comparison.csv, external_summary.json), o yüzden
# yalnız alt dizinleri dışlıyoruz.
HARIC=(-path '*/external/resfinder' -o -path '*/external/amrfinder' -o -path '*/external/rgi'
       -o -path '*_blast' -o -path '*/genomes' -o -path '*/tmp')
# Tek bir tablo bu boyutu aşıyorsa şekil verisi değildir (ör. unitig matrisi dökümü).
AZAMI_MB=${AZAMI_MB:-80}

bul() {
    find "$KOK" \( "${HARIC[@]}" \) -prune -o \
         -type f \( -name '*.csv' -o -name '*.json' -o -name '*.tsv' \) \
         -size -${AZAMI_MB}M -print
}

case "$KOMUT" in

olc)
    echo "kaynak : $KOK"
    echo "sinir  : dosya basina < ${AZAMI_MB} MB"
    echo
    mapfile -t F < <(bul)
    echo "secilen dosya : ${#F[@]}"
    printf '%s\0' "${F[@]}" | du -ch --files0-from=- 2>/dev/null | tail -1 | sed 's/^/toplam boyut  : /'
    echo
    echo "--- en cok yer tutan 15 dizin (organizma/antibiyotik/alt) ---"
    printf '%s\n' "${F[@]}" | sed "s|$KOK/||" | awk -F/ '{print (NF>3 ? $1"/"$2"/"$3 : $0)}' \
        | sort | uniq -c | sort -rn | head -15
    echo
    echo "--- en buyuk 10 dosya ---"
    printf '%s\0' "${F[@]}" | du -h --files0-from=- 2>/dev/null | sort -rh | head -10 \
        | sed "s|$KOK/||"
    echo
    echo "--- dislanan dizinlerde kalan dosya sayisi ---"
    find "$KOK" \( "${HARIC[@]}" \) -type d -printf '%p\n' 2>/dev/null \
        | while read -r d; do printf '%8d  %s\n' "$(find "$d" -type f | wc -l)" "${d/#$KOK\//}"; done \
        | sort -rn | head -8
    echo
    echo "boyut uygunsa: bash $0 paketle <cikti.tgz>   (debug dugumunde)"
    ;;

paketle)
    CIKTI=${2:?cikti dosyasi gerekli, or. $AMR_WORK/sekil_verisi.tgz}
    case "$(hostname)" in arf-ui*) echo "HATA: giris dugumunde tar SIGXCPU alir; debug dugumunde kos." >&2; exit 1;; esac
    LISTE=$(mktemp); trap 'rm -f "$LISTE"' EXIT
    bul | sed "s|^$KOK/||" > "$LISTE"
    echo "$(wc -l < "$LISTE") dosya paketleniyor -> $CIKTI"
    tar -czf "$CIKTI" -C "$KOK" --files-from="$LISTE" --transform 's|^|results/|'
    echo "bitti: $(du -h "$CIKTI" | cut -f1)  $CIKTI"
    echo
    echo "yerelde ac ve sekilleri uret:"
    echo "  mkdir -p ~/kanit_veri && tar -xzf $(basename "$CIKTI") -C ~/kanit_veri"
    echo "  bash makale/figures/hepsini_uret.sh ~/kanit_veri <cikti-dizini>"
    ;;

*) echo "bilinmeyen komut: $KOMUT (olc | paketle)" >&2; exit 2;;
esac
