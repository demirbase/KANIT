#!/usr/bin/env python3
"""makale-yazimi denetleyicisi.

Kullanım:
    python3 denetle.py DOSYA.md --dil tr|en [--bolum "Introduction"] [--yalniz-unlem]

Yalnız Python standart kütüphanesiyle çalışır. Kural verilerini skill klasöründen okur:
    kurallar/bos_kaliplar.md  (boş kalıp, muaf, nicelik tabloları)
    sozluk.csv                (not sütunundaki 'kaçın-tr:', 'kaçın-en:', 'muaf-en:' alanları)

Düzeyler:  !  düzelt      ?  bak (bağlama göre karar ver)
"""
import argparse
import csv
import os
import re
import statistics
import sys
from collections import Counter

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ESIK = {  # (uyarı, bölünmeli)
    "tr": (25, 32),
    "en": (30, 35),
}
TEKDUZE_CV = 0.20        # paragrafta cümle uzunluklarının değişkenlik katsayısı alt sınırı
TEKDUZE_MIN_CUMLE = 4

# ---------------------------------------------------------------- kural verisi

def tablo_oku(yol, ad):
    """bos_kaliplar.md içindeki <!-- DENETIM-TABLOSU: ad --> bloğunu satır sözlükleri olarak döndürür."""
    metin = open(yol, encoding="utf-8").read()
    m = re.search(r"<!-- DENETIM-TABLOSU: %s -->(.*?)<!-- /DENETIM-TABLOSU -->" % ad, metin, re.S)
    if not m:
        return []
    satirlar = [s for s in m.group(1).strip().split("\n") if s.startswith("|")]
    def hucre(s):
        s = s.strip().strip("|").replace("\\|", "\x00")
        return [h.strip().replace("\x00", "|") for h in s.split("|")]
    basliklar = hucre(satirlar[0])
    sonuc = []
    for s in satirlar[2:]:
        h = hucre(s)
        d = dict(zip(basliklar, h))
        if "regex" in d:
            d["regex"] = d["regex"].strip("`")
        sonuc.append(d)
    return sonuc


def sozluk_oku(yol):
    kacin = {"tr": [], "en": []}
    muaf = {"tr": [], "en": []}
    with open(yol, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            for parca in (r.get("not") or "").split(";"):
                parca = parca.strip()
                m = re.match(r"(kaçın|muaf)-(tr|en):\s*(.+)$", parca)
                if not m:
                    continue
                tur, dil, rx = m.groups()
                if tur == "kaçın":
                    dogru = r["turkce"] if dil == "tr" else r["ingilizce"]
                    kacin[dil].append((re.compile(r"(?<![\w-])(?:%s)(?![\w-])" % rx, re.I), dogru))
                else:
                    muaf[dil].append(re.compile(rx, re.I))
    return kacin, muaf

# ---------------------------------------------------------------- dile özgü kurallar

KAYNAK_ISTEYEN = {
    "en": r"\b(studies|reports|work|research|analyses|authors)\s+(have\s+|has\s+)?(shown|demonstrated|reported|suggested|found|indicated)\b"
          r"|\bit\s+(is|has\s+been)\s+(well\s+|widely\s+)?(known|established|reported|shown|demonstrated|accepted)\b"
          r"|\b(widely|generally|commonly)\s+(used|accepted|known|reported|believed|recognised|recognized)\b"
          r"|\bprevious(ly)?\s+(studies|work|reported|shown|described)\b"
          r"|\bin\s+the\s+literature\b|\brecent\s+(studies|work|reports)\b"
          r"|\b(has|have)\s+been\s+(shown|reported|described|proposed|demonstrated)\b",
    "tr": r"çalışmalar\w*\s+(\w+\s+){0,6}göster\w*|(bildiril|gösteril|rapor\s+edil|öne\s+sürül|tanımlan)mişti?r?\b"
          r"|literatürde|(yaygın|genel)\s+olarak\s+(kullanıl|kabul|bilin)\w*|son\s+yıllarda|önceki\s+çalışmalar",
}
NEDENSEL = {
    "en": r"\b(caus(e|es|ed|ing)|confer(s|red|ring)?|prov(e|es|ed|en|ing)|is\s+responsible\s+for|are\s+responsible\s+for)\b",
    "tr": r"neden\s+ol\w*|direnç\s+kazandır\w*|sorumlu(dur|dür)|kanıtla(r|maktadır|mıştır)\b",
}
YAZAR_KALIBI = {
    "en": r"\b[A-Z][\w'-]+\s+(and\s+colleagues|and\s+co-workers|et\s+al\.?)",
    "tr": r"\b[A-ZÇĞİÖŞÜ][\w'-]+\s+(ve\s+(arkadaşları|ark\.|diğ\.|meslektaşları)|et\s+al\.?)",
}
ONCELIK = {
    "en": r"\b(for\s+the\s+first\s+time|the\s+first\s+(study|report|time|method|framework|tool|to)"
          r"|(our|this|a|an)\s+(\w+\s+)?novel\s+(method|approach|framework|tool|pipeline|strategy|resource|knowledge\s+base))\b",
    "tr": r"\bilk\s+(kez|defa|kez\s+olarak|çalışma)\b|\byeni\s+bir\s+(yöntem|yaklaşım|çerçeve)\b",
}
EKLEME_BAGLACI = {
    "en": r"^(furthermore|moreover|additionally|in\s+addition|also|besides)\b",
    "tr": r"^(ayrıca|bunun\s+yanı\s+sıra|buna\s+ek\s+olarak|ek\s+olarak|öte\s+yandan|bu\s+bağlamda|bunun\s+yanında)\b",
}
OZET_KAPANIS = {
    "en": r"^(in\s+summary|overall|taken\s+together|collectively|in\s+conclusion|in\s+short|together,\s+these)\b",
    "tr": r"^(özetle|kısacası|sonuç\s+olarak|genel\s+olarak|birlikte\s+değerlendirildiğinde|toparlarsak)\b",
}
SIFAT_UCLU = {
    "en": r"\b(\w+(?:ive|al|ous|ble|ent|ant|ic|ful|less|ary|ate))\s*,\s*(\w+(?:ive|al|ous|ble|ent|ant|ic|ful|less|ary|ate))\s*,?\s+(?:and|or)\s+(\w+(?:ive|al|ous|ble|ent|ant|ic|ful|less|ary|ate))\b",
    "tr": r"\b(\w+(?:lı|li|lu|lü|sal|sel|ik|ık|ebilir|abilir|ilir|ılır))\s*,\s*(\w+(?:lı|li|lu|lü|sal|sel|ik|ık|ebilir|abilir|ilir|ılır))\s*,?\s+(?:ve|ya\s+da|veya)\s+(\w+(?:lı|li|lu|lü|sal|sel|ik|ık|ebilir|abilir|ilir|ılır))\b",
}
BIRINCI_SAHIS = {
    "en": r"\b(we|our|us|here|this\s+study|this\s+work)\b",
    "tr": r"\b(biz|bizim|bu\s+çalışma\w*|burada)\b|\w+(d|t)(ık|ik|uk|ük)\b|\w+(ımız|imiz|umuz|ümüz)\w*",
}
# İngiliz yazımı: ABD biçimleri
ABD_YAZIM = [
    (r"\b\w*(?<!s)(?<!pr)(?<!se)(?<!ca)i[z](e|es|ed|ing|ation|ations|er|ers)\b", "-ise/-isation"),
    (r"\b(analyz|paralyz|catalyz)\w*", "-yse (analyse)"),
    (r"\b(colou?r|behavio|neighbo|favo|labo|tumo|harbo|hono)r(s|ed|ing|al)?\b", "-our"),
    (r"\b(model|label|travel|cancel|signal)(ed|ing)\b", "-lled/-lling"),
    (r"\bcatalogs?\b", "catalogue"),
    (r"\bcent(er|ers|ered)\b", "centre"),
]
ABD_MUAF = {"size", "sizes", "sized", "prize", "seize", "capsize", "resize", "resized", "citizen",
            "citizens", "colour", "behaviour", "neighbour", "favour", "labour", "tumour", "harbour",
            "honour", "colours", "behaviours", "neighbours", "neighbouring", "modelled", "modelling",
            "labelled", "labelling", "travelled", "cancelled", "signalled", "signalling", "parallel"}

# "aynı başlangıç" kuralında sayılmayan işlev sözcükleri
ISLEV_SOZCUK = {
    "en": {"the", "a", "an", "in", "of", "to", "for", "on", "at", "by", "with", "as", "from"},
    "tr": {"bir", "ve", "ile", "de", "da"},
}

KISALTMA_MUAF = {"SD", "CC", "BY", "SHA", "TODO", "QQ", "DNA", "RNA", "WHO", "PCR", "USA", "UK", "EU", "ID", "CSV", "URL", "DOI", "API",
                 "SI", "AND", "OR", "NOT", "IN", "SQL", "PDF", "NOTE", "AUC", "TR", "EN",
                 "KAYNAK", "SAYI", "KARAR", "NOT"}      # taslaktaki açık iş işaretleri

NOKTA_KORUMA = [
    r"\be\.g\.", r"\bi\.e\.", r"\bet al\.", r"\bspp?\.", r"\bFigs?\.", r"\bvs\.", r"\bcf\.",
    r"\bapprox\.", r"\bca\.", r"\bno\.", r"\bDr\.", r"\bProf\.", r"\bSt\.", r"\bsubsp\.", r"\bvar\.",
    r"\bör\.", r"\bvb\.", r"\bvd\.", r"\bbkz\.", r"\bŞek\.", r"\bsy\.", r"\bed\.", r"\bEq\.",
    r"\b[A-ZÇĞİÖŞÜ]\.",          # tek harfli baş harf: E. coli, K. pneumoniae
]

# ---------------------------------------------------------------- markdown → paragraflar

def bolum_kes(satirlar, bolum):
    """--bolum verilmişse başlığı eşleşen bölümün satır aralığını döndürür (1 tabanlı)."""
    if not bolum:
        return 1, len(satirlar)
    bas = None
    seviye = None
    for i, s in enumerate(satirlar, 1):
        m = re.match(r"^(#+)\s+(.*)", s)
        if not m:
            continue
        if bas is None and bolum.lower() in m.group(2).lower():
            bas, seviye = i, len(m.group(1))
        elif bas is not None and len(m.group(1)) <= seviye:
            return bas, i - 1
    if bas is None:
        sys.exit("Bölüm bulunamadı: %r" % bolum)
    return bas, len(satirlar)


def paragraflar(satirlar, bas, son):
    """(ilk_satir, [satır metinleri], başlık) üçlüleri."""
    sonuc, cur, cur_bas = [], [], None
    kod = False
    atla = False
    baslik = ""
    def kapat():
        nonlocal cur, cur_bas
        if cur:
            sonuc.append((cur_bas, cur, baslik))
        cur, cur_bas = [], None
    for i in range(bas, son + 1):
        s = satirlar[i - 1]
        if s.strip().startswith("```"):
            kod = not kod
            kapat()
            continue
        if kod:
            continue
        m = re.match(r"^#+\s+(.*)", s)
        if m:
            kapat()
            baslik = m.group(1)
            atla = bool(re.match(ATLANAN_BOLUM, baslik, re.I))
            continue
        if atla:
            continue
        if (not s.strip() or s.lstrip().startswith(">") or s.lstrip().startswith("|")
                or re.match(r"^\s*(---+|\*\*\*+)\s*$", s) or s.lstrip().startswith("<!--")):
            kapat()
            continue
        if re.match(r"^\s*([-*+]|\d+\.)\s+", s):
            kapat()
            s = re.sub(r"^\s*([-*+]|\d+\.)\s+", lambda x: " " * len(x.group(0)), s)
        if cur_bas is None:
            cur_bas = i
        cur.append(s)
    kapat()
    return sonuc


ATLANAN_BOLUM = r"^(\d+\.?\s*)?(references|reference list|bibliography|kaynakça|kaynaklar)\b"
SAYI_KAYNAK_BOLUMU = r"introduction|discussion|background|giriş|tartışma|arka plan|impact"

ATIF_RX = re.compile(r"\[[^\[\]]*?(\d{4}|KAYNAK\?|\d+(\s*[,–-]\s*\d+)*)[^\[\]]*?\]")
YER_TUTUCU_RX = re.compile(r"\{[\w.]+\}")
KOD_RX = re.compile(r"`[^`]*`")


def maskele(metin):
    """Uzunluğu koruyarak: kod ve yer tutucu → 'N___', atıf → boşluk, vurgu işaretleri → boşluk.

    Kod tek sözcük sayılır ve büyük harfle başlar, böylece `kod` ile başlayan bir cümlenin
    başlangıcı görülür (boşluğa çevrilseydi önceki cümleyle birleşirdi)."""
    m = KOD_RX.sub(lambda x: "N" + "_" * (len(x.group(0)) - 1), metin)
    m = ATIF_RX.sub(lambda x: " " * len(x.group(0)), m)
    m = YER_TUTUCU_RX.sub(lambda x: "N" + "_" * (len(x.group(0)) - 1), m)
    m = re.sub(r"(?<!\w)[*_]{1,2}|[*_]{1,2}(?!\w)", lambda x: " " * len(x.group(0)), m)
    return m


def cumlelere_bol(maskeli):
    korunan = set()
    for rx in NOKTA_KORUMA:
        for m in re.finditer(rx, maskeli):
            korunan.add(m.end() - 1)
    sinirlar, bas = [], 0
    for m in re.finditer(r"[.!?](?=\s+[\"'(\[A-ZÇĞİÖŞÜ0-9N]|\s*$)", maskeli):
        if m.start() in korunan:
            continue
        sinirlar.append((bas, m.end()))
        bas = m.end()
    if maskeli[bas:].strip():
        sinirlar.append((bas, len(maskeli)))
    # baştaki boşlukları at
    sonuc = []
    for b, s in sinirlar:
        while b < s and maskeli[b].isspace():
            b += 1
        if b < s:
            sonuc.append((b, s))
    return sonuc


SOZCUK_RX = re.compile(r"[\w][\w'’\-–]*", re.U)


def sozcuk_say(t):
    return len([w for w in SOZCUK_RX.findall(t) if re.search(r"[^\W_]", w) or w.startswith("N_")])

# ---------------------------------------------------------------- denetim

class Rapor:
    def __init__(self):
        self.bulgular = []          # (paragraf_no, satir, duzey, kural, mesaj)
        self.sayac = Counter()

    def ekle(self, pno, satir, duzey, kural, mesaj):
        self.bulgular.append((pno, satir, duzey, kural, mesaj))
        self.sayac[(kural, duzey)] += 1


def tr_kucuk(s, dil):
    if dil == "tr":
        s = s.replace("I", "ı").replace("İ", "i")
    return s.lower()


def denetle(yol, dil, bolum=None):
    satirlar = open(yol, encoding="utf-8").read().split("\n")
    bas, son = bolum_kes(satirlar, bolum)
    bos = [r for r in tablo_oku(os.path.join(KOK, "kurallar", "bos_kaliplar.md"), "bos") if r["dil"] == dil]
    muaf = [re.compile(r["regex"], re.I) for r in tablo_oku(os.path.join(KOK, "kurallar", "bos_kaliplar.md"), "muaf") if r["dil"] == dil]
    nicelik = [r for r in tablo_oku(os.path.join(KOK, "kurallar", "bos_kaliplar.md"), "nicelik") if r["dil"] == dil]
    kacin, sozluk_muaf = sozluk_oku(os.path.join(KOK, "sozluk.csv"))
    bos_rx = [(re.compile(r["regex"], re.I), r["kalıp"], r["yerine"]) for r in bos]
    nic_rx = [(re.compile(r["regex"], re.I), r["yerine"], r.get("düzey", "!") or "!") for r in nicelik]
    uyari, bol = ESIK[dil]
    rapor = Rapor()
    paras = paragraflar(satirlar, bas, son)
    tum_cumle = 0
    tanimli_kisaltma = set()
    goruldu_kisaltma = set()

    for pno, (ilk, parca, baslik) in enumerate(paras, 1):
        orijinal = " ".join(parca)
        # satır eşlemesi
        ofsetler, o = [], 0
        for k, s in enumerate(parca):
            ofsetler.append((o, ilk + k))
            o += len(s) + 1
        def satir_bul(poz):
            sat = ilk
            for o_, l_ in ofsetler:
                if poz >= o_:
                    sat = l_
            return sat
        maskeli = maskele(orijinal)
        # muaf ifadeleri boş kalıp/sözlük aramasından gizle
        arama = maskeli
        for rx in muaf + sozluk_muaf[dil]:
            arama = rx.sub(lambda x: " " * len(x.group(0)), arama)

        cumleler = cumlelere_bol(maskeli)
        tum_cumle += len(cumleler)
        uzunluklar = []
        onceki_ilk = None
        ekleme_sayisi = 0

        for ci, (cb, cs) in enumerate(cumleler):
            c_mask = maskeli[cb:cs]
            c_ara = arama[cb:cs]
            c_orj = orijinal[cb:cs]
            sat = satir_bul(cb)
            n = sozcuk_say(c_mask)
            uzunluklar.append(n)
            atif_var = bool(ATIF_RX.search(c_orj))
            sayi_var = bool(re.search(r"\d", c_mask) or YER_TUTUCU_RX.search(c_orj))
            c_kucuk = tr_kucuk(c_mask.strip(), dil)

            # 1 uzun cümle
            if n > bol:
                rapor.ekle(pno, sat, "!", "uzun-cümle", "%d sözcük (> %d): tek iddiaya göre böl" % (n, bol))
            elif n > uyari:
                rapor.ekle(pno, sat, "?", "uzun-cümle", "%d sözcük (> %d): bölünebilir mi?" % (n, uyari))

            # 2 aynı sözcükle başlama
            ilk_soz = (SOZCUK_RX.findall(c_kucuk) or [""])[0]
            if ilk_soz and ilk_soz == onceki_ilk and ilk_soz not in ISLEV_SOZCUK[dil]:
                rapor.ekle(pno, sat, "?", "aynı-başlangıç", "art arda iki cümle '%s' ile başlıyor: öznelerden birini değiştir" % ilk_soz)
            onceki_ilk = ilk_soz

            # 3 ekleme bağlacı zinciri
            if re.match(EKLEME_BAGLACI[dil], c_kucuk, re.I):
                ekleme_sayisi += 1
                if ekleme_sayisi >= 2:
                    rapor.ekle(pno, sat, "!", "bağlaç-zinciri", "paragrafta %d. ekleme bağlacı: ilişkiyi söyleyen bağlaç kullan ya da sil" % ekleme_sayisi)

            # 4 retorik soru
            if c_mask.rstrip().endswith("?") and "KAYNAK?" not in c_orj:
                rapor.ekle(pno, sat, "!", "retorik-soru", "soruyu düz cümleyle kur ('We asked whether…' / '… olup olmadığını test ettik')")

            # 5 özet kapanış
            if ci == len(cumleler) - 1 and len(cumleler) >= 3 and re.match(OZET_KAPANIS[dil], c_kucuk, re.I):
                rapor.ekle(pno, sat, "?", "özet-kapanış", "paragrafı tekrar mı ediyor? yeni sonuç yoksa sil")

            # 6 nicelik
            if not sayi_var:
                for rx, yerine, duzey in nic_rx:
                    m = rx.search(c_ara)
                    if m:
                        rapor.ekle(pno, satir_bul(cb + m.start()), duzey, "sayısız-niteleme", "'%s': %s" % (m.group(0), yerine))

            # 7 kaynak bekleyen iddia
            if not atif_var:
                m = re.search(KAYNAK_ISTEYEN[dil], c_ara, re.I)
                if m:
                    rapor.ekle(pno, satir_bul(cb + m.start()), "!", "kaynak-yok", "'%s' atıf istiyor; [Yazar Yıl] ekle (yoksa [KAYNAK?])" % m.group(0).strip())
                elif ((not baslik or re.search(SAYI_KAYNAK_BOLUMU, baslik, re.I))
                      and re.search(r"\d", re.sub(r"\b(?:Fig|Figure|Table|Tablo|Şekil|Section|SM|H|S)\.?\s*\d[\w.]*|\bk\s*=\s*\d+", "", c_mask))
                      and not re.search(BIRINCI_SAHIS[dil], c_kucuk, re.I)
                      and not YER_TUTUCU_RX.search(c_orj)):
                    rapor.ekle(pno, sat, "?", "kaynak-yok", "sayı var, atıf yok: literatürdense kaynak ekle")

            # 8 nedensel dil
            m = re.search(NEDENSEL[dil], c_ara, re.I)
            if m and dil == "en" and re.search(r"\b(not|no|never|nor|without|rather than|cannot)\b[^.;]{0,50}$", c_ara[:m.start()], re.I):
                m = None
            if m and dil == "tr" and re.search(r"^\S*\s*(değil|olmaz|olmadığ|olmuyor)", c_ara[m.end():], re.I):
                m = None
            if m:
                rapor.ekle(pno, satir_bul(cb + m.start()), "?", "nedensel", "'%s': ilişkisel dil kullan (rehber §8), mekanizma kanıtı yoksa" % m.group(0))

            # 9 öncelik iddiası
            m = re.search(ONCELIK[dil], c_ara, re.I)
            if m:
                rapor.ekle(pno, satir_bul(cb + m.start()), "?", "öncelik", "'%s': öncelik iddiası; taramayla sınırla ya da sil" % m.group(0))

            # 10 üçlü sıfat
            m = re.search(SIFAT_UCLU[dil], c_ara, re.I)
            if m:
                rapor.ekle(pno, satir_bul(cb + m.start()), "?", "refleks-üçlü", "'%s': her sıfat kanıtlı mı?" % m.group(0))

            # 11 Türkçe: -mektedir ve kendi işimiz için edilgen
            if dil == "tr":
                for m in re.finditer(r"\w+m(e|a)kt(e|a)(dır|dir|ydı|ydi)\b", c_ara, re.I):
                    rapor.ekle(pno, satir_bul(cb + m.start()), "!", "mektedir", "'%s': geniş/şimdiki zaman ya da etken geçmiş kullan" % m.group(0))
                if not atif_var:
                    m = re.search(r"\w+(ıl|il|ul|ül|ın|in|un|ün|[aeıioöuü]n)((mış|miş|muş|müş)(tır|tir|tur|tür)?|dı|di|du|dü)\b", c_ara, re.I)
                    if m:
                        rapor.ekle(pno, satir_bul(cb + m.start()), "?", "edilgen", "'%s': kendi işimizse 'biz' ile etken çatı" % m.group(0))
                # isim tamlaması zinciri: 5 sözcüklük pencerede 3 tamlayan eki
                soz = list(re.finditer(r"\w+", c_ara))
                genitif = [bool(re.search(r"(n[ıiuü]n|[ıiuü]n)$", w.group(0), re.I))
                           and w.group(0).lower() not in {"için", "bütün", "gün", "dün", "ün", "kadın", "olan"}
                           for w in soz]
                for k in range(len(soz) - 4):
                    if sum(genitif[k:k + 5]) >= 3:
                        rapor.ekle(pno, satir_bul(cb + soz[k].start()), "?", "tamlama-zinciri", "'%s…': ikiden fazla tamlama üst üste; yeniden kur" % " ".join(w.group(0) for w in soz[k:k + 5]))
                        break
            else:
                # İngilizce: "of" zinciri
                if len(re.findall(r"\bof\s+the\b", c_ara, re.I)) >= 3:
                    rapor.ekle(pno, sat, "?", "of-zinciri", "üç ya da daha çok 'of the': isim zincirini yeniden kur")

        # paragraf düzeyi
        if len(uzunluklar) >= TEKDUZE_MIN_CUMLE:
            ort = statistics.mean(uzunluklar)
            cv = statistics.pstdev(uzunluklar) / ort if ort else 0
            if cv < TEKDUZE_CV:
                rapor.ekle(pno, ilk, "?", "tekdüze-uzunluk", "cümleler hep %d±%d sözcük (CV %.2f): birini kısalt ya da birleştir" % (round(ort), round(statistics.pstdev(uzunluklar)), cv))

        # boş kalıplar
        eslesmeler = [(m.start(), m.end(), m.group(0), yerine)
                      for rx, kalip, yerine in bos_rx for m in rx.finditer(arama)]
        kaplanan = []
        for b_, s_, metin_, yerine in sorted(eslesmeler, key=lambda x: x[0] - x[1]):  # uzun önce
            if any(b_ >= kb and s_ <= ks for kb, ks in kaplanan):
                continue
            kaplanan.append((b_, s_))
            rapor.ekle(pno, satir_bul(b_), "!", "boş-kalıp", "'%s' → %s" % (metin_, yerine))

        # sözlük
        for rx, dogru in kacin[dil]:
            for m in rx.finditer(arama):
                rapor.ekle(pno, satir_bul(m.start()), "?", "sözlük", "'%s' → sözlükte '%s'" % (m.group(0), dogru))

        # tireler
        for m in re.finditer(r"—|\s–\s", maskeli):
            rapor.ekle(pno, satir_bul(m.start()), "!", "uzun-tire", "virgül, noktalı virgül ya da yeni cümle")

        # yazar adı kalıbı
        for m in re.finditer(YAZAR_KALIBI[dil], maskeli):
            rapor.ekle(pno, satir_bul(m.start()), "!", "yazar-kalıbı", "'%s': yazar adını kaldır, atıf cümle sonunda" % m.group(0))

        # sayı biçimi (yer tutucu, kod ve atıflar maskelenmiş metinde)
        sayi_metin = re.sub(r"\b(Fig|Figure|Table|Tablo|Şekil|Section|Bölüm|SM|Ek|version|sürüm|v)\.?\s*[\dS][\d.]*", lambda x: " " * len(x.group(0)), maskeli, flags=re.I)
        sayi_metin = re.sub(r"\b\d+(\.\d+){2,}\b", lambda x: " " * len(x.group(0)), sayi_metin)  # sürüm numaraları 1.3.5
        if dil == "tr":
            for m in re.finditer(r"(?<![\d.,])\d+\.\d+(?![\d.,])", sayi_metin):
                tam, ond = m.group(0).split(".")
                if len(ond) == 3 and tam != "0":
                    continue  # binlik ayırıcı olabilir (78.556)
                rapor.ekle(pno, satir_bul(m.start()), "!", "sayı-biçimi", "'%s' → ondalık virgül: '%s'" % (m.group(0), m.group(0).replace(".", ",")))
            for m in re.finditer(r"\d+(,\d+)?\s?%", sayi_metin):
                rapor.ekle(pno, satir_bul(m.start()), "!", "sayı-biçimi", "'%s' → yüzde işareti önde: '%%%s'" % (m.group(0), m.group(0).rstrip("% ")))
        else:
            for m in re.finditer(r"(?<![\d.,])\d+,\d+(?![\d.,])", sayi_metin):
                tam, ond = m.group(0).split(",")
                if len(ond) == 3 and tam != "0":
                    continue  # binlik ayırıcı (78,556)
                rapor.ekle(pno, satir_bul(m.start()), "!", "sayı-biçimi", "'%s' → ondalık nokta: '%s'" % (m.group(0), m.group(0).replace(",", ".")))
            for m in re.finditer(r"%\s?\d", sayi_metin):
                rapor.ekle(pno, satir_bul(m.start()), "!", "sayı-biçimi", "yüzde işareti sayıdan sonra: 39.5%")
            # ABD yazımı
            for rx, oneri in ABD_YAZIM:
                for m in re.finditer(rx, arama, re.I):
                    if m.group(0).lower() in ABD_MUAF:
                        continue
                    # büyük harfli ve cümle başında değilse özel ad / başlık say (ör. "Grid Computing Center")
                    if m.group(0)[0].isupper() and not re.search(r"(^|[.!?:]\s+)$", arama[:m.start()]):
                        continue
                    rapor.ekle(pno, satir_bul(m.start()), "!", "ingiliz-yazımı", "'%s' → %s" % (m.group(0), oneri))

        # kısaltmalar: ilk kullanımda tanım
        # tanım: "(ABBR)" ya da bizim kullandığımız "Türkçe terim (English term, ABBR)" biçimi —
        # parantezin içindeki her kısaltma görünümlü belirteç tanımlanmış sayılır
        for m in re.finditer(r"\(([^()]*)\)", maskeli):
            for t in re.finditer(r"\b[A-Z][A-Za-z0-9-]*[A-Z][A-Za-z0-9/-]*\b", m.group(1)):
                tanimli_kisaltma.add(t.group(0))
                # "ACMG/AMP", "ROC-AUC" gibi birleşik tanımların parçaları da tanımlı sayılır
                tanimli_kisaltma.update(re.split(r"[/-]", t.group(0)))
        for m in re.finditer(r"\b([A-Z][A-Z0-9-]*[A-Z])s?\b(\s*\()?", maskeli):
            k = m.group(1)
            if k in KISALTMA_MUAF or k in goruldu_kisaltma or len(k) < 2:
                continue
            goruldu_kisaltma.add(k)
            if k not in tanimli_kisaltma and not m.group(2):
                rapor.ekle(pno, satir_bul(m.start()), "?", "kısaltma", "'%s' ilk geçişte açılmamış" % k)

        # özette atıf
        if re.search(r"abstract|özet", baslik, re.I) and ATIF_RX.search(orijinal):
            rapor.ekle(pno, ilk, "!", "özette-atıf", "özette kaynak gösterilmez (bolum_yapisi.md §2)")

    return rapor, paras, tum_cumle, (bas, son)


def yazdir(rapor, paras, tum_cumle, aralik, yol, dil, bolum, yalniz_unlem):
    print("denetle.py | %s | dil: %s%s | satır %d–%d" % (os.path.basename(yol), dil,
          (" | bölüm: %s" % bolum) if bolum else "", *aralik))
    print("düzey: ! düzelt   ? bak\n")
    gruplu = {}
    for b in sorted(rapor.bulgular, key=lambda x: (x[0], x[1])):
        if yalniz_unlem and b[2] != "!":
            continue
        gruplu.setdefault(b[0], []).append(b)
    for pno, (ilk, parca, _) in enumerate(paras, 1):
        if pno not in gruplu:
            continue
        print("¶%d (satır %d–%d)" % (pno, ilk, ilk + len(parca) - 1))
        for _, sat, duzey, kural, mesaj in gruplu[pno]:
            print("  L%-4d %s %-16s %s" % (sat, duzey, kural, mesaj))
        print()
    toplam_u = sum(v for (k, d), v in rapor.sayac.items() if d == "!")
    toplam_s = sum(v for (k, d), v in rapor.sayac.items() if d == "?")
    print("ÖZET: %d paragraf, %d cümle | ! %d   ? %d" % (len(paras), tum_cumle, toplam_u, toplam_s))
    kurallar = sorted({k for k, _ in rapor.sayac})
    for k in kurallar:
        u, s = rapor.sayac.get((k, "!"), 0), rapor.sayac.get((k, "?"), 0)
        print("  %-16s !%-3d ?%-3d" % (k, u, s))


def main():
    ap = argparse.ArgumentParser(description="Makale taslağı denetleyicisi")
    ap.add_argument("dosya")
    ap.add_argument("--dil", choices=["tr", "en"], required=True)
    ap.add_argument("--bolum", help="yalnız başlığı bu metni içeren bölümü denetle")
    ap.add_argument("--yalniz-unlem", action="store_true", help="yalnız '!' düzeyini göster")
    a = ap.parse_args()
    rapor, paras, n, aralik = denetle(a.dosya, a.dil, a.bolum)
    yazdir(rapor, paras, n, aralik, a.dosya, a.dil, a.bolum, a.yalniz_unlem)


if __name__ == "__main__":
    main()
