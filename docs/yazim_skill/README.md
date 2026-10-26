# Yazım skill'inin depo kopyası

Canlı sürüm yerel skill dizininde duruyor ve oturumlar **oradan** yüklüyor. Burası yalnız sürüm
kontrolündeki yedek: o dizin git dışında, yani yanlış bir düzenleme ya da silinme geri alınamaz.

Düzenleme yönü **canlıdan depoya**, tersi değil:

```bash
cp -R "$SKILL_DIZINI/makale-yazimi/." docs/yazim_skill/
```

Kopyaladıktan sonra **`SKILL.md`'yi temizle.** Canlı dosya, yükleme sırasında yerine konan bir
araç değişkeni taşıyor; depo herkese açık olduğu için araç ve sağlayıcı adları burada durmamalı.
Değişkeni `${SKILL_DIZINI}` ile değiştir, sonra kopyayı araç adlarına karşı bir kez tara —
commit kancası (`.git/hooks/commit-msg`) hangi ibarelerin yasak olduğunu zaten listeliyor.

Depodaki kopyayı düzenlemek hiçbir şeyi değiştirmez; önce canlı sürümü düzelt, sonra kopyala.
Bu dosya o yüzden var — iki yönlü doğruluk kaynağı olduğunda hangisinin hangisinden türediğini
unutmak pahalıya patlıyor (`scripts/_kaynaklar_uret.py` bir koşuşta 246 düzeltilmiş künyeyi silmişti).

## İçerik

- `SKILL.md` — iş akışları: Türkçe taslak yazma, gözden geçirme, yazarın İngilizcesini denetleme
- `kurallar/turkce_taslak.md` — cümle, çatı, iddia disiplini, **çekince kuralı**, künye doğrulama
- `kurallar/bos_kaliplar.md` — denetleyicinin okuduğu kalıp tabloları
- `kurallar/bolum_yapisi.md` — dergi bölüm kuralları
- `kurallar/ingilizce_ceviri.md` — çeviri denetim ölçütleri
- `sozluk.csv` — terim sözlüğü
- `betikler/denetle.py` — denetleyici (yalnız standart kütüphane)

## Kullanım

```bash
python3 docs/yazim_skill/betikler/denetle.py <dosya> --dil tr --bolum "<bölüm>"
```
