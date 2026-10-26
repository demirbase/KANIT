---
name: makale-yazimi
description: KANIT makalesi (bakteriyel genomlarda AMR, ML ile aday gösterilen biyobelirteçlerin kanıt derecesine göre sınıflanması) ve ondan türeyen tez bölümleri için yazım düzeni. Bir makale bölümünü (Abstract, Introduction, Methods, Results, Discussion, Impact statement) Türkçe taslak olarak yazarken, Türkçe taslağı gözden geçirirken ya da yazarın Türkçe taslaktan kurduğu İngilizce metni denetlerken MUTLAKA kullan. "Türkçe taslak yaz", "bu paragrafı düzelt", "çevirimi kontrol et", "denetle", "draft section", "check my translation" gibi isteklerde, KANIT, ESKAPEE, unitig, evidence grade, lineage-aware gibi terimler geçtiğinde de kullan. Terim sözlüğü, boş kalıp listesi, dergi bölüm kuralları ve yalnız standart kütüphaneyle çalışan bir denetleyici betik içerir.
---

# makale-yazimi

Akış: **yapay zekâ Türkçe taslak yazar → yazar İngilizceye kendisi aktarır → son metni yazar
düzeltir.** Bu skill iki yerde devreye girer: Türkçe taslağı üretirken ve yazarın İngilizcesini
denetlerken.

## Değişmez kurallar

1. **Bilimsel içerik, sayılar ve kaynaklar değişmez.** Kaynak metinde ne yazıyorsa o aktarılır.
2. **Kaynak uydurulmaz.** Eldeki kaynakçada olmayan bir atıf gerekiyorsa `[KAYNAK?]` yazılır.
3. **Sayı uydurulmaz.** Bu çalışmanın henüz hesaplanmamış sayısı `{anahtar}` yer tutucusuyla kalır;
   yer tutucunun adı değiştirilmez.
4. **Terimler `sozluk.csv`'den alınır.** Derece adları (`confirmed`, `strong_novel`, `candidate`,
   `weak`, `none`) iki dilde de ters tırnak içinde kod olarak kalır.
5. Atıf biçimi kaynaktaki gibi kalır: `[Yazar Yıl]`, `[Yazar Yıl; Yazar Yıl]`.

Kuralların amacı metnin **açık, özgül ve doğru** olması. Bir kural bu amaca ters düşüyorsa
(ör. bölünemeyecek bir tanım cümlesi) kural bırakılır ve nedeni raporda söylenir.

## Hangi dosyayı ne zaman oku

| Dosya | Ne zaman |
|---|---|
| `kurallar/turkce_taslak.md` | Türkçe taslak yazmadan ya da düzeltmeden önce, her zaman |
| `kurallar/bolum_yapisi.md` | Bir bölümü baştan yazarken (o bölümün kısmını) |
| `kurallar/bos_kaliplar.md` | Taslak yazarken üst kısmı; tablolar denetleyici içindir |
| `kurallar/ingilizce_ceviri.md` | İngilizce metni denetlerken; yazara geri bildirim verirken |
| `sozluk.csv` | Her terim kararında. `grep -i "terim" sozluk.csv` ile ara |

## İş akışı A: Türkçe taslak yazma

1. **Girdiyi topla.** Hangi bölüm; dayanak metin (mevcut İngilizce taslak, protokol, notlar);
   kullanılabilir kaynaklar ve sayılar. Dayanak yoksa uydurma, kullanıcıya sor.
2. **Yaz.** `kurallar/turkce_taslak.md` ve o bölüm için `kurallar/bolum_yapisi.md`'ye göre yaz. Çıktı
   bir `.md` dosyası olsun. Başlık numaraları İngilizce taslakla aynı kalır. Her terimin ilk
   geçişinde İngilizcesini parantezde ver.
3. **Denetle.**
   ```bash
   python3 ${SKILL_DIZINI}/betikler/denetle.py taslak.md --dil tr
   ```
   `${SKILL_DIZINI}` yerine skill klasörünün yolu gelir; betik ona göre `betikler/denetle.py`
   yolundadır. Bulunamazsa `find / -name denetle.py -path "*makale-yazimi*" 2>/dev/null`.
4. **Raporu göster.** Kullanıcıya raporun özetini ve `!` satırlarını göster. Uzun raporu olduğu
   gibi yapıştırma.
5. **Düzelt.** Bütün `!` maddelerini düzelt. `?` maddelerinin her biri için karar ver: ya düzelt ya
   tek cümleyle neden bıraktığını söyle. Yeniden denetle. `!` sıfır olana ya da kalanların her birinin
   gerekçesi yazılana kadar tekrarla.
6. **Teslim et.** Son taslak; denetimin son özet satırı; açık kalan `[KAYNAK?]`, `[SAYI?]`,
   `[KARAR: …]` işaretlerinin listesi.

## İş akışı B: Türkçe taslağı gözden geçirme

Kullanıcı bir Türkçe taslak verip "gözden geçir" derse: önce `--dil tr` ile denetle, raporu göster,
sonra düzeltme öner. Önerilen değişiklik bilimsel içeriğe dokunuyorsa (sayı, iddia gücü, atıf)
uygulamadan önce sor.

## İş akışı C: İngilizce çeviriyi denetleme

Yazar Türkçe taslaktan kendi İngilizcesini kurar. Bu metin **yazarındır**: yeniden yazma.

1. `python3 ${SKILL_DIZINI}/betikler/denetle.py bolum_en.md --dil en` çalıştır. Uzun dosyada
   `--bolum "Introduction"` ile tek bölüm seç.
2. Raporu göster.
3. `kurallar/ingilizce_ceviri.md`'ye göre okuyarak denetle. Betiğin göremediği şeylere bak:
   - Türkçe taslakla anlam farkı (Türkçe taslak elindeyse cümle cümle karşılaştır),
   - sayıların taslakla ve kaynakla aynı olması,
   - zaman tutarlılığı,
   - iddia gücünün değişip değişmediği.
4. Geri bildirimi **satır numarası + sorun + öneri** biçiminde ver. Tam paragraf yeniden yazımı
   yalnız istenirse.

## Denetleyici hakkında

- Düzeyler: `!` düzelt, `?` bak. `--yalniz-unlem` yalnız `!` satırlarını gösterir.
- Kaynakça bölümleri, kod blokları, `>` ile başlayan not blokları ve tablolar denetlenmez.
- Kural verisi betikte değil, dosyalardadır. Boş kalıp ya da nicelik ifadesi eklemek için
  `kurallar/bos_kaliplar.md` tablolarını düzenle. Eş anlamlı yasaklamak için `sozluk.csv`'nin `not`
  sütununa `kaçın-tr: a|b` ya da `kaçın-en: c|d` ekle. Özel adı korumak için `muaf-en: regex` ekle.
- Sezgisel kurallar (`refleks-üçlü`, `tamlama-zinciri`, `edilgen`, sayı-atıf) bazı durumları kaçırır
  ve bazen yanlış alarm verir. Bu yüzden `?` düzeyindedir.

## Son kontrol (teslimden önce kullanıcıya hatırlat)

- Sayılar kaynağıyla karşılaştırıldı mı? (Ondalık ayıracı değişimi sayıyı sessizce bozar.)
- Bütün `[KAYNAK?]` işaretleri kapandı mı?
- Hedef derginin güncel yazar kılavuzuna bakıldı mı?
- Yapay zekâ kullanımı Teşekkür bölümünde beyan edildi mi? Microbiology Society bunu taslak yazımı
  ve çeviri için istiyor (`kurallar/bolum_yapisi.md` §9). Beyan metni yazarın kararıdır.
