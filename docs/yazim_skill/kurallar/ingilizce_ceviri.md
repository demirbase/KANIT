# İngilizce çeviri kuralları (yazarın uygulayacağı)

Türkçe taslaktan İngilizce makale metnini yazar kendisi kurar. Bu dosya o adımın kurallarıdır.
Denetleyici `--dil en` ile bunların mekanik kısmını yakalar; geri kalanı okuyarak denetlenir.

## 1. Anlamı aktar, sözcüğü değil

- Türkçe cümleyi oku, **iddiasını** tek cümleyle İngilizce söyle, sonra Türkçeyle karşılaştır.
  Sözcük sözcük çeviri Türkçe söz dizimini İngilizceye taşır ("The stability of the ranking of the
  importance of the features…").
- Türkçede iki cümle olan şey İngilizcede de iki cümle kalır. Birleştirmeyi yalnız iki cümle aynı
  iddiayı taşıyorsa yap.
- Çeviriden sonra paragrafı bir kez **yalnız İngilizce** oku. Akış, Türkçeye bakmadan anlaşılmalı.

## 2. Cümle

- Hedef 15–25 sözcük; 30 üstü uyarı; 35 üstü bölünür (rehber §2.1).
- Bir cümle, bir iddia.
- **Etken çatı, "we".** Nature dergileri yazarların etken çatıyı kullanmasını açıkça tercih ediyor;
  Microbial Genomics ve benzerlerinde de "Here we…" kalıbı standart. Kendi işimiz: *We trained,
  we found, we graded.*
- **Ama her cümle "We" ile başlamasın.** Arka arkaya üç "We" tekdüzedir. Bulguyu özne yap:
  *"Twelve of the {n} models…"*, *"The grade depends on…"*.
- Yöntemlerde edilgen yalnız eylemi yapanın önemsiz olduğu yerde kalabilir: *"Genomes were
  downloaded from BV-BRC."*

## 3. Zaman

| Ne | Zaman |
|---|---|
| Yaptığımız iş, bulduğumuz sonuç (Methods, Results) | Geçmiş: *we trained, models reached* |
| Yerleşik bilgi, genel ilke | Geniş: *Many pathogens are clonal [..]* |
| Literatürdeki tek bir çalışmanın bulgusu | Geçmiş: *…fell from 64% to 25–33% [Hu 2024]* |
| Şekil/tablo ne gösteriyor | Şimdiki: *Fig. 2 shows…* |
| Makalenin kendisi (giriş sonu, tartışma) | Şimdiki: *Here we present…* |

## 4. Yazım: İngiliz İngilizcesi, tutarlı

Mevcut taslak İngiliz yazımıyla yazılmış (*optimisation, generalisability, neighbouring,
catalogue, licence*). Aynı yazım sürdürülür:

| Kullan | Kullanma |
|---|---|
| -ise, -isation (optimise, characterise) | -ize, -ization |
| -yse (analyse) | -yze |
| -our (behaviour, colour, neighbour) | -or |
| modelling, labelled, travelled | modeling, labeled |
| catalogue | catalog |
| licence (isim), license (fiil) | — |
| centre | center |

Özel adlar ve alıntılanan başlıklar olduğu gibi kalır (ör. bir yazılımın ya da kılavuzun adı).

## 5. Sayılar

| Türkçe taslak | İngilizce |
|---|---|
| 0,842 | **0.842** |
| 78.556 | **78,556** |
| %39,5 | **39.5%** |
| p = 5,7 × 10⁻¹⁴ | p = 5.7 × 10⁻¹⁴ |
| 0,7–79,3 | 0.7–79.3 |

⚠️ Ondalık ayıracı değişimi sayıyı **sessizce** bozar. Çeviriden sonra her sayıyı Türkçe
taslaktakiyle ve kaynağıyla (`results/tables/`, `tez_sayilari.csv`) karşılaştır. `{anahtar}`
yer tutucularına dokunma. Ondalık basamak sayısını değiştirme.

Cümle başındaki sayı sözcükle yazılır (*"Twelve models…"*) ya da cümle yeniden kurulur.

## 6. Çapraz atıf etiketleri

| Türkçe | İngilizce |
|---|---|
| Şekil 2 | Fig. 2 (cümle başında: Figure 2) |
| Şekil 2a | Fig. 2(a) |
| Tablo 1 | Table 1 |
| Ek Tablo S3 | Supplementary Table S3 |
| Ek Şekil S1 | Supplementary Fig. S1 |
| Ek Yöntem SM2 | Supplementary Methods SM2 |
| Ek / Ek Dosya | Supplementary File |

Çok panelli şekillerde panel harfleri küçük ve parantez içinde (a), (b) yazılır (Microbiology
Society biçimi).

## 7. Adlar

- Tür adları italik: *Klebsiella pneumoniae*, sonra *K. pneumoniae*. Aynı harfle başlayan cinsler
  (*Enterococcus*, *Enterobacter*, *Escherichia*) aynı bağlamdaysa kısaltma belirsizleşir; tam yaz.
- Cins adı tek başına üyeler için çoğul yapılmaz: *"strains of Salmonella"*, *"the Salmonellae"*
  değil.
- Gen italik (*mecA*), protein düz (PBP2a).
- Kısaltmalar her bölümde değil, metinde ilk geçişte açılır; **özet kendi içinde ayrıca açar.**
- Derece adları ters tırnak içinde kod olarak kalır: `confirmed`, `strong_novel`, `candidate`,
  `weak`, `none`.

## 8. Kaynaklar

- Atıf biçimi **değişmez**: `[Yazar Yıl]`, `[Yazar Yıl; Yazar Yıl]`. Sırası, yazımı, noktalaması
  Türkçe taslaktakiyle aynı kalır. Dergi biçimine (Vancouver vb.) dönüşüm en son, kaynak
  yöneticisiyle yapılır.
- `[KAYNAK?]` işareti çeviride de kalır; kaynak bulunmadan metin teslim edilmez.
- Yazar adı kalıbı (*X and colleagues*, *X et al. showed*) yok.

## 9. İddia dili

Türkçe taslaktaki ilişkisel dil İngilizcede de korunur (rehber §8):
*is associated with, is predictive of, co-occurs with, is consistent with.* Çeviri sırasında
"causes", "confers", "proves" yazılmaz; istisna, bağımsız mekanizma kanıtı olan determinantlar.

## 10. Bitirme kontrolü

1. `python3 betikler/denetle.py bolum.md --dil en` çalıştır; `!` düzeyindeki her maddeyi kapat.
2. Sayıları kaynakla karşılaştır.
3. Terimleri `sozluk.csv`'ye göre tara (denetleyici kaçınılacak eş anlamlıları yakalar, ama sözlükte
   olmayan bir sapmayı yakalayamaz).
4. Yapay zekâ kullanımını beyan et (bkz. `bolum_yapisi.md` §9).
