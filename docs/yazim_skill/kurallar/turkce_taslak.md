# Türkçe taslak kuralları

Kaynak: `YAZIM_REHBERI.md` §2 (cümle ve paragraf), §8 (iddia disiplini), §9 (sayı ve birim).
Rehberdeki akış tez içindi: İngilizce taslak yazılıp Türkçeye çevriliyordu. Makalede akış ters.
Türkçe taslak yazılır, yazar onu İngilizceye kendisi aktarır. Bu yüzden tek tasarım ilkesi şu:

> **Türkçe taslak, İngilizce düşünülerek yazılır.** Her cümlenin sorusu: *"Bunu İngilizceye
> aktarırken bölmem ya da yeniden kurmam gerekir mi?"* Cevap evetse cümle şimdiden bölünür.

Amaç güzel Türkçe değil. Amaç anlamı tek ve açık olan, cümle cümle aktarılabilen Türkçe.

---

## 1. Cümle

**Bir cümle, bir iddia.** İkinci bir iddia varsa nokta koy, yeni cümle kur.

| | Sözcük |
|---|---|
| Hedef | 12–22 |
| Uyarı | 25'in üstü |
| Bölünmeli | 32'nin üstü |

Türkçe sözcükler eklerle uzadığı için 25 Türkçe sözcük kabaca 30–35 İngilizce sözcüğe denk gelir.
Köşeli parantezli atıflar ve `{anahtar}` yer tutucuları (bir sözcük sayılır) sayıma girmez.

**En çok bir yan cümle.** "-dığı için … -ınca … -arak" gibi iç içe yan cümleler İngilizcede üç
cümleye dağılır. Baştan üç cümle yaz.

Kötü (36 sözcük, üç iddia):
> *Teslim edilen çalışma unitig gösterimini kullandığı ve graf, çağıranın varsayılan değeri olan
> k = 31 ile kurulduğu hâlde 03u herhangi bir --kmer argümanı geçmediği için, her unitig satırına
> grafın kurulmadığı bir k değeri yazıldı.*

İyi (üç cümle, her biri tek iddia):
> *Teslim edilen çalışma unitig gösterimini kullanıyor. Graf k = 31 ile kuruluyor; bu, çağıranın
> varsayılan değeri ve 03u onu değiştirmiyor. Bu yüzden her unitig satırına yanlış k yazıldı.*

**Cümle uzunluklarını çeşitlendir.** Bir paragraftaki bütün cümleler aynı boydaysa metin tekdüze
okunur. Kısa bir cümle bir sonucu vurgular; uzun cümle bağlamı taşır.

## 2. Fiil ve çatı

| Ne anlatılıyor | Biçim | Örnek |
|---|---|---|
| Kendi yaptığımız iş | Birinci çoğul, etken, geçmiş zaman | *eğittik, hesapladık, bulduk, test ettik* |
| Kendi bulgumuz | Birinci çoğul ya da bulgu özne | *… olduğunu bulduk.* / *Modellerin {n} tanesi … ile ilişkiliydi.* |
| Literatürdeki bulgu | Geçmiş zaman, edilgen kabul, **atıf zorunlu** | *… olduğu bildirilmiştir [Hu 2024].* |
| Genel geçer bilgi | Geniş zaman | *… ile ilişkilidir [Kaynak].* |
| Şekil/tablo içeriği | Şimdiki ya da geniş zaman | *Şekil 2 … gösteriyor.* |

**"-mektedir/-maktadır" zinciri yok.** İngilizcede karşılığı genelde basit geniş zamandır ve
Türkçede ağırlık yapar: "göstermektedir" yerine "gösteriyor" ya da "gösterir".

**Kendi işimiz için edilgen yok.** "Modeller eğitilmiştir" İngilizceye "Models were trained"
diye geçer ve kimin yaptığı kaybolur. "Modelleri eğittik" yaz.

**Boş fiil yok.** "gerçekleştirmek", "yapmak", "sağlamak" asıl eylemi gizler: "analiz
gerçekleştirilmiştir" → "analiz ettik"; "karşılaştırma yapıldı" → "karşılaştırdık".

## 3. İsim tamlaması zinciri

Üst üste ikiden fazla tamlama yapma. *"modelin özelliklerinin önem sıralamasının kararlılığı"*
İngilizcede "the stability of the importance ranking of the features of the model" olur.
Yeniden kur: *"Özelliklerin önem sıralaması, alt örneklemler arasında ne kadar kararlı?"* ya da
*"Önem sıralamasının kararlılığını ölçtük."*

## 4. Paragraf

- Bir paragraf, **bir iş**: bir bağlam, bir yöntem adımı, bir bulgu ya da bir karşılaştırma.
- 4–8 cümle. İlk cümle paragrafın ne yapacağını söyler.
- **Özet kapanış cümlesi yok.** Son cümle ya yeni bir sonuç bağlar ya da sonraki paragrafa geçer.
- Bağlaçlar ilişkiyi söyler: karşıtlık *ancak, buna karşılık*; sonuç *bu nedenle, bu yüzden*.
  "Ayrıca / Bunun yanı sıra / Bu bağlamda" zinciri yok (`bos_kaliplar.md`).
- Uzun tire (—) ve boşluklu orta tire ( – ) yok. Virgül, noktalı virgül ya da yeni cümle kullan.
  Sayı aralığı (0,7–79,3) ve birleşik ifade (tür–ilaç) muaf.

## 5. İddia disiplini (rehber §8)

**Tüm KB iddiaları ilişkiseldir.**

| Kullan | Kullanma |
|---|---|
| … ile **ilişkilidir** | **neden olur** |
| … için **öngörücüdür** | **direnç kazandırır** |
| … ile **birlikte görülür** | **sorumludur** |
| … ile **tutarlıdır** | **kanıtlar** |

Tek istisna: ortogonal, mekanizma-farkında bağımsız kanıtla doğrulanmış determinantlar (ör. bilinen
bir nokta mutasyonunun mekanizması, atıfla).

**Üç kalıp, ezberlenecek:**
- Az gözlemden (birkaç model, birkaç organizma) çıkan ilişki → **"eğilim"**; asla "tahmin" ya da
  "kestirim".
- Geçmeyen bir katman → **"eşiği geçmedi"** ve eşiği ver (ör. *"q < 0,05 eşiğini geçmedi"*); asla
  "önemsiz" ya da "etkisi yok".
- Taramada bulunamayan bir şey → **"bu taramada bulunmadı"**; asla "yoktur".

**Çekince kanıta göre ayarlanır.** Yapay zekâ metninin en güçlü izi sözcükler değil, **tekdüze
temkinliliktir**: kanıtın gücü ne olursa olsun her cümle aynı tonda çekinir ("may contribute to",
"could potentially play a role"). İnsan yazısında çekincenin derecesi kanıta göre değişir. Üç
düzeyi ayır ve her cümlede hangisinde olduğunu bil:

| Kanıt | Dil |
|---|---|
| Ölçtük, sayısı var | **düz**: "68 modelde en az bir kararlı örüntü var" |
| Tutarlı ama az gözlem | **sınırlı**: "bu altı organizmada eğilim şu yönde" |
| Gerçekten belirsiz | **çekinceli**: "bunu açıklayabilecek bir neden, ama sınamadık" |

Ölçtüğün bir şeyi çekinceyle yazmak, belirsiz bir şeyi düz yazmak kadar yanlıştır. Çekince nezaket
değil, kanıtın büyüklüğünün bildirimidir.

**Öncelik iddiası yok.** "ilk kez", "ilk defa", "yeni bir yöntem" yazma. Boşluk iddiası
gerekiyorsa taramayla sınırla: *"Bildiğimiz kadarıyla … tanımlanmadı [taranan kaynaklar]."*

## 6. Kaynaklar

- Literatürden gelen her iddia aynı cümlede köşeli parantezli atıf taşır. Biçim taslaktaki
  gibidir: `[Yazar Yıl]`, birden çoksa `[Yazar Yıl; Yazar Yıl]`. **Biçim değiştirilmez.**
- "X ve arkadaşları", "X et al." kalıbı yok (rehber §2.3). *"Bir çalışmada … gösterildi [Hu 2024]."*
- **Kaynak uydurma yok.** Gereken atıf eldeki kaynakçada yoksa cümleye `[KAYNAK?]` koy ve raporla.
- "Çalışmalar göstermiştir", "literatürde", "yaygın olarak kabul edilir" gibi her ifade atıf ister.

### Künye doğrulama — zorunlu

**Künye hiçbir zaman çıkarsanmaz.** Ne DOI sonekinden yıl, ne dergi adından cilt, ne arama
sonucunun özetinden sayfa. Her yeni künye Crossref'ten alınır:

```bash
curl -s "https://api.crossref.org/works/<DOI>" -H "User-Agent: <proje>/1.0 (mailto:<e-posta>)"
```

DOI bilinmiyorsa başlıkla aranır (`?query.bibliographic=<başlık>&rows=3`), **DOI tahmin edilmez** —
tahmin edilen bir CSBJ DOI'si bambaşka bir makale getirdi (2026-10-26).

İki tuzak:
- **Ön baskı ile yayımlanmış sürüm karışır.** Elindeki kayıt bioRxiv olabilir ama makale dergide
  çıkmıştır; Crossref'te başlık biraz değişmiş olarak durur. Yayımlanmış sürüm tercih edilir.
  (CABBAGE: elde bioRxiv 2025 kaydı vardı, gerçeği *Nature Microbiology* 10(9):2338–2352, 2025.)
- **Aynı soyadı başka makale getirir.** Dönen başlığın beklediğinle aynı olduğunu gör.

Kaynakçada zaten var sanılan bir künyeyi eklemeden önce ara: aynı anahtarın iki kez girmesi,
eksik olmasından daha kötüdür.

## 7. Sayı ve birim (rehber §9)

| | Türkçe taslak | (İngilizcesi) |
|---|---|---|
| Ondalık | **0,842** | 0.842 |
| Binlik | **78.556** | 78,556 |
| p-değeri | p = 5,7 × 10⁻¹⁴ | p = 5.7 × 10⁻¹⁴ |
| Yüzde | **%39,5** (işaret önde) | 39.5% |
| Aralık | 0,7–79,3 | 0.7–79.3 |

- **Her niteleme sayıya dayanır.** "önemli ölçüde arttı" değil, "0,79'dan 0,91'e çıktı".
- **Sayı uydurma yok.** Bu çalışmanın sayısı henüz yoksa `{anahtar}` yer tutucusu kullan
  (taslaktaki gibi `{panel.n_included}`). Yer tutucunun adı değiştirilmez.
- Ondalık basamak sayısı kaynaktakiyle aynıdır. Literatür sayıları kaynakta yazdığı gibi aktarılır.
- Birimler SI; baz çifti **bp**.
- Yüzde verirken mümkünse payı ve paydayı da ver: *"{n} modelin {k}'i (%{p})"*.

## 8. Terimler ve adlar

- Terimler `sozluk.csv`'den alınır. Bir terim **ilk geçtiğinde** İngilizcesi parantezde verilir:
  *"kanıt derecesi (evidence grade)"*. Bu hem okuru hem çeviriyi sabitler.
- Kısaltma ilk geçişte açılır: *"antimikrobiyal direnç (AMR)"*. Makale İngilizce olacağı için
  İngilizce kısaltma kullanılır (AMR, WGS, MIC); sözlükteki Türkçe kısaltmalar (TGD, MİK) yalnız tez
  içindir.
- Derece adları kod olarak kalır ve ters tırnak içinde yazılır: `confirmed`, `strong_novel`,
  `candidate`, `weak`, `none`.
- Tür adları italik, ilk geçişte tam (*Escherichia coli*), sonra kısa (*E. coli*). Aynı harfle
  başlayan iki cins aynı paragraftaysa (*Enterococcus*, *Enterobacter*, *Escherichia*) her seferinde
  tam yazılır. Gen adları italik (*mecA*), protein adları düz (PBP2a).
- Çapraz atıf etiketleri: **Şekil 2**, **Tablo 1**, **Ek Tablo S3**, **Ek Yöntem SM2**.

## 9. Biçim

- Markdown. Paragraflar arasında bir boş satır. Başlıklar İngilizce taslaktaki numaralarla aynı.
- Her paragraf kendi başına aktarılabilir olsun: önceki paragrafa "bu" diye gönderme yapıyorsan
  neye gönderme yaptığını adıyla söyle.
- Taslakta açık kalan her şey köşeli parantezle işaretlenir: `[KAYNAK?]`, `[SAYI?]`, `[KARAR: …]`.
