# Bölüm yapısı (dergi kılavuzlarından)

Bu dosya, alanın dergilerinin yazarlardan ne istediğini bölüm bölüm toplar. Taslakta *Impact
statement* ve *Data summary* bölümleri olduğu için Microbial Genomics (Microbiology Society) biçimi
temel alındı. Nature dergileri ve ML raporlama önerileri (DOME) eklendi.

Kaynaklar (Ekim 2026'da bakıldı):
- Microbiology Society, *How to prepare an article for submission*:
  https://www.microbiologyresearch.org/prepare-an-article
- Nature, *Formatting guide*: https://www.nature.com/nature/authors/gta
- DOME önerileri (denetimli ML raporlama): Walsh ve ark., Nature Methods 2021; DOME Registry,
  GigaScience 2024 (PMC11633452)
- Örnek alınan güncel AMR/ML makalesi: Yu, Wheeler, Barquist, PLOS Biology 2025,
  doi:10.1371/journal.pbio.3003539

⚠️ Hedef dergi kesinleşince o derginin güncel kılavuzuna bakılır; bu dosya o kılavuzun yerine geçmez.

---

## 1. Başlık

Ana sonucu vurgulayan kısa bir ifade (Microbiology Society). Yöntemin adı değil, bulgu öne çıksın.

## 2. Özet (Abstract)

- İlk cümle konuyu tanıtır; **son cümle ana sonucu verir** (Microbiology Society).
- **Özette kaynak gösterilmez.** Kullanılan her kısaltma özetin içinde açılır.
- Alanda yerleşik akış (ör. Yu ve ark. 2025 özeti) beş adımdır:
  1. bağlam (1–2 cümle),
  2. boşluk, genelde *"However, …"* (1 cümle),
  3. *"Here we …"* ile ne yapıldığı (1 cümle),
  4. bulgular, **sayılarla** (2–4 cümle),
  5. sonuç ve anlamı (1 cümle).
- Nature'ın özet paragrafı ayrıca sayı ve kısaltmadan kaçınmayı ister. Microbial Genomics için
  bu kural yoktur, bulgular sayıyla verilir.

## 3. Impact statement ve Data summary

- **Impact statement:** makalenin önemini ve alana ne eklediğini anlatır. Uzman olmayan okura
  yazılır. Ne kadar geniş bir okur kitlesini ilgilendirdiği ve katkının adım mı yoksa artış mı
  olduğu açıkça söylenir. Abartı burada en çok göze batar: `bos_kaliplar.md` burada da sıkı uygulanır.
- **Data summary:** bütün veri, kod ve yazılım; her biri için DOI ya da erişim numarası ve URL.
  Microbial Genomics'te açık veri zorunludur. Zenodo için kavram DOI'si verilir (rehber §10.4).

## 4. Giriş

- Kısa tutulur; amacı söyler ve ilgili çalışmaları anarak zemini kurar. Okurun makalenin geri
  kalanını yorumlayabileceği kadar ayrıntı verilir, fazlası verilmez.
- Paragraf sırası: sorun → mevcut yaklaşım → boşluk → *"Here we…"* paragrafı (ne yaptık, hangi
  hipotezleri test ettik).
- Boşluk cümlesi öncelik iddiası değil, tarama sonucu olarak kurulur (bkz. `bos_kaliplar.md`,
  öncelik iddiaları).

## 5. Yöntemler

- Başkasının tekrar edebileceği kadar ayrıntı. **Her yazılım sürüm numarası ve parametreleriyle
  anılır** (Microbiology Society). Örnek: *Bifrost 1.3.5, k = 31*.
- ML kısmı DOME'un dört başlığını karşılar. Her biri bir alt başlık ya da en az bir paragraf olur:
  - **Veri:** kaynak, sayılar, dışlama ölçütleri, eğitim/test bölmesi ve bölmeler arasındaki
    bağımsızlık (bizde soy-farkında bölme ve neden gerektiği).
  - **Optimizasyon:** özellik kodlaması, hiperparametreler, nasıl seçildikleri (iç içe tasarım),
    erken durdurma.
  - **Model:** model türü, yorumlanabilirlik (özellik önemi, permütasyon önemi), modelin ve kodun
    erişilebilirliği.
  - **Değerlendirme:** metrikler (birincil: soy-farkında kolda havuzlanmış katlar dışı ROC-AUC, beş tekrarın
    ortalaması; protokol §6.3), karşılaştırma temeli
    (şans düzeyi, etiket permütasyonu), değerlendirme çıktılarının erişilebilirliği.
- Önceden belirleme (protokol) burada, tarihiyle yazılır.

## 6. Bulgular

- Yeterli sayıda alt başlık; **alt başlık bulguyu söyler** (*"Lineage-aware evaluation lowers
  AUC in five of six species"*), konuyu değil (*"Model performance"*).
- Her sonuçta: kaç kez tekrarlandığı (tohum sayısı, kat sayısı), verilen özetin türü (ortalama,
  medyan) ve değişkenliği (SD, %95 GA) yazılır.
- Kısa bir cümleyle söylenebilecek sonuç için şekil ya da tablo kullanılmaz.
- Yorum yapılabilir ama karşılaştırma ve genelleme Tartışmaya kalır.

## 7. Tartışma

- Bulgular önceki çalışmalarla **kısaca** karşılaştırılır; sonuçlar baştan anlatılmaz.
- Her paragraf bir soruya cevap verir (taslaktaki 4.2 *What a grade means*, 4.3 *Why several lines
  of evidence* bu yapıda; korunmalı).
- Sınırlılıklar somut olur: neyin ölçülmediği ve bunun sonucu nasıl etkileyebileceği.
- Tartışmada yeni sonuç verilmez.

## 8. Şekil ve tablolar

- Metne bakmadan anlaşılır; kullanılan kısaltma ve simgeler lejantta açılır.
- Yöntem ayrıntısı lejantta tekrarlanmaz.
- Panel harfleri küçük ve parantez içinde: (a), (b).
- Her şekle en çok 255 karakterlik alternatif metin (alt-text) önerilir.
- Tekrar sayısı ve değişkenlik lejantta belirtilir.

## 9. Beyanlar

- **Yapay zekâ kullanımı beyan edilir.** Microbiology Society, taslak yazımı, düzenleme, çeviri,
  argüman kurma ve üslup değişikliği için üretken yapay zekâ kullanımının Teşekkür bölümünde
  belirtilmesini istiyor. Bu iş akışı (Türkçe taslağı yapay zekâ yazıyor, yazar çeviriyor) tam
  olarak bu kapsamda. Beyan metni yazarın kararıdır; skill bunu yalnızca hatırlatır.
- Yazar katkıları CRediT biçiminde; çıkar çatışması ve fon cümleleri zorunlu.
- Kaynak listesi gönderimde tutarlı olduğu sürece herhangi bir biçimde olabilir. Kabulde Vancouver'a
  çevrilir; DOI'ler eklenmelidir.

## 10. Adlandırma

- Prokaryot adları ICNP'ye uyar. Cins + tür ilk geçişte tam yazılır, sonra cins baş harfiyle
  kısaltılır. Aynı baş harfli birkaç cins varsa her seferinde tam yazılır.
- Gen ile protein ayrılır (*gag* / Gag).
- SI birimleri.
