# Boş ve şişiren kalıplar

Bu dosya iki iş görür. İnsan için: hangi kalıbın neden zayıf olduğunu ve yerine ne geleceğini söyler.
Denetleyici (`betikler/denetle.py`) için: aşağıdaki tablolardaki `regex` sütununu okur.
Bir kalıp eklemek ya da çıkarmak için yalnız tabloyu düzenle; betiğe dokunmana gerek yok.

**Ortak ilke.** Bir sözcük ya da kalıp şu üç durumdan birindeyse zayıftır:
(1) cümleden çıkarıldığında anlam değişmiyor; (2) bir büyüklük iddia ediyor ama sayı vermiyor;
(3) okura ne düşüneceğini söylüyor ama nedenini göstermiyor. Yerine konacak şey hemen her zaman
ya **sayı**, ya **somut fiil**, ya da **hiçbir şey**dir (silmek).

**Not: bu liste dergi kuralı değil, ev kuralıdır.** Yüksek etkili dergilerde 2025'te yayımlanmış
AMR/ML makalelerinin özetlerinde bile "comprehensive", "underscore", "highlight" geçiyor. Bu sözcükler
reddedilme nedeni değil. Listeyi kullanma nedenimiz metnin özgül kalması ve yapay zekâ taslağının
kalıplarına düşmemek.

Tablo biçimi: `dil | kalıp | regex | neden zayıf | yerine`. Regex büyük-küçük harfe duyarsızdır;
tablo içinde `|` işareti `\|` olarak yazılır.

<!-- DENETIM-TABLOSU: bos -->
| dil | kalıp | regex | neden zayıf | yerine |
|---|---|---|---|---|
| en | delve | `\bdelv(e\|es\|ed\|ing)\b` | Süslü "inceledik"; yapay zekâ metninin en bilinen izi. | examine, analyse, test, compare (ne yapıldıysa o) |
| en | pivotal | `\bpivotal\b` | Önemi iddia eder, göstermez. | Etkiyi sayıyla ver ya da sil. |
| en | crucial | `\bcrucial(ly)?\b` | Aynı; okura neden önemli olduğunu söylemez. | "X is required for Y" gibi somut ilişki; ya da sil. |
| en | landscape | `\blandscapes?\b` | Mecaz; neyin kapsandığını belirsiz bırakır. | Kapsananı adlandır: "the 83 species–drug combinations", "published tools". |
| en | realm | `\brealms?\b` | Mecaz, içi boş. | field, area ya da sil. |
| en | underscore | `\bunderscor(e\|es\|ed\|ing)\b` | Bulgunun ne gösterdiğini söylemeden "vurgular". | show, indicate + ne gösterdiği. |
| en | highlight | `\bhighlight(s\|ed\|ing)?\b` | Underscore ile aynı işlev. | show, indicate, report. |
| en | leverage | `\bleverag(e\|es\|ed\|ing)\b` | İşletme dili; "kullanmak"tan fazlasını söylemez. | use |
| en | harness | `\bharness(es\|ed\|ing)?\b` | Mecaz, "kullanmak". | use |
| en | seamless | `\bseamless(ly)?\b` | Doğrulanamaz övgü. | Sil ya da somut özelliği söyle ("runs without manual steps"). |
| en | robust (istatistik dışı) | `\brobust(ly\|ness)?\b` | Neye karşı dayanıklı olduğu söylenmezse boş övgü. İstatistiksel anlamı ("robust standard errors") muaftır. | "stable under resampling (Supplementary Table S3)", "unchanged when X". |
| en | comprehensive | `\bcomprehensive(ly)?\b` | Kapsamı iddia eder, ölçmez. Özel ad (CARD) muaf. | Kapsamı sayıyla ver: "all {n} models", "six species". |
| en | notably | `\bnotably\b` | Okura neye dikkat edeceğini söyler, nedenini değil. | Sil; cümle önemliyse kendisi gösterir. |
| en | interestingly / importantly | `\b(interestingly\|importantly\|remarkably\|strikingly)\b` | Aynı; yazarın tepkisi, bulgu değil. | Sil. |
| en | plays a key role | `\bplay(s\|ed\|ing)? an? (key\|crucial\|pivotal\|central\|critical\|vital\|important\|major) role\b` | Ne yaptığını söylemez. | Mekanizmayı ya da ilişkiyi yaz: "is associated with", "is required for". |
| en | it is important to note | `\bit is (important\|worth\|interesting) (to note\|noting)\b` | Dolgu. | Sil, cümleye doğrudan başla. |
| en | paving the way | `\bpav(e\|es\|ed\|ing) the way\b` | Kanıtsız gelecek iddiası. | Neyi mümkün kıldığını somut söyle ya da sil. |
| en | shed light | `\bsh(ed\|eds\|edding) (new )?light\b` | Mecaz. | show, explain, test. |
| en | state-of-the-art / cutting-edge | `\b(state-of-the-art\|cutting-edge\|groundbreaking\|game-chang\w*)\b` | Pazarlama dili. | Karşılaştırmayı sayıyla ver. |
| en | utilise | `\butili[sz](e\|es\|ed\|ing\|ation)\b` | "use"un uzun hali. | use |
| en | in order to | `\bin order to\b` | "to" yeter. | to |
| en | due to the fact that | `\bdue to the fact that\b` | Dolambaçlı. | because |
| en | additionally | `\badditionally\b` | Dolgu bağlaç; Türkçedeki "ayrıca" zincirinin İngilizce ikizi. Kobak 2025'in fazla-sözcük listesinde üst sıralarda. | Sil; cümle kendi başına dursun. İlişki varsa onu söyle. |
| en | showcase | `\bshowcas(e\|es\|ed\|ing)\b` | Pazarlama dili; "göstermek"ten fazlasını söylemez. | show, report, present |
| en | a wide range of | `\ba (wide\|broad\|diverse) (range\|array\|variety) of\b` | Ne kadar geniş olduğunu söylemez. | Sayı ver: "{n} models", "six species". |
| tr | önem arz etmektedir | `önem arz\w*` | Önemi iddia eder; ağır ve kalıp. | Neden önemli olduğunu söyle; ya da sil. |
| tr | kritik bir rol oynamaktadır | `(kritik\|önemli\|kilit\|merkezi\|büyük) (bir )?rol oyn\w*` | "plays a key role" karşılığı; işlevi söylemez. | "… ile ilişkilidir", "… için gereklidir". |
| tr | büyük önem taşımaktadır | `(büyük\|önemli\|kritik) (bir )?önem taşı\w*` | Aynı. | Sil ya da somut gerekçe. |
| tr | ışık tutmaktadır | `ışık tut\w*` | Mecaz. | "gösterir", "açıklar", "test eder". |
| tr | gerçekleştirilmiştir | `gerçekleştir\w*` | Boş fiil; asıl eylemi gizler, edilgen çatıyı çağırır. | Somut fiil: "eğittik", "hesapladık", "taradık". |
| tr | yapılmıştır (boş fiil) | `\b(analiz\|değerlendirme\|karşılaştırma\|inceleme\|hesaplama\|tarama) (işlemi )?yapıl\w*` | "analiz yapıldı" = "analiz ettik"; isimleştirme. | "analiz ettik", "karşılaştırdık". |
| tr | söz konusu | `söz konusu` | Gönderimi bulanıklaştırır; İngilizcesi ("the said") yok. | "bu", "bu model" ya da adı tekrar et. |
| tr | bu bağlamda | `bu bağlamda` | Bağlamı söylemez; zincirlenince içi boşalır. | İlişkiyi söyleyen bağlaç: "bu nedenle", "buna karşılık" ya da sil. |
| tr | son derece | `son derece` | Büyüklük iddiası, sayısız. | Sayı ver ya da sil. |
| tr | kapsamlı bir şekilde | `kapsamlı (bir )?(şekilde\|biçimde)` | Kapsamı ölçmez. | Kapsamı sayıyla ver. |
| tr | dikkat çekici biçimde | `dikkat çek(ici\|mektedir\|en)\w*` | "Notably" karşılığı. | Sil. |
| tr | önemli bir yere sahiptir | `önemli bir (yere\|yer) (sahip\|tut)\w*` | Dolgu. | Sil ya da somut gerekçe. |
| tr | katkı sağlamaktadır | `katkı sağla\w*` | Boş fiil. | "katkıda bulunur" ya da ne yaptığını söyle. |
| tr | mevcut olan | `mevcut olan` | "olan" fazlalık. | mevcut, var olan |
<!-- /DENETIM-TABLOSU -->

## Tabloya girmeyen, ama dikkat edilecek sözcükler

Kobak ve ark. (2025, *Science Advances*, 15 milyon PubMed özeti) LLM sonrası fazla kullanılan
sözcükleri çıkardı; fazlalık isimlerden **fiil ve sıfatlara** kaydı. Listenin üst sıralarındaki
`across`, `particularly` ve `enhance` bizim tablomuzda **yok**, çünkü üçü de yerinde kullanıldığında
doğru sözcük ("across the six organisms", "particularly in *E. coli*", "overexpression enhances
efflux"). Regexle yasaklamak yanlış alarm üretir. Yazarken kendine sor: bu sözcük bir şey mi
söylüyor, yoksa cümleyi mi dolduruyor?

Aynı çalışmanın bir bulgusu daha var: sonraki yıllarda işaretler **zayıflıyor**, çünkü yazarlar
metni daha çok düzenliyor. Yani sözcük listesi tek başına yetmez; asıl iş aşağıdaki çekince
kuralında.

## Muaf tutulanlar

Bu ifadelerin içindeki sözcükler yukarıdaki tablo yüzünden işaretlenmez (özel adlar ve terim anlamları).

<!-- DENETIM-TABLOSU: muaf -->
| dil | regex | neden |
|---|---|---|
| en | `Comprehensive Antibiotic Resistance Database` | CARD'ın özel adı. |
| en | `robust (standard errors?\|regression\|estimat\w+\|statistic\w*\|covariance)` | İstatistiksel terim. |
| en | `highlighted in (Fig\|Figure\|Table)` | Şekilde görsel vurgu, mecaz değil. |
<!-- /DENETIM-TABLOSU -->

## Sayısal karşılığı olmayan nitelemeler

Cümlede sayı (rakam ya da `{anahtar}` yer tutucusu) **yoksa** işaretlenir. Sayı varsa niteleme
o sayıya dayanıyor kabul edilir. "significant" sözcüğü yalnız istatistiksel test anlamında kullanılır.
Kullanılıyorsa test ve değeri aynı cümlede ya da hemen öncesinde verilmelidir.
Bu yüzden "significant" ve "anlamlı" yalnız `?` (bak) düzeyindedir. Tartışmada bir kuralı anlatırken
("remains significant after correction") sayı olmadan da doğru kullanılabilir.

<!-- DENETIM-TABLOSU: nicelik -->
| dil | regex | yerine | düzey |
|---|---|---|---|
| en | `\bsignificant(ly)?\b` | Test + değer ver (p, q, FDR) ya da "markedly" yerine sayı. | ? |
| en | `\b(substantial(ly)?\|considerabl[ey]\|marked(ly)?\|dramatic(ally)?\|greatly\|vastly\|largely)\b` | Farkı sayıyla ver. | ! |
| en | `\bhighly (accurate\|predictive\|effective\|variable\|conserved\|similar)\b` | Değeri ver (AUC, %, kimlik). | ! |
| en | `\b(the vast majority\|a (large\|great\|substantial\|small) (number\|proportion\|fraction\|share\|majority))\b` | Oranı ver. | ! |
| en | `\b(much\|far) (higher\|lower\|better\|worse\|larger\|smaller\|more\|fewer)\b` | Farkı ver. | ! |
| tr | `önemli ölçüde` | Farkı sayıyla ver. | ! |
| tr | `anlamlı (derecede\|ölçüde\|biçimde\|şekilde\|olarak\|bir)` | Test ve değeri ver; yoksa "anlamlı" deme. | ? |
| tr | `büyük (oranda\|ölçüde\|çoğunluğu?)` | Oranı ver. | ! |
| tr | `(yüksek\|düşük) oranda` | Oranı ver. | ! |
| tr | `belirgin (biçimde\|şekilde\|olarak\|bir)` | Farkı ver. | ! |
| tr | `(kayda\|dikkate) değer` | Değeri ver. | ! |
| tr | `\boldukça\b` | Sayı ver ya da sil. | ! |
| tr | `çok sayıda\|\bbirçok\b\|\bçoğu\b` | Kaç tane olduğunu söyle. | ! |
<!-- /DENETIM-TABLOSU -->

## Yapısal kalıplar (betik bunları kendi kurallarıyla arar)

**Cümle başında zincirleme ekleme bağlaçları.** *Furthermore, Moreover, Additionally, In addition,
Also* / *Ayrıca, Bunun yanı sıra, Buna ek olarak, Ek olarak, Öte yandan, Bu bağlamda*. Bunlar
cümleler arasındaki ilişkiyi söylemez, yalnızca "bir şey daha" der. Bir paragrafta ikinci kez
geçtiğinde paragraf bir liste haline gelmiştir. Yerine ilişkiyi söyleyen bağlaç gelir (*However*,
*Therefore*, *In contrast*). Ya da hiç bağlaç gelmez: iki cümle yan yana duruyorsa ilişki zaten
"ve"dir.

**Refleks üçlüler.** *"fast, accurate and scalable"*, *"hızlı, doğru ve ölçeklenebilir"*. Üç sıfat
ritim için seçilmişse en az biri iddiayı taşımıyordur. Her sıfat ayrı kanıt istiyorsa ayrı cümle
olsun; istemiyorsa çıkar. Tür, ilaç ya da yöntem listeleri bu kurala girmez. Betik yalnız sıfat
benzeri üçlüleri "bak" düzeyinde işaretler.

**Uzun tire bolluğu.** Uzun tire (—) ve boşluklu orta tire ( – ) düz metinde kullanılmaz (rehber
§2.2). Araya giren cümle ya virgülle ya da ayrı cümleyle yazılır. Sayı aralıkları (0.7–79.3) ve
birleşik ifadeler (species–drug) muaftır.

**Paragrafı özetleyen kapanış cümleleri.** *In summary, Overall, Taken together, Collectively, Thus
these results* / *Özetle, Kısacası, Sonuç olarak, Genel olarak, Birlikte değerlendirildiğinde*.
Paragraf az önce söylediğini tekrar ediyorsa bu cümle silinir. Cümle yeni bir sonuç bağlıyorsa
kalır, ama o zaman sonucun kendisiyle başlar.

**Retorik sorular.** Bilimsel metinde soru, cevabı ertelemekten başka iş görmez. Soru bir
araştırma sorusuysa düz cümleyle yazılır: *"We asked whether…"* / *"… olup olmadığını test ettik."*

**Öncelik iddiaları.** *first, for the first time, novel* / *ilk kez, ilk defa, yeni bir*. Dergiler
kanıtlanmamış öncelik iddialarına karşı uyarıyor; hakemler genelde bunlara takılıyor. Yerine
**neyin taranıp bulunamadığı** yazılır: *"was not identified in the corpus surveyed here"*.
Kod olan `strong_novel` muaftır.

**Yazar adı kalıbı.** *X and colleagues, X et al. showed* / *X ve arkadaşları*. Rehber §2.3 uyarınca
kullanılmaz; atıf cümle sonundaki köşeli paranteze gider.
