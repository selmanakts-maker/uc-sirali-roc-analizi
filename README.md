# Üç sıralı tanı grubunda ROC analizi

**ROC yüzeyi, VUS ve eşik seçimi — elektronik ek**

Selman Aktaş'ın *Sağlık Bilimlerinde Biyoistatistik: Teoriden Uygulamaya* kitabı için hazırladığı bölüme eşlik eden analiz kodu ve sonuçlarıdır. Kod, bölümde açıklanan yöntemlerden yeniden oluşturulmuş; R 4.4.3 ile çalıştırılarak bölümdeki sonuçlarla karşılaştırılmıştır.

## Dosyalar

- `roc_chapter_analysis.R`: veriyi indirir, doğrular ve analizleri yürütür.
- `sonuclar/`: hesaplanmış özetler, bootstrap sonuçları, ROC yüzeyi ve çalışma ortamı bilgisi.
- `CITATION.cff`: elektronik ekin atıf bilgileri.

## Çalıştırma

1. Dosyaları bilgisayarınıza indirin; ZIP indirdiyseniz tamamını bir klasöre çıkarın.
2. R veya RStudio'da çalışma klasörünü `roc_chapter_analysis.R` dosyasının bulunduğu klasöre ayarlayın.
3. Aşağıdaki komutu çalıştırın:

```r
source("roc_chapter_analysis.R", encoding = "UTF-8")
```

R 3.6.0 veya üstü gerekir; bu sürüm R 4.4.3 / Windows ortamında doğrulanmıştır. Ek R paketi kurulmaz. İlk çalıştırmada internet bağlantısı gerekir. Analiz 2.000 VUS ve 1.000 eşik bootstrap tekrarı içerdiği için tamamlanması zaman alabilir.

Veri, [CRAN arşivindeki DiagTest3Grp 1.6 paketinden](https://cran.r-project.org/src/contrib/Archive/DiagTest3Grp/DiagTest3Grp_1.6.tar.gz) alınır. Beklenen MD5 özeti kodda sabittir. Paket kurulmaz; `AL.rda` veri nesnesi yüklenir. İndirilen dosyalar `veri_onbellegi/` klasöründe saklanır; sonraki çalıştırmalar doğrulanmış önbelleği kullanır. Ham birey kayıtları bu depoda dağıtılmaz.

Hesaplanan dosyalar `sonuclar/` altında aynı adlarla yenilenir. CSV dosyalarında ondalık ayırıcı nokta, alan ayırıcı virgüldür.

## Bölümdeki Kod 2: analizi çalıştırma ve sonuçları okuma

Önce bu kod çalıştırılır. Böylece Kod 1'in kullandığı `v` nesnesi de oluşur. Kod numaraları bölümdeki anlatım sırasını izler.

```r
source("roc_chapter_analysis.R", encoding = "UTF-8")
s <- read.csv(file.path(out_dir, "AL_cutoff_summary.csv"))
s[, c("method", "c1", "c2", "tcf1", "tcf2", "tcf3")]
read.csv(file.path(out_dir, "AL_vus_summary.csv"))
```

Tam analiz dosyasını daha önce aynı oturumda çalıştırdıysanız `source` komutunu tekrarlamanız gerekmez; sonuçları okuyan satırlardan devam edebilirsiniz.

## Bölümdeki Kod 1: VUS'nin doğrudan hesabı

`v[[1]]`, `v[[2]]` ve `v[[3]]`, CDR=0, 0,5 ve 1 gruplarındaki T=−kfront değerlerini içerir.

```r
z <- expand.grid(a = v[[1]], b = v[[2]], c = v[[3]])
w <- with(z,
  1.0 * (a < b & b < c) +
  0.5 * (a == b & b < c) +
  0.5 * (a < b & b == c) +
  (1/6) * (a == b & b == c)
)
mean(w) # 0.6565276
```

## Kontrol değerleri

| Ölçüt | Sonuç |
|---|---:|
| Tam gözlem sayısı | 109 (45, 43, 21) |
| VUS | 0,6565276 |
| VUS %95 yüzdelik güven aralığı | 0,5397908–0,7611000 |
| Uç grupların AUC'si | 0,9936508 |
| Youden eşikleri | −0,6902921 ve 0,6112543 |
| Uzaklık eşikleri | −1,6310137 ve 0,8582474 |
| Youden düzeltilmiş dengeli doğruluk | 0,6520801 |
| Uzaklık düzeltilmiş dengeli doğruluk | 0,6572869 |

## Yeniden üretim ayrıntıları

- Grup sırası CDR 0 → 0,5 → 1; T=−kfront; yalnızca tam gözlemler.
- Sınıflama: T≤c₁; c₁<T≤c₂; T>c₂.
- VUS eşitlik ağırlıkları 1/2 ve 1/6.
- Rastgele sayı tohumu 20260909; Mersenne-Twister / Inversion / Rejection.
- Önce 2.000 VUS, ardından tohum yenilenmeden 1.000 eşik bootstrap tekrarı.
- Yüzdelikler R `type=7` kuralıyla hesaplanır.
- Birincil amaçta 10⁻¹² tolerans; ardından en yüksek minimum TCF ve sıralı çiftlerin alt orta satırı.
- Kod, bölümdeki Şekil 3'ü üretir. Diğer şekillerin dayandığı sayısal özetleri hesaplar.

Eşik yüzdelikleri betimleyici marjinal özetlerdir; ortak güven bölgesi değildir. İyimserlik düzeltmesi iç doğrulamadır. Eğitim amaçlı bu analizdeki eşikler klinik karar sınırı önerisi değildir.

Kaynak ve yöntem kaydı `sonuclar/analysis_provenance.txt`, doğrulama ortamı `sonuclar/R_sessionInfo.txt` dosyasındadır. Veri ve paket kaynağı: Luo, J., & Xiong, C. (2012). *DiagTest3Grp: An R package for analyzing diagnostic tests with three ordinal groups*. Journal of Statistical Software, 51(3), 1–24. https://doi.org/10.18637/jss.v051.i03
