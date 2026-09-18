# Hesaplanmış sonuçlar

Bu klasör `roc_chapter_analysis.R` dosyası R 4.4.3 ile çalıştırılarak elde edilen çıktıları içerir. Kod yeniden çalıştırıldığında sonuç dosyaları aynı adlarla yenilenir.

- `AL_descriptive.csv`: grup bazında tanımlayıcı istatistikler.
- `AL_vus_summary.csv` ve `AL_vus_bootstrap.csv`: VUS özeti ve 2.000 bootstrap tekrarı.
- `AL_cutoff_summary.csv`: iki yöntemin eşikleri ve sınıflama ölçütleri.
- `AL_confusion_*.csv`: sınıflama tabloları.
- `AL_primary_optima_*.csv`: birincil amaçta eşdeğer optimumlar.
- `AL_cutoff_bootstrap.csv` ve `AL_cutoff_percentiles.csv`: 1.000 eşik bootstrap tekrarı ve marjinal yüzdelikler.
- `AL_optimism_summary.csv`: dengeli doğruluk için bootstrap iyimserlik düzeltmesi.
- `normal_scenarios.csv`: normal dağılım senaryolarının sayısal sonuçları.
- `Sekil_3_ROC_yuzeyi.png`: bölümdeki ROC yüzeyi.
- `R_sessionInfo.txt` ve `analysis_provenance.txt`: çalışma ortamı ve yeniden üretim ayrıntıları.

CSV dosyalarında ondalık ayırıcı nokta, alan ayırıcı virgüldür. Bu klasör ham birey kayıtları içermez.
