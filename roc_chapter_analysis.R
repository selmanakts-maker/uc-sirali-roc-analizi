# Uc sirali tani grubunda ROC analizi - yeniden uretilebilir elektronik ek
# Ilgili kitap bolumunun yazari: Selman Aktas.
# Bu ek, bolumde aciklanan yontemlere gore yeniden kurulmustur.
# Yalnizca R ile gelen paketler kullanilir. Ilk calistirmada internet gerekir.
# R >= 3.6.0. Dogrulanan ortam: R 4.4.3, Windows.
# Calistirma: source("roc_chapter_analysis.R", encoding = "UTF-8")

if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))

script_frames <- sys.frames()
source_files <- lapply(script_frames, function(x) x$ofile)
source_files <- Filter(Negate(is.null), source_files)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(source_files)) tail(source_files, 1)[[1]] else
  if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "roc_chapter_analysis.R"
root_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
out_dir <- file.path(root_dir, "sonuclar")
cache_dir <- file.path(root_dir, "veri_onbellegi")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(x, name) write.csv(x, file.path(out_dir, name), row.names = FALSE,
                                      fileEncoding = "UTF-8", na = "")

archive_url <- "https://cran.r-project.org/src/contrib/Archive/DiagTest3Grp/DiagTest3Grp_1.6.tar.gz"
archive_file <- file.path(cache_dir, "DiagTest3Grp_1.6.tar.gz")
expected_md5 <- "416c690f42846c2e3c4e1797c4b04c03"
if (!file.exists(archive_file)) {
  utils::download.file(archive_url, archive_file, mode = "wb", quiet = FALSE)
}
if (!identical(unname(tools::md5sum(archive_file)), expected_md5)) {
  stop("Paket butunluk kontrolu basarisiz. Farkli dosyayla analize devam edilmedi.")
}
utils::untar(archive_file, files = "DiagTest3Grp/data/AL.rda", exdir = cache_dir)
data_env <- new.env(parent = emptyenv())
load(file.path(cache_dir, "DiagTest3Grp/data/AL.rda"), envir = data_env)
AL <- data_env$AL
group_levels <- c("D-", "D0", "D+")
cdr_labels <- c("0", "0.5", "1")
stopifnot(nrow(AL) == 118, all(c("group", "kfront") %in% names(AL)))
complete <- !is.na(AL$group) & !is.na(AL$kfront)
v <- lapply(group_levels, function(g) -AL$kfront[complete & AL$group == g])
names(v) <- cdr_labels
stopifnot(identical(lengths(v), setNames(c(45L, 43L, 21L), cdr_labels)))

vus_direct <- function(v) {
  z <- expand.grid(a = v[[1]], b = v[[2]], c = v[[3]])
  with(z, mean(1 * (a < b & b < c) +
                 0.5 * (a == b & b < c) +
                 0.5 * (a < b & b == c) +
                 (1 / 6) * (a == b & b == c)))
}
vus_fast <- function(v) {
  a <- v[[1]]; b <- v[[2]]; c <- v[[3]]
  less <- colMeans(outer(a, b, "<"))
  equal1 <- colMeans(outer(a, b, "=="))
  greater <- colMeans(outer(c, b, ">"))
  equal3 <- colMeans(outer(c, b, "=="))
  mean(less * greater + 0.5 * equal1 * greater +
         0.5 * less * equal3 + (1 / 6) * equal1 * equal3)
}
auc_pair <- function(a, b) mean(outer(a, b, "<") + 0.5 * outer(a, b, "=="))
cf_matrix <- function(v, c1, c2) {
  do.call(rbind, lapply(v, function(x) c(sum(x <= c1), sum(x > c1 & x <= c2), sum(x > c2))))
}
evaluate <- function(v, c1, c2) {
  mat <- cf_matrix(v, c1, c2)
  tcf <- unname(diag(mat) / rowSums(mat))
  c(tcf1 = tcf[1], tcf2 = tcf[2], tcf3 = tcf[3],
    BA = mean(tcf), accuracy = sum(diag(mat)) / sum(mat),
    J3 = (sum(tcf) - 1) / 2, distance = sqrt(sum((1 - tcf)^2)))
}

# Her farkli gozlem araligi bir siniflama sinirini temsil eder.
# Ic araliklarda orta noktalar; dis araliklarda en yakin gozlemden 1 birim
# uzakliktaki temsilci kullanilir. Ayni araliktaki iki esik icin o araligin
# 1/3 ve 2/3 noktalari secilir. Dis araliklar bu durumda 2 birimle sinirlanir.
candidate_grid <- function(v) {
  u <- sort(unique(unlist(v, use.names = FALSE)))
  m <- length(u)
  lower <- c(u[1] - 2, u)
  upper <- c(u, u[m] + 2)
  representative <- (lower + upper) / 2
  idx <- expand.grid(i = seq_along(representative), j = seq_along(representative))
  idx <- idx[idx$i <= idx$j, ]
  c1 <- representative[idx$i]; c2 <- representative[idx$j]
  same <- idx$i == idx$j
  c1[same] <- lower[idx$i[same]] + (upper[idx$i[same]] - lower[idx$i[same]]) / 3
  c2[same] <- lower[idx$i[same]] + 2 * (upper[idx$i[same]] - lower[idx$i[same]]) / 3
  f <- lapply(v, function(x) findInterval(representative, sort(x)) / length(x))
  t1 <- f[[1]][idx$i]
  t2 <- f[[2]][idx$j] - f[[2]][idx$i]
  t3 <- 1 - f[[3]][idx$j]
  data.frame(c1, c2, tcf1 = t1, tcf2 = t2, tcf3 = t3,
             J3 = (t1 + t2 + t3 - 1) / 2,
             distance = sqrt((1 - t1)^2 + (1 - t2)^2 + (1 - t3)^2),
             min_tcf = pmin(t1, t2, t3), empty_middle = same)
}
select_pair <- function(grid, method, tolerance = 1e-12) {
  score <- if (method == "Youden") grid$J3 else -grid$distance
  candidates <- which(abs(score - max(score)) <= tolerance)
  n_optima <- length(candidates)
  mm <- grid$min_tcf[candidates]
  candidates <- candidates[abs(mm - max(mm)) <= tolerance]
  candidates <- candidates[order(grid$c1[candidates], grid$c2[candidates])]
  selected <- candidates[ceiling(length(candidates) / 2)]
  list(row = grid[selected, ], n_optima = n_optima,
       all_primary = grid[which(abs(score - max(score)) <= tolerance), ])
}

# Bag cozumunun ve hizli hesaplamanin dogrulanmasi.
toy <- list(c(1, 4), c(2, 5), c(3, 6))
checks <- list(toy, list(c(1, 1), c(1, 1), c(1, 1)),
               list(3, 2, 1), list(c(1, 2), c(2, 3), c(3, 4)), v)
stopifnot(all(vapply(checks, function(x) abs(vus_direct(x) - vus_fast(x)) < 1e-12, logical(1))),
          abs(vus_fast(toy) - 0.5) < 1e-12,
          abs(vus_fast(list(1, 1, 1)) - 1/6) < 1e-12,
          vus_fast(list(3, 2, 1)) == 0)

descriptive <- do.call(rbind, lapply(seq_along(v), function(g) {
  raw <- -v[[g]]
  n_all <- sum(AL$group == group_levels[g])
  q <- quantile(raw, c(.25, .5, .75), type = 7)
  data.frame(CDR = cdr_labels[g], N = n_all, n = length(raw), missing = n_all - length(raw),
             mean = mean(raw), sd = sd(raw), Q1 = q[1], median = q[2], Q3 = q[3])
}))
save_csv(descriptive, "AL_descriptive.csv")
vus <- vus_fast(v)
auc13 <- auc_pair(v[[1]], v[[3]])
grid <- candidate_grid(v)
stopifnot(nrow(grid) == 5253, sum(grid$empty_middle) == 102,
          all(grid$c1 < grid$c2), all(grid$tcf2 >= 0))
# Butun kurallarin gercek siniflamayla tutarliligini kontrol et.
all_tcf <- t(vapply(seq_len(nrow(grid)), function(i)
  evaluate(v, grid$c1[i], grid$c2[i])[1:3], numeric(3)))
stopifnot(max(abs(all_tcf - as.matrix(grid[, c("tcf1", "tcf2", "tcf3")]))) < 1e-12)
methods <- c("Youden", "Distance")
fitted <- setNames(lapply(methods, function(method) select_pair(grid, method)), methods)
cutoff_summary <- do.call(rbind, lapply(methods, function(method) {
  p <- fitted[[method]]$row
  data.frame(method, c1 = p$c1, c2 = p$c2,
             as.list(evaluate(v, p$c1, p$c2)), n_optima = fitted[[method]]$n_optima)
}))
save_csv(cutoff_summary, "AL_cutoff_summary.csv")
for (method in methods) {
  p <- fitted[[method]]$row
  mat <- cf_matrix(v, p$c1, p$c2)
  save_csv(data.frame(reference_CDR = cdr_labels, pred_0 = mat[,1], pred_05 = mat[,2], pred_1 = mat[,3]),
           paste0("AL_confusion_", method, ".csv"))
  save_csv(fitted[[method]]$all_primary, paste0("AL_primary_optima_", method, ".csv"))
}

# R'nin rastgele ornekleme algoritmasi surumler arasinda sabitlenir.
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(20260909)
resample_groups <- function(v) lapply(v, function(x) x[sample.int(length(x), length(x), replace = TRUE)])
B_vus <- 2000L
boot_vus <- replicate(B_vus, vus_fast(resample_groups(v)))
ci_vus <- quantile(boot_vus, c(.025, .975), type = 7)
vus_summary <- data.frame(VUS = vus, lower95 = ci_vus[1], upper95 = ci_vus[2],
                          AUC13 = auc13, n = sum(lengths(v)), triplets = prod(lengths(v)),
                          bootstrap_B = B_vus, seed = 20260909)
save_csv(vus_summary, "AL_vus_summary.csv")
save_csv(data.frame(replicate = seq_len(B_vus), VUS = boot_vus), "AL_vus_bootstrap.csv")

# RNG sifirlanmaz: esik bootstrap'i VUS bootstrap'inin ardindan devam eder.
B_cut <- 1000L
boot_cut <- vector("list", B_cut * length(methods))
k <- 0L
for (b in seq_len(B_cut)) {
  vb <- resample_groups(v)
  gb <- candidate_grid(vb)
  for (method in methods) {
    k <- k + 1L
    p <- select_pair(gb, method)
    train <- evaluate(vb, p$row$c1, p$row$c2)
    original <- evaluate(v, p$row$c1, p$row$c2)
    boot_cut[[k]] <- data.frame(replicate = b, method, c1 = p$row$c1, c2 = p$row$c2,
                                n_optima = p$n_optima, empty_middle = p$row$empty_middle,
                                BA_boot = unname(train["BA"]), BA_original = unname(original["BA"]),
                                optimism = unname(train["BA"] - original["BA"]))
  }
}
boot_cut <- do.call(rbind, boot_cut)
save_csv(boot_cut, "AL_cutoff_bootstrap.csv")
percentiles <- do.call(rbind, lapply(methods, function(method) {
  x <- boot_cut[boot_cut$method == method, ]
  do.call(rbind, lapply(c("c1", "c2"), function(cut) {
    q <- quantile(x[[cut]], c(.025, .975), type = 7)
    data.frame(method, cutoff = cut, lower025 = q[1], upper975 = q[2])
  }))
}))
save_csv(percentiles, "AL_cutoff_percentiles.csv")
optimism <- do.call(rbind, lapply(methods, function(method) {
  x <- boot_cut[boot_cut$method == method, ]
  apparent <- cutoff_summary$BA[cutoff_summary$method == method]
  data.frame(method, apparent_BA = apparent, mean_optimism = mean(x$optimism),
             corrected_BA = apparent - mean(x$optimism),
             multiple_optima_percent = 100 * mean(x$n_optima > 1),
             empty_middle_count = sum(x$empty_middle), bootstrap_B = B_cut)
}))
save_csv(optimism, "AL_optimism_summary.csv")

# Normal dagilim senaryolari: integrasyon ve kisitli ortak izgara.
normal_results <- do.call(rbind, lapply(c(1, 2), function(sd2) {
  cuts <- seq(-6, 6, by = .01)
  ij <- which(outer(seq_along(cuts), seq_along(cuts), "<"), arr.ind = TRUE)
  c1 <- cuts[ij[,1]]; c2 <- cuts[ij[,2]]
  t1 <- pnorm(c1, -2, 1)
  t2 <- pnorm(c2, 0, sd2) - pnorm(c1, 0, sd2)
  t3 <- 1 - pnorm(c2, 2, 1)
  normal_vus <- integrate(function(t) pnorm(t, -2, 1) * (1 - pnorm(t, 2, 1)) * dnorm(t, 0, sd2),
                          -Inf, Inf, rel.tol = 1e-10)$value
  sel <- c(which.max(t1 + t2 + t3), which.min((1-t1)^2 + (1-t2)^2 + (1-t3)^2))
  data.frame(sd2, method = methods, c1 = c1[sel], c2 = c2[sel],
             tcf1 = t1[sel], tcf2 = t2[sel], tcf3 = t3[sel],
             VUS = normal_vus, AUC13 = pnorm(4 / sqrt(2)))
}))
save_csv(normal_results, "normal_scenarios.csv")

# Sekil 3: ampirik esik aginin TCF uzayinda izdusu mu.
# VUS bu geometrik cizimden hesaplanmaz; yukaridaki gozlem ucluleri kullanilir.
draw_roc_surface <- function(v, fitted, output_file) {
  png(output_file, width = 2100, height = 1850, res = 220, pointsize = 18, type = "cairo")
  on.exit(dev.off())
  par(mar = c(1, 1, 2.4, 1), family = "serif", xpd = NA)
  project <- function(p) {
    p <- matrix(p, ncol = 3)
    cbind(.92 * (p[,1] - p[,2]), .36 * (p[,1] + p[,2]) + 1.08 * p[,3])
  }
  line3 <- function(a, b, ...) {
    p <- project(rbind(a, b)); segments(p[1,1], p[1,2], p[2,1], p[2,2], ...)
  }
  point3 <- function(p, ...) { q <- project(p); points(q[,1], q[,2], ...) }
  text3 <- function(p, ...) { q <- project(p); text(q[,1], q[,2], ...) }
  plot(NA, xlim = c(-1.24, 1.24), ylim = c(-.23, 2.02), asp = 1,
       axes = FALSE, xlab = "", ylab = "")
  cube <- as.matrix(expand.grid(x = c(0,1), y = c(0,1), z = c(0,1)))
  for (i in 1:8) for (j in i:8) if (sum(abs(cube[i,] - cube[j,])) == 1)
    line3(cube[i,], cube[j,], col = "#d7dce0", lwd = 1.6)
  u <- sort(unique(unlist(v, use.names = FALSE)))
  representatives <- c(u[1]-1, (head(u,-1)+tail(u,-1))/2, tail(u,1)+1)
  f <- lapply(v, function(x) findInterval(representatives, sort(x))/length(x))
  xyz <- function(i,j) c(f[[1]][i], f[[2]][j]-f[[2]][i], 1-f[[3]][j])
  n <- length(representatives)
  faces <- list()
  for (i in 1:(n-1)) for (j in i:(n-1)) {
    faces[[length(faces)+1L]] <- if (i == j)
      rbind(xyz(i,j), xyz(i,j+1), xyz(i+1,j+1)) else
      rbind(xyz(i,j), xyz(i+1,j), xyz(i+1,j+1), xyz(i,j+1))
  }
  depth <- vapply(faces, function(p) mean(p[,1]+p[,2]), numeric(1))
  for (i in order(depth, decreasing = TRUE)) {
    p <- project(faces[[i]])
    polygon(p, col = "#99bdcf", border = "#dce8ee", lwd = .3)
  }
  for (axis in 1:3) {
    e <- c(0,0,0); e[axis] <- 1
    origin <- if (axis == 3) c(0,1,0) else c(0,0,0)
    line3(origin, origin+e, col = "#202020", lwd = 2)
    for (t in c(0,.5,1)) {
      p <- origin+e*t; q <- project(p)
      points(q[1,1], q[1,2], pch = 3, cex = .45)
      if (t != 0) {
        dx <- if (axis == 1) .075 else -.075
        text(q[1,1]+dx, q[1,2]-.015, gsub("\\.", ",", as.character(t)), cex = .8)
      }
    }
  }
  text(.66, .11, expression(TCF[1]~"(CDR 0)"), srt = 22, cex = 1.05)
  text(-.66, .11, expression(TCF[2]~"(CDR 0,5)"), srt = -22, cex = 1.05)
  text(-1.13, .96, expression(TCF[3]~"(CDR 1)"), srt = 90, cex = 1.05)
  text(0, -.04, "0", cex = .8)
  point3(c(1,1,1), pch = 8, cex = 1.5, col = "#783f8e", lwd = 2)
  text3(c(1,1,1), labels = "\u0130deal nokta (1,1,1)", pos = 3, cex = 1)
  colors <- c("#222222", "#d55e00")
  for (k in seq_along(fitted)) {
    r <- fitted[[k]]$row
    point3(as.numeric(r[,c("tcf1","tcf2","tcf3")]), pch = c(21,24)[k],
           bg = colors[k], col = "white", cex = 1.65, lwd = 1.7)
  }
  legend("top", inset = c(0,-.01), legend = c("Youden", "\u00d6klid uzakl\u0131\u011f\u0131"),
         pch = c(21,24), pt.bg = colors, col = colors, pt.cex = 1.3,
         bty = "n", horiz = TRUE, cex = 1.05)
}
draw_roc_surface(v, fitted, file.path(out_dir, "Sekil_3_ROC_yuzeyi.png"))

capture.output(sessionInfo(), file = file.path(out_dir, "R_sessionInfo.txt"))
writeLines(c(paste("Source:", archive_url), paste("Archive MD5:", expected_md5),
             "Direction: T = -kfront; complete cases only; CDR order 0, 0.5, 1",
             "TCF decision: T<=c1; c1<T<=c2; T>c2",
             "RNG: Mersenne-Twister / Inversion / Rejection; seed 20260909",
             "VUS bootstrap B=2000, then cutoff bootstrap B=1000 without resetting seed",
             "Percentiles: type=7; primary objective tolerance=1e-12",
             "Tie break: maximum minimum TCF, then lower middle pair in c1,c2 sorted order",
             "Threshold percentiles are descriptive marginal summaries, not a joint confidence region.",
             "Outputs summarize the complete-case training sample and bootstrap internal validation.",
             "Raw person records are downloaded from CRAN, not redistributed in the electronic supplement."),
           file.path(out_dir, "analysis_provenance.txt"))
print(vus_summary, row.names = FALSE)
print(cutoff_summary, row.names = FALSE)
print(percentiles, row.names = FALSE)
print(optimism, row.names = FALSE)
message("Analiz tamamlandi. Sonuclar: ", out_dir)
