# Class 16 homework plots
Emily Chen (PID:A16925878)

``` r
library(ggplot2)
```

``` r
results <- read.delim("result copy.tsv")
colnames(results) <-c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                      "qstart", "qend", "sstart", "send", "evalue", "bitscore")
```

``` r
library(ggplot2)
ggplot(results, aes(pident, bitscore)) + geom_point(alpha=0.1) + geom_smooth()
```

    `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

![](hw-16-graphs.markdown_strict_files/figure-markdown_strict/unnamed-chunk-3-1.png)

``` r
library(ggplot2)

ggplot(results, aes((pident * (qend - qstart)), bitscore))+ 
  geom_point(alpha=0.1) + 
  geom_smooth()
```

    `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

![](hw-16-graphs.markdown_strict_files/figure-markdown_strict/unnamed-chunk-4-1.png)

``` r
library(ggplot2)

ggplot(results, aes(bitscore) )+ 
  geom_histogram(bins=30)
```

![](hw-16-graphs.markdown_strict_files/figure-markdown_strict/unnamed-chunk-5-1.png)

``` r
plot(results$pident  * (results$qend - results$qstart), results$bitscore)
```

![](hw-16-graphs.markdown_strict_files/figure-markdown_strict/unnamed-chunk-6-1.png)
