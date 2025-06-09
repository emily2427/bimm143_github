# Class 7: Machine Learning 1
Emily Chen (PID: A16925878)

-   [Clustering](#clustering)
    -   [K-means](#k-means)
    -   [Hierarchical Clustering](#hierarchical-clustering)
-   [Principle Component Analysis
    (PCA)](#principle-component-analysis-pca)
    -   [Data import](#data-import)
    -   [PCA to the rescue](#pca-to-the-rescue)

Today we will explore unsupervised machine learning methods starting
with clustering and dimensionality reductions.

## Clustering

To start, let’s make up some data to cluster where we know what the
answer should be. The `rnorm` functions will help us here.

``` r
hist(rnorm(10000, mean=3))
```

![](class07.markdown_strict_files/figure-markdown_strict/unnamed-chunk-1-1.png)

Return 30 numbers centered on -3 So we are going to do rnorm() with 30
inside, and then for the mean we put -3 so that it will be centered at
-3.

``` r
rnorm(30, mean=-3)
```

     [1] -1.892397 -3.135979 -3.564512 -2.559882 -2.090461 -4.076803 -3.151369
     [8] -3.389699 -2.842163 -2.876741 -2.812701 -4.121359 -2.361715 -3.405383
    [15] -2.584585 -2.154837 -4.383828 -1.384394 -2.707144 -3.365970 -3.415209
    [22] -4.020589 -1.144613 -3.520380 -2.971186 -3.285926 -2.359362 -1.401008
    [29] -1.751466 -1.886200

``` r
rnorm(30, mean=+3)
```

     [1] 2.026991 2.286506 3.305365 3.335182 2.818854 3.111344 3.166591 4.127781
     [9] 1.666207 3.252620 3.489044 4.576076 3.656446 3.520015 3.515030 1.848468
    [17] 2.883103 3.993563 4.719987 4.441387 2.899704 1.599343 5.294527 2.904015
    [25] 2.842897 3.035425 4.580256 4.855444 3.174126 2.775521

To combine the two, we will create a vector using the c() function.
Below we have 60 points in total

``` r
tmp<-c(rnorm(30, mean=-3),
  rnorm(30, mean=+3) )

x<- cbind(x=tmp, y=rev(tmp))

x
```

                  x         y
     [1,] -4.695592  3.081494
     [2,] -2.664620  3.305911
     [3,] -3.856196  1.683888
     [4,] -3.549683  2.271456
     [5,] -2.700320  3.387088
     [6,] -1.876462  3.468711
     [7,] -3.599468  4.495104
     [8,] -4.077324  3.734322
     [9,] -2.496261  4.612891
    [10,] -2.555539  1.347806
    [11,] -2.788996  3.330045
    [12,] -3.283969  2.383132
    [13,] -1.097882  4.885257
    [14,] -2.239879  3.447379
    [15,] -3.481751  3.404251
    [16,] -4.142593  1.596865
    [17,] -3.599175  3.196438
    [18,] -3.493054  3.008510
    [19,] -3.189386  4.606236
    [20,] -1.615082  2.877253
    [21,] -5.367229  1.651307
    [22,] -3.083278  1.971678
    [23,] -3.457198  4.146373
    [24,] -4.591451  3.608344
    [25,] -3.739405  2.961757
    [26,] -2.881218  3.148562
    [27,] -2.214021  3.166248
    [28,] -1.039581  4.117091
    [29,] -3.188470  3.195898
    [30,] -3.637721  2.673002
    [31,]  2.673002 -3.637721
    [32,]  3.195898 -3.188470
    [33,]  4.117091 -1.039581
    [34,]  3.166248 -2.214021
    [35,]  3.148562 -2.881218
    [36,]  2.961757 -3.739405
    [37,]  3.608344 -4.591451
    [38,]  4.146373 -3.457198
    [39,]  1.971678 -3.083278
    [40,]  1.651307 -5.367229
    [41,]  2.877253 -1.615082
    [42,]  4.606236 -3.189386
    [43,]  3.008510 -3.493054
    [44,]  3.196438 -3.599175
    [45,]  1.596865 -4.142593
    [46,]  3.404251 -3.481751
    [47,]  3.447379 -2.239879
    [48,]  4.885257 -1.097882
    [49,]  2.383132 -3.283969
    [50,]  3.330045 -2.788996
    [51,]  1.347806 -2.555539
    [52,]  4.612891 -2.496261
    [53,]  3.734322 -4.077324
    [54,]  4.495104 -3.599468
    [55,]  3.468711 -1.876462
    [56,]  3.387088 -2.700320
    [57,]  2.271456 -3.549683
    [58,]  1.683888 -3.856196
    [59,]  3.305911 -2.664620
    [60,]  3.081494 -4.695592

Make a plot of `x`

``` r
plot(x)
```

![](class07.markdown_strict_files/figure-markdown_strict/unnamed-chunk-4-1.png)

### K-means

The main function in “base” R for K-means clustering is called
`kmeans()`: There are two required arguments in this function and these
are x which is the data and centers which is the number of cluster you
want.

``` r
km<- kmeans(x,centers=2)
km
```

    K-means clustering with 2 clusters of sizes 30, 30

    Cluster means:
              x         y
    1 -3.140093  3.158810
    2  3.158810 -3.140093

    Clustering vector:
     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

    Within cluster sum of squares by cluster:
    [1] 54.08501 54.08501
     (between_SS / total_SS =  91.7 %)

    Available components:

    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

The `kmean()` function returns a “list” with 9 components , YOu can see
named components of any list with the `attributes()` function.

``` r
attributes(km)
```

    $names
    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

    $class
    [1] "kmeans"

> Q. How many points are in each cluster

The number you get is the number of points in each cluster.

``` r
km$size
```

    [1] 30 30

> Q. Cluster assignment/membership vector?

``` r
km$cluster
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

> Q. Cluster center

``` r
km$cluster
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

> Q. Make a plot of our `kmeans()` results showing cluster assignments
> using different color for each cluster/group of points and cluster
> center.

When we write the code `col=c("purple", "skyblue")`, it will tell R to
alternate the color so one point purple and then the next skyblue, and
we can’t differentiate the two clusters. To do this, we will write the
code `col=km$cluster` and that will color code our two cluster. pch
refers to the shape of the point and then cex is the size of the point.

``` r
plot(x, col=km$cluster)
points(km$centers, col="skyblue", pch=15, cex=1.5)
```

![](class07.markdown_strict_files/figure-markdown_strict/unnamed-chunk-10-1.png)

> Q. Run `kmeans()` again on `x` and this cluster into 4 groups/clusters
> and plot the smae results like the figure as above.

``` r
km2<- kmeans(x,centers=4)
km2
```

    K-means clustering with 4 clusters of sizes 7, 15, 8, 30

    Cluster means:
              x         y
    1 -1.797024  3.796404
    2 -3.476495  3.507356
    3 -3.684526  1.947392
    4  3.158810 -3.140093

    Clustering vector:
     [1] 2 2 3 3 2 1 2 2 1 3 2 3 1 1 2 3 2 2 2 1 3 3 2 2 2 2 1 1 2 3 4 4 4 4 4 4 4 4
    [39] 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4

    Within cluster sum of squares by cluster:
    [1]  5.387102  9.244921  6.349174 54.085014
     (between_SS / total_SS =  94.2 %)

    Available components:

    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

``` r
plot(x, col=km$cluster)
points(km$centers, col="skyblue", pch=15, cex=1.4)
```

![](class07.markdown_strict_files/figure-markdown_strict/unnamed-chunk-12-1.png)

> **key-point**: K-means clustering is super popular but can be misused.
> One big limitation is that it can impose a clustering pattern on your
> data even if clear natural groupings do not exist. i.e its does what
> you tell it to do in terms of `center`

### Hierarchical Clustering

The main function in “base” R for hierarchical clustering is called
`hclust()`. This functions has one required input d.You can’t just pass
our data set as is into `hclust()` you must give “distance matrix” as
input. We can get this fro the `dist()` function in R.

``` r
d<- dist(x)
hc <- hclust(d)
hc
```


    Call:
    hclust(d = d)

    Cluster method   : complete 
    Distance         : euclidean 
    Number of objects: 60 

The results of \``hclust` don’t have a useful `print()` method but do
have a special `plot()` method. When we plot it we get this speical
graph that looks like a phylogentic tree.

``` r
plot(hc)
abline(h=8, col="pink")
```

![](class07.markdown_strict_files/figure-markdown_strict/unnamed-chunk-14-1.png)

To get our main cluster assignment (membership vector) we need to “cut”
tree at the big goal posts. The code for this would be cutree()

``` r
grps<- cutree(hc, h=8)
grps
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
table(grps)
```

    grps
     1  2 
    30 30 

``` r
plot(x,col=grps)
```

![](class07.markdown_strict_files/figure-markdown_strict/unnamed-chunk-17-1.png)

Hierarchical Clustering is distinct in that the dendrogram (tree figure)
can reveal the potential grouping in your data (unlike K-means)

## Principle Component Analysis (PCA)

PCA is a common and highly useful dimensionality reduction technique
used in many fields, particularly bioinformatics.

Here we will anlayze some data from the UK on food consumption

### Data import

``` r
url <- "https://tinyurl.com/UK-foods"
x <- read.csv(url)
head(x)
```

                   X England Wales Scotland N.Ireland
    1         Cheese     105   103      103        66
    2  Carcass_meat      245   227      242       267
    3    Other_meat      685   803      750       586
    4           Fish     147   160      122        93
    5 Fats_and_oils      193   235      184       209
    6         Sugars     156   175      147       139

This way is not ideal because if we keep running it, we will get rid of
one of the columns after every run…

``` r
rownames(x) <-x[,1]
x<- x[,-1]
head(x)
```

                   England Wales Scotland N.Ireland
    Cheese             105   103      103        66
    Carcass_meat       245   227      242       267
    Other_meat         685   803      750       586
    Fish               147   160      122        93
    Fats_and_oils      193   235      184       209
    Sugars             156   175      147       139

``` r
rownames(x) <-x[,1]
x<- x[,-1]
head(x)
```

        Wales Scotland N.Ireland
    105   103      103        66
    245   227      242       267
    685   803      750       586
    147   160      122        93
    193   235      184       209
    156   175      147       139

A way to fix it is do what is this is to write this code

``` r
x <- read.csv(url, row.names=1)
head(x)
```

                   England Wales Scotland N.Ireland
    Cheese             105   103      103        66
    Carcass_meat       245   227      242       267
    Other_meat         685   803      750       586
    Fish               147   160      122        93
    Fats_and_oils      193   235      184       209
    Sugars             156   175      147       139

``` r
barplot(as.matrix(x), beside=F, col=rainbow(nrow(x)))
```

![](class07.markdown_strict_files/figure-markdown_strict/unnamed-chunk-22-1.png)

One conventional plot that can be useful is called a “pairs” plot

``` r
pairs (x, col=rainbow(nrow(x)), pch=16)
```

![](class07.markdown_strict_files/figure-markdown_strict/unnamed-chunk-23-1.png)

For the second graph in the diagram, the y-axis is England, and the
x-axis is Wales. Each dot on the graph represent the differnt food we
are looking at in the data

### PCA to the rescue

The main functions in the base R for PCA is called `prcomp()`

``` r
pca <- prcomp(t(x))
summary(pca)
```

    Importance of components:
                                PC1      PC2      PC3       PC4
    Standard deviation     324.1502 212.7478 73.87622 2.921e-14
    Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
    Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00

The `prcomp()` function returns a list of options for our results with
fice attributes/components.

``` r
attributes(pca)
```

    $names
    [1] "sdev"     "rotation" "center"   "scale"    "x"       

    $class
    [1] "prcomp"

The two main “results” inhere are `pca$x` adn `pca$rotation`. The first
of these (`pca$x`) contains the scores of the data on the new PC axis-
we use these to make our “PCA plot”

``` r
pca$x
```

                     PC1         PC2        PC3           PC4
    England   -144.99315   -2.532999 105.768945 -9.152022e-15
    Wales     -240.52915 -224.646925 -56.475555  5.560040e-13
    Scotland   -91.86934  286.081786 -44.415495 -6.638419e-13
    N.Ireland  477.39164  -58.901862  -4.877895  1.329771e-13

``` r
library(ggplot2)
library (ggrepel)

#Make a plot of pca$x with PC1 vs PC2
ggplot(pca$x) +
  aes(PC1, PC2, label=rownames(pca$x)) +
  geom_point()+
  geom_text_repel()
```

![](class07.markdown_strict_files/figure-markdown_strict/unnamed-chunk-27-1.png)

The plot above is showing a plot that is comparing two different
dimensions, PC1 and PC2. we are looking at the variance.

The second major result is curated in the `pca$rotaation` object or
component. Lets plot this to see what PCA is picking up

``` r
pca$rotation
```

                                 PC1          PC2         PC3          PC4
    Cheese              -0.056955380  0.016012850  0.02394295 -0.409382587
    Carcass_meat         0.047927628  0.013915823  0.06367111  0.729481922
    Other_meat          -0.258916658 -0.015331138 -0.55384854  0.331001134
    Fish                -0.084414983 -0.050754947  0.03906481  0.022375878
    Fats_and_oils       -0.005193623 -0.095388656 -0.12522257  0.034512161
    Sugars              -0.037620983 -0.043021699 -0.03605745  0.024943337
    Fresh_potatoes       0.401402060 -0.715017078 -0.20668248  0.021396007
    Fresh_Veg           -0.151849942 -0.144900268  0.21382237  0.001606882
    Other_Veg           -0.243593729 -0.225450923 -0.05332841  0.031153231
    Processed_potatoes  -0.026886233  0.042850761 -0.07364902 -0.017379680
    Processed_Veg       -0.036488269 -0.045451802  0.05289191  0.021250980
    Fresh_fruit         -0.632640898 -0.177740743  0.40012865  0.227657348
    Cereals             -0.047702858 -0.212599678 -0.35884921  0.100043319
    Beverages           -0.026187756 -0.030560542 -0.04135860 -0.018382072
    Soft_drinks          0.232244140  0.555124311 -0.16942648  0.222319484
    Alcoholic_drinks    -0.463968168  0.113536523 -0.49858320 -0.273126013
    Confectionery       -0.029650201  0.005949921 -0.05232164  0.001890737

``` r
ggplot(pca$rotation)+
  aes(PC1, rownames(pca$rotation))+ 
  geom_col()
```

![](class07.markdown_strict_files/figure-markdown_strict/unnamed-chunk-29-1.png)

The barplot above shows the variance in the different food
consumptions.PC1 is a dimension, and it will show the greatest amount of
variance.

The PC number can be no more than the the number of food and also the
countries.
