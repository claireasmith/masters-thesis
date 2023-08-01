Testing for intraspecific pollen size-number tradeoffs
================
Claire Smith
2023-07-28

### Goal

Testing for a within-species tradeoff in pollen size and number using
linear regression: pollen number ~ pollen volume.

If n = number of grains per anther, s = pollen grain size, and e is the
average resources for pollen production per flower, then we should
expect

$$n \propto e s^a $$ where a is a constant determining the shape of the
relationship between pollen number and size.

If we let c be some constant converting pollen size to number, we have

$$ n = c s^a $$

If pollen size and number trade off, a should be negative. Since pollen
diameter is our measure for s, and pollen from wind-pollinated plants
tends to be spherical, we should expect s to be related to n by a power
of 3. So if pollen diameter and number trade off, we should expect a to
be equal to (or less than) -3.

I’ll test for a tradeoff between pollen diameter and number by
log-transforming n and s,

$$log(n) = log(c s^a) = log(c) + a log(s) $$ and estimating a and the
intercept log(c) using linear regression.

``` r
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
source("theme_cs.R")
```

``` r
# Read in individual-level pollen size and number data
prodfull <- read.csv("processed-data/size-prod-norep.csv", stringsAsFactors = T) 
# head(prodfull)
# str(prodfull)
# summary(prodfull)

# Keep only entries with both size and prod values
prod <- prodfull %>% 
  filter(!is.na(Avg_diam) & !is.na(Avg_pol_anth))

# Add log-transformed pollen per anther to dataframe
prod$Log_avg_pol_anth <- log(prod$Avg_pol_anth)
prod$Log_avg_diam <- log(prod$Avg_diam)

# Sort species alphabetically
prod <- prod %>% arrange(Species)

# Remove outliers identified below
prod <- prod %>% 
  filter(Species == "Plantago lanceolata" & Ind!="q" | Species !="Plantago lanceolata") %>% 
  filter(Species == "Phleum pratense" & Ind!="C17" | Species !="Phleum pratense") %>% 
  filter(Species == "Hierochloe odorata" & Ind!="08P" | Species !="Hierochloe odorata") %>% 
  filter(Species == "Festuca campestris" & Ind!="FF07" | Species !="Festuca campestris")

# Nest data by species to make it easier to handle
by_spp <- prod %>% group_by(Species, Sex_sys) %>% 
  nest() %>% 
  arrange(Species) %>% 
  # Keep only species with at least 10 entries %>% 
  mutate(nind = purrr::map(data, nrow)) %>% 
  filter(nind>=10)
# 23 species remaining
```

To test for tradeoffs in pollen size and number within species, I will
build a linear model for each species predicting pollen production by
pollen size.

First I’ll go through each species’ data and look at the distributions
of pollen size and number to see how symmetric they are and if there are
any obvious outliers.

``` r
#inspect distributions of mean pollen diam and pollen production for each species - are they symmetric? any obvious outliers? 
#define function for this
checkfun <- function(df,i){
  spec = as.character(df[i,]$Species)
  n=length(df[i,]$data[[1]]$Avg_diam)
  par(mfrow = c(2, 2))
  hist(df[i,]$data[[1]]$Avg_diam, main=paste(spec, "n=",n), xlab="Mean pollen diameter")
  hist(log(df[i,]$data[[1]]$Avg_diam), main=paste(spec, "n=",n), xlab="Log mean pollen diameter")
  hist(df[i,]$data[[1]]$Avg_pol_anth, main=paste(spec, "n=",n), xlab="Mean pollen per anther")
  hist(log(df[i,]$data[[1]]$Avg_pol_anth), main=paste(spec, "n=",n), xlab="Log mean pollen per anther")
  par(mfrow=c(1,1))
}

# Run through all species
# by_spp$Species
checkfun(by_spp,13)
```

![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/inspect%20data-1.png)<!-- -->

``` r
# Pollen production per anther tends to be right-skewed 
```

I’ll plot all the data as pollen production vs size to get an idea of
what it looks like, and if there’s any clear outliers or potential
grouping in the data.

``` r
prod %>% 
  filter(Species %in% by_spp$Species) %>% 
  ggplot(aes(x=Avg_diam, y=Avg_pol_anth)) + 
  geom_point() + 
  geom_smooth(method = "glm", formula = y~log(x),
              method.args = list(family = gaussian(link = 'log')),
              se=F) +
  facet_wrap(~Species, scales = "free", ncol=4) + 
  scale_y_continuous(name="Pollen production per anther") + 
  scale_x_continuous(name=expression(paste("Pollen diameter (", mu, "m",")"))) + 
  theme_cs(font = "sans", fontsize=20) +
  theme(strip.background = element_rect(linewidth = NULL,
                                        linetype = NULL,
                                        colour = "white"),
        strip.text = element_text(size = 18, face="italic"))
```

![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/all%20tradeoff%20plots-1.png)<!-- -->

Now I’ll run linear models predicting pollen production by pollen size

``` r
#define function to run model 
spp_model <- function(df) {
  lm(Log_avg_pol_anth ~ Log_avg_diam, data = df)
}
#run models and add model to by_spp dataframe
by_spp <- by_spp %>% 
  mutate(model = purrr::map(data, spp_model)) %>% 
  arrange(Sex_sys, Species) # arrange in order of sex sys, then alphabetical
#get model summary using glance() from "broom"
glance <- by_spp %>% 
  mutate(glance = purrr::map(model, broom::glance)) %>% 
  unnest(glance, .drop=T)
```

    ## Warning: The `.drop` argument of `unnest()` is deprecated as of tidyr 1.0.0.
    ## ℹ All list-columns are now preserved.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
# put this in a viewable format (remove nested data and model info)
glancetab <- glance %>% select(-c(data, model)) %>% 
  arrange(Sex_sys, Species)
#save this in a new table without the ugly list elements, which I can print to a file
glancetab_write <- glance %>% select(-c(data, model)) %>% 
  arrange(Sex_sys, Species) %>% apply(2, as.character) # coercing all columns to characters so that I can write this to a file
write.csv(glancetab_write, "processed-data/tradeoff-reg-model-sum.csv", row.names = F)
#take a look
print(glancetab, n=Inf)
```

    ## # A tibble: 23 × 15
    ## # Groups:   Species, Sex_sys [23]
    ##    Sex_sys   Species nind  r.squared adj.r.squared sigma statistic p.value    df
    ##    <fct>     <fct>   <lis>     <dbl>         <dbl> <dbl>     <dbl>   <dbl> <dbl>
    ##  1 dioecious Rumex … <int>  0.000544      -0.0521  0.290    0.0103 9.20e-1     1
    ##  2 dioecious Thalic… <int>  0.0321        -0.0140  0.493    0.697  4.13e-1     1
    ##  3 hermaphr… Bromus… <int>  0.326          0.302   0.496   13.6    9.79e-4     1
    ##  4 hermaphr… Chenop… <int>  0.117          0.0727  0.261    2.65   1.20e-1     1
    ##  5 hermaphr… Elymus… <int>  0.0154        -0.0163  0.460    0.486  4.91e-1     1
    ##  6 hermaphr… Elymus… <int>  0.177          0.0740  0.190    1.72   2.26e-1     1
    ##  7 hermaphr… Festuc… <int>  0.142          0.127   0.321    9.45   3.24e-3     1
    ##  8 hermaphr… Festuc… <int>  0.410          0.361   0.312    8.34   1.36e-2     1
    ##  9 hermaphr… Hieroc… <int>  0.00714       -0.0140  0.596    0.338  5.64e-1     1
    ## 10 hermaphr… Koeler… <int>  0.275          0.243   0.417    8.72   7.12e-3     1
    ## 11 hermaphr… Leymus… <int>  0.201          0.176   0.307    8.28   6.97e-3     1
    ## 12 hermaphr… Phleum… <int>  0.00342       -0.0335  0.372    0.0928 7.63e-1     1
    ## 13 hermaphr… Planta… <int>  0.284          0.250   0.203    8.35   8.78e-3     1
    ## 14 hermaphr… Schiza… <int>  0.0158        -0.0826  0.380    0.161  6.97e-1     1
    ## 15 hermaphr… Stipa … <int>  0.273          0.182   0.646    3.00   1.22e-1     1
    ## 16 monoecio… Amaran… <int>  0.369          0.324   0.209    8.18   1.26e-2     1
    ## 17 monoecio… Ambros… <int>  0.140          0.0896  0.506    2.77   1.14e-1     1
    ## 18 monoecio… Carex … <int>  0.111          0.0475  0.138    1.75   2.07e-1     1
    ## 19 monoecio… Carex … <int>  0.0456        -0.0141  0.243    0.764  3.95e-1     1
    ## 20 monoecio… Carex … <int>  0.126          0.0710  0.214    2.30   1.49e-1     1
    ## 21 monoecio… Carex … <int>  0.0267        -0.00278 0.179    0.906  3.48e-1     1
    ## 22 monoecio… Rumex … <int>  0.0227        -0.0316  0.219    0.418  5.26e-1     1
    ## 23 monoecio… Scirpu… <int>  0.00756       -0.0751  0.342    0.0914 7.68e-1     1
    ## # ℹ 6 more variables: logLik <dbl>, AIC <dbl>, BIC <dbl>, deviance <dbl>,
    ## #   df.residual <int>, nobs <int>

``` r
#similar - but using broom::tidy(), gives me model coefficients
tidy <- by_spp %>% 
  mutate(tidy = purrr::map(model, broom::tidy)) %>% 
  unnest(tidy, .drop=T)
```

    ## Warning: The `.drop` argument of `unnest()` is deprecated as of tidyr 1.0.0.
    ## ℹ All list-columns are now preserved.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
tidytab <- tidy %>% select(-c(data, model)) %>% 
  arrange(Sex_sys, Species)
tidytab_write <- tidy %>% select(-c(data, model)) %>% 
  arrange(Sex_sys, Species) %>% apply(2, as.character)
write.csv(tidytab_write, "processed-data/tradeoff-reg-model-coefs.csv", row.names = F)
#take a look
# print(tidytab, n=Inf)

# Which species have estimates <= -0? 
tidytab[which(tidytab$term == "Log_avg_diam" & tidytab$estimate <= 0 & tidytab$p.value<0.05),]
```

    ## # A tibble: 6 × 8
    ## # Groups:   Species, Sex_sys [6]
    ##   Sex_sys        Species        nind  term  estimate std.error statistic p.value
    ##   <fct>          <fct>          <lis> <chr>    <dbl>     <dbl>     <dbl>   <dbl>
    ## 1 hermaphroditic Bromus inermis <int> Log_…    -4.51     1.23      -3.68 9.79e-4
    ## 2 hermaphroditic Festuca campe… <int> Log_…    -2.36     0.767     -3.07 3.24e-3
    ## 3 hermaphroditic Festuca prate… <int> Log_…    -6.51     2.25      -2.89 1.36e-2
    ## 4 hermaphroditic Koeleria cris… <int> Log_…    -4.86     1.65      -2.95 7.12e-3
    ## 5 hermaphroditic Plantago lanc… <int> Log_…    -2.14     0.740     -2.89 8.78e-3
    ## 6 monoecious     Amaranthus re… <int> Log_…    -2.30     0.803     -2.86 1.26e-2

``` r
# 6 species: Bromus inermis, Chenopodium album, Festuca pratensis, Koeleria cristata,
# Ambrosia artemisiifolia, and Carex hirtifolia
# Any >=0?
tidytab[which(tidytab$term == "Log_avg_diam" & tidytab$estimate >= 3 & tidytab$p.value<0.05),]
```

    ## # A tibble: 1 × 8
    ## # Groups:   Species, Sex_sys [1]
    ##   Sex_sys        Species        nind  term  estimate std.error statistic p.value
    ##   <fct>          <fct>          <lis> <chr>    <dbl>     <dbl>     <dbl>   <dbl>
    ## 1 hermaphroditic Leymus innova… <int> Log_…     3.25      1.13      2.88 0.00697

``` r
# 2 species: Leymus innovatus and Stipa columbiana
```

Check model assumptions, look for outliers:

``` r
## Test assumptions for each linear model
#define function for this
checkassum <- function(df,spec){
  dat = df %>% filter(Species == spec)
  mod = lm(Log_avg_pol_anth ~ Log_avg_diam, data = dat)
  par(mfrow = c(2, 2))
  plot(mod, main=as.character(spec))
  par(mfrow=c(1,1))
}

prod2 <- prod %>% filter(Species %in% by_spp$Species) %>% droplevels()
spvec <- levels(prod2$Species)
for (spec in spvec) {
  checkassum(prod, spec)
}
```

![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-1.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-2.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-3.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-4.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-5.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-6.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-7.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-8.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-9.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-10.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-11.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-12.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-13.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-14.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-15.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-16.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-17.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-18.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-19.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-20.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-21.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-22.png)<!-- -->![](Q1-2_within-spp-sizenum-tradeoff_files/figure-gfm/test%20model%20assumptions-23.png)<!-- -->

``` r
# Quite a few species look like there's deviations from normality (from QQ line) near extremes of residuals

# point 20 in Plantago lanceolata has very negative standardized residual (<-3), almost 0.5 cook's d
# will remove
# point 5 Phleum pratense also has standardized residual <-3
# Hierochloe odorata 5 has standardized residual <-4!! 
# Festuca campestris point 7 also has standardized residual <-4
# Also seems to be a point with unusually high leverage in A artemisiifolia

# test = prod %>% filter(Species == "Plantago lanceolata")
# # point 20 Ind="q"
# test = prod %>% filter(Species == "Phleum pratense")
# # point 5 Ind = C17
# test = prod %>% filter(Species == "Hierochloe odorata")
# # point 5 Ind=08P -- 50 grains per anther! seems like an error
# test = prod %>% filter(Species == "Festuca campestris")
# # point 7 Ind=FF07, also very low grains per anther
# 
# filter(Species == "Plantago lanceolata" & Ind!="q" | Species !="Plantago lanceolata") %>% 
# filter(Species == "Phleum pratense" & Ind!="C17" | Species !="Phleum pratense") %>% 
# filter(Species == "Hierochloe odorata" & Ind!="08P" | Species !="Hierochloe odorata") %>% 
# filter(Species == "Festuca campestris" & Ind!="FF07" | Species !="Festuca campestris")
# 
# 
# test = prod %>% filter(Species == "Ambrosia artemisiifolia")
# library(ggfortify)
# autoplot(tmod)
# 
# model.diag.metrics <- broom::augment(tmod)
# head(model.diag.metrics)
# 
# model.diag.metrics <- model.diag.metrics %>%
#   mutate(index = 1:nrow(model.diag.metrics)) %>%
#   select(index, everything(), -.sigma)
# # point 4 has very high hat values relative to everything else - 0.5
# test
# # Individual 4 has the largest pollen measured for this species: 19 um, with high sd. It also has
# # relatively few N diam measurements (~6000), which isn't the lowest reported but is a lot lower
# # than other inds. It has relatively low pollen per anther but again not the lowest. 
# 
# 
```
