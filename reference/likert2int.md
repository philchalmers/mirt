# Convert ordered Likert-scale responses (character or factors) to integers

Given a matrix or data.frame object consisting of Likert responses
return an object of the same dimensions with integer values.

## Usage

``` r
likert2int(x, levels = NULL)
```

## Arguments

- x:

  a matrix of character values or data.frame of character/factor vectors

- levels:

  a named character vector indicating which integer values should be
  assigned to which elements. If omitted, the order of the elements will
  be determined after converting each column in `x` to a factor variable

## References

Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## See also

[`key2binary`](https://philchalmers.github.io/mirt/reference/key2binary.md),
[`poly2dich`](https://philchalmers.github.io/mirt/reference/poly2dich.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

# simulate data

dat1 <- matrix(sample(c('Disagree', 'Strongly Disagree', 'Agree',
                        'Neutral', 'Strongly Agree'), 1000*5, replace=TRUE),
               nrow=1000, ncol=5)
dat1[2,2] <- dat1[3,3] <- dat1[1,3] <- NA # NAs added for flavour
dat2 <- matrix(sample(c('D', 'SD', 'A', 'N', 'SA'), 1000*5, replace=TRUE),
               nrow=1000, ncol=5)
dat <- cbind(dat1, dat2)

# separately
intdat1 <- likert2int(dat1)
head(dat1)
#>      [,1]                [,2]                [,3]               
#> [1,] "Strongly Disagree" "Strongly Disagree" NA                 
#> [2,] "Agree"             NA                  "Disagree"         
#> [3,] "Strongly Agree"    "Strongly Disagree" NA                 
#> [4,] "Neutral"           "Strongly Disagree" "Strongly Agree"   
#> [5,] "Disagree"          "Agree"             "Agree"            
#> [6,] "Neutral"           "Strongly Disagree" "Strongly Disagree"
#>      [,4]             [,5]               
#> [1,] "Agree"          "Agree"            
#> [2,] "Neutral"        "Agree"            
#> [3,] "Disagree"       "Agree"            
#> [4,] "Strongly Agree" "Strongly Agree"   
#> [5,] "Neutral"        "Strongly Disagree"
#> [6,] "Strongly Agree" "Strongly Agree"   
head(intdat1)
#>   V1 V2 V3 V4 V5
#> 1 NA NA NA NA NA
#> 2 NA NA NA NA NA
#> 3 NA NA NA NA NA
#> 4 NA NA NA NA NA
#> 5 NA NA NA NA NA
#> 6 NA NA NA NA NA

# more useful with explicit levels
lvl1 <- c('Strongly Disagree'=1, 'Disagree'=2, 'Neutral'=3, 'Agree'=4,
          'Strongly Agree'=5)
intdat1 <- likert2int(dat1, levels = lvl1)
head(dat1)
#>      [,1]                [,2]                [,3]               
#> [1,] "Strongly Disagree" "Strongly Disagree" NA                 
#> [2,] "Agree"             NA                  "Disagree"         
#> [3,] "Strongly Agree"    "Strongly Disagree" NA                 
#> [4,] "Neutral"           "Strongly Disagree" "Strongly Agree"   
#> [5,] "Disagree"          "Agree"             "Agree"            
#> [6,] "Neutral"           "Strongly Disagree" "Strongly Disagree"
#>      [,4]             [,5]               
#> [1,] "Agree"          "Agree"            
#> [2,] "Neutral"        "Agree"            
#> [3,] "Disagree"       "Agree"            
#> [4,] "Strongly Agree" "Strongly Agree"   
#> [5,] "Neutral"        "Strongly Disagree"
#> [6,] "Strongly Agree" "Strongly Agree"   
head(intdat1)
#>   V1 V2 V3 V4 V5
#> 1  1  1 NA  4  4
#> 2  4 NA  2  3  4
#> 3  5  1 NA  2  4
#> 4  3  1  5  5  5
#> 5  2  4  4  3  1
#> 6  3  1  1  5  5

# second data
lvl2 <- c('SD'=1, 'D'=2, 'N'=3, 'A'=4, 'SA'=5)
intdat2 <- likert2int(dat2, levels = lvl2)
head(dat2)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] "SD" "N"  "D"  "SD" "SA"
#> [2,] "SD" "N"  "D"  "N"  "N" 
#> [3,] "D"  "A"  "SD" "SA" "A" 
#> [4,] "N"  "D"  "A"  "SA" "SD"
#> [5,] "D"  "SA" "D"  "D"  "SD"
#> [6,] "SA" "N"  "A"  "A"  "SD"
head(intdat2)
#>   V1 V2 V3 V4 V5
#> 1  1  3  2  1  5
#> 2  1  3  2  3  3
#> 3  2  4  1  5  4
#> 4  3  2  4  5  1
#> 5  2  5  2  2  1
#> 6  5  3  4  4  1

# full dataset (using both mapping schemes)
intdat <- likert2int(dat, levels = c(lvl1, lvl2))
head(dat)
#>      [,1]                [,2]                [,3]               
#> [1,] "Strongly Disagree" "Strongly Disagree" NA                 
#> [2,] "Agree"             NA                  "Disagree"         
#> [3,] "Strongly Agree"    "Strongly Disagree" NA                 
#> [4,] "Neutral"           "Strongly Disagree" "Strongly Agree"   
#> [5,] "Disagree"          "Agree"             "Agree"            
#> [6,] "Neutral"           "Strongly Disagree" "Strongly Disagree"
#>      [,4]             [,5]                [,6] [,7] [,8] [,9] [,10]
#> [1,] "Agree"          "Agree"             "SD" "N"  "D"  "SD" "SA" 
#> [2,] "Neutral"        "Agree"             "SD" "N"  "D"  "N"  "N"  
#> [3,] "Disagree"       "Agree"             "D"  "A"  "SD" "SA" "A"  
#> [4,] "Strongly Agree" "Strongly Agree"    "N"  "D"  "A"  "SA" "SD" 
#> [5,] "Neutral"        "Strongly Disagree" "D"  "SA" "D"  "D"  "SD" 
#> [6,] "Strongly Agree" "Strongly Agree"    "SA" "N"  "A"  "A"  "SD" 
head(intdat)
#>   V1 V2 V3 V4 V5 V6 V7 V8 V9 V10
#> 1  1  1 NA  4  4  1  3  2  1   5
#> 2  4 NA  2  3  4  1  3  2  3   3
#> 3  5  1 NA  2  4  2  4  1  5   4
#> 4  3  1  5  5  5  3  2  4  5   1
#> 5  2  4  4  3  1  2  5  2  2   1
#> 6  3  1  1  5  5  5  3  4  4   1


#####
# data.frame as input with ordered factors

dat1 <- data.frame(dat1)
dat2 <- data.frame(dat2)
dat.old <- cbind(dat1, dat2)
colnames(dat.old) <- paste0('Item_', 1:10)
str(dat.old) # factors are leveled alphabetically by default
#> 'data.frame':    1000 obs. of  10 variables:
#>  $ Item_1 : chr  "Strongly Disagree" "Agree" "Strongly Agree" "Neutral" ...
#>  $ Item_2 : chr  "Strongly Disagree" NA "Strongly Disagree" "Strongly Disagree" ...
#>  $ Item_3 : chr  NA "Disagree" NA "Strongly Agree" ...
#>  $ Item_4 : chr  "Agree" "Neutral" "Disagree" "Strongly Agree" ...
#>  $ Item_5 : chr  "Agree" "Agree" "Agree" "Strongly Agree" ...
#>  $ Item_6 : chr  "SD" "SD" "D" "N" ...
#>  $ Item_7 : chr  "N" "N" "A" "D" ...
#>  $ Item_8 : chr  "D" "D" "SD" "A" ...
#>  $ Item_9 : chr  "SD" "N" "SA" "SA" ...
#>  $ Item_10: chr  "SA" "N" "A" "SD" ...

# create explicit ordering in factor variables
for(i in 1:ncol(dat1))
   levels(dat1[[i]]) <- c('Strongly Disagree', 'Disagree', 'Neutral', 'Agree',
                          'Strongly Agree')

for(i in 1:ncol(dat2))
   levels(dat2[[i]]) <- c('SD', 'D', 'N', 'A', 'SA')

dat <- cbind(dat1, dat2)
colnames(dat) <- colnames(dat.old)
str(dat) # note ordering
#> 'data.frame':    1000 obs. of  10 variables:
#>  $ Item_1 : chr  "Strongly Disagree" "Agree" "Strongly Agree" "Neutral" ...
#>   ..- attr(*, "levels")= chr [1:5] "Strongly Disagree" "Disagree" "Neutral" "Agree" ...
#>  $ Item_2 : chr  "Strongly Disagree" NA "Strongly Disagree" "Strongly Disagree" ...
#>   ..- attr(*, "levels")= chr [1:5] "Strongly Disagree" "Disagree" "Neutral" "Agree" ...
#>  $ Item_3 : chr  NA "Disagree" NA "Strongly Agree" ...
#>   ..- attr(*, "levels")= chr [1:5] "Strongly Disagree" "Disagree" "Neutral" "Agree" ...
#>  $ Item_4 : chr  "Agree" "Neutral" "Disagree" "Strongly Agree" ...
#>   ..- attr(*, "levels")= chr [1:5] "Strongly Disagree" "Disagree" "Neutral" "Agree" ...
#>  $ Item_5 : chr  "Agree" "Agree" "Agree" "Strongly Agree" ...
#>   ..- attr(*, "levels")= chr [1:5] "Strongly Disagree" "Disagree" "Neutral" "Agree" ...
#>  $ Item_6 : chr  "SD" "SD" "D" "N" ...
#>   ..- attr(*, "levels")= chr [1:5] "SD" "D" "N" "A" ...
#>  $ Item_7 : chr  "N" "N" "A" "D" ...
#>   ..- attr(*, "levels")= chr [1:5] "SD" "D" "N" "A" ...
#>  $ Item_8 : chr  "D" "D" "SD" "A" ...
#>   ..- attr(*, "levels")= chr [1:5] "SD" "D" "N" "A" ...
#>  $ Item_9 : chr  "SD" "N" "SA" "SA" ...
#>   ..- attr(*, "levels")= chr [1:5] "SD" "D" "N" "A" ...
#>  $ Item_10: chr  "SA" "N" "A" "SD" ...
#>   ..- attr(*, "levels")= chr [1:5] "SD" "D" "N" "A" ...

intdat <- likert2int(dat)
head(dat)
#>              Item_1            Item_2            Item_3         Item_4
#> 1 Strongly Disagree Strongly Disagree              <NA>          Agree
#> 2             Agree              <NA>          Disagree        Neutral
#> 3    Strongly Agree Strongly Disagree              <NA>       Disagree
#> 4           Neutral Strongly Disagree    Strongly Agree Strongly Agree
#> 5          Disagree             Agree             Agree        Neutral
#> 6           Neutral Strongly Disagree Strongly Disagree Strongly Agree
#>              Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> 1             Agree     SD      N      D     SD      SA
#> 2             Agree     SD      N      D      N       N
#> 3             Agree      D      A     SD     SA       A
#> 4    Strongly Agree      N      D      A     SA      SD
#> 5 Strongly Disagree      D     SA      D      D      SD
#> 6    Strongly Agree     SA      N      A      A      SD
head(intdat)
#>   Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> 1      1      1     NA      4      4      1      3      2      1       5
#> 2      4     NA      2      3      4      1      3      2      3       3
#> 3      5      1     NA      2      4      2      4      1      5       4
#> 4      3      1      5      5      5      3      2      4      5       1
#> 5      2      4      4      3      1      2      5      2      2       1
#> 6      3      1      1      5      5      5      3      4      4       1

# }
```
