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

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
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
#>      [,1]             [,2]                [,3]             [,4]               
#> [1,] "Strongly Agree" "Strongly Disagree" NA               "Neutral"          
#> [2,] "Disagree"       NA                  "Agree"          "Disagree"         
#> [3,] "Disagree"       "Neutral"           NA               "Disagree"         
#> [4,] "Disagree"       "Neutral"           "Disagree"       "Neutral"          
#> [5,] "Neutral"        "Strongly Agree"    "Agree"          "Strongly Disagree"
#> [6,] "Strongly Agree" "Agree"             "Strongly Agree" "Strongly Agree"   
#>      [,5]               
#> [1,] "Strongly Agree"   
#> [2,] "Strongly Agree"   
#> [3,] "Strongly Disagree"
#> [4,] "Neutral"          
#> [5,] "Strongly Agree"   
#> [6,] "Agree"            
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
#>      [,1]             [,2]                [,3]             [,4]               
#> [1,] "Strongly Agree" "Strongly Disagree" NA               "Neutral"          
#> [2,] "Disagree"       NA                  "Agree"          "Disagree"         
#> [3,] "Disagree"       "Neutral"           NA               "Disagree"         
#> [4,] "Disagree"       "Neutral"           "Disagree"       "Neutral"          
#> [5,] "Neutral"        "Strongly Agree"    "Agree"          "Strongly Disagree"
#> [6,] "Strongly Agree" "Agree"             "Strongly Agree" "Strongly Agree"   
#>      [,5]               
#> [1,] "Strongly Agree"   
#> [2,] "Strongly Agree"   
#> [3,] "Strongly Disagree"
#> [4,] "Neutral"          
#> [5,] "Strongly Agree"   
#> [6,] "Agree"            
head(intdat1)
#>   V1 V2 V3 V4 V5
#> 1  5  1 NA  3  5
#> 2  2 NA  4  2  5
#> 3  2  3 NA  2  1
#> 4  2  3  2  3  3
#> 5  3  5  4  1  5
#> 6  5  4  5  5  4

# second data
lvl2 <- c('SD'=1, 'D'=2, 'N'=3, 'A'=4, 'SA'=5)
intdat2 <- likert2int(dat2, levels = lvl2)
head(dat2)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] "SA" "A"  "SA" "N"  "D" 
#> [2,] "D"  "D"  "SD" "SD" "D" 
#> [3,] "N"  "N"  "SD" "A"  "N" 
#> [4,] "SD" "A"  "SA" "SD" "SA"
#> [5,] "SA" "N"  "N"  "SD" "A" 
#> [6,] "A"  "D"  "SD" "SD" "SA"
head(intdat2)
#>   V1 V2 V3 V4 V5
#> 1  5  4  5  3  2
#> 2  2  2  1  1  2
#> 3  3  3  1  4  3
#> 4  1  4  5  1  5
#> 5  5  3  3  1  4
#> 6  4  2  1  1  5

# full dataset (using both mapping schemes)
intdat <- likert2int(dat, levels = c(lvl1, lvl2))
head(dat)
#>      [,1]             [,2]                [,3]             [,4]               
#> [1,] "Strongly Agree" "Strongly Disagree" NA               "Neutral"          
#> [2,] "Disagree"       NA                  "Agree"          "Disagree"         
#> [3,] "Disagree"       "Neutral"           NA               "Disagree"         
#> [4,] "Disagree"       "Neutral"           "Disagree"       "Neutral"          
#> [5,] "Neutral"        "Strongly Agree"    "Agree"          "Strongly Disagree"
#> [6,] "Strongly Agree" "Agree"             "Strongly Agree" "Strongly Agree"   
#>      [,5]                [,6] [,7] [,8] [,9] [,10]
#> [1,] "Strongly Agree"    "SA" "A"  "SA" "N"  "D"  
#> [2,] "Strongly Agree"    "D"  "D"  "SD" "SD" "D"  
#> [3,] "Strongly Disagree" "N"  "N"  "SD" "A"  "N"  
#> [4,] "Neutral"           "SD" "A"  "SA" "SD" "SA" 
#> [5,] "Strongly Agree"    "SA" "N"  "N"  "SD" "A"  
#> [6,] "Agree"             "A"  "D"  "SD" "SD" "SA" 
head(intdat)
#>   V1 V2 V3 V4 V5 V6 V7 V8 V9 V10
#> 1  5  1 NA  3  5  5  4  5  3   2
#> 2  2 NA  4  2  5  2  2  1  1   2
#> 3  2  3 NA  2  1  3  3  1  4   3
#> 4  2  3  2  3  3  1  4  5  1   5
#> 5  3  5  4  1  5  5  3  3  1   4
#> 6  5  4  5  5  4  4  2  1  1   5


#####
# data.frame as input with ordered factors

dat1 <- data.frame(dat1)
dat2 <- data.frame(dat2)
dat.old <- cbind(dat1, dat2)
colnames(dat.old) <- paste0('Item_', 1:10)
str(dat.old) # factors are leveled alphabetically by default
#> 'data.frame':    1000 obs. of  10 variables:
#>  $ Item_1 : chr  "Strongly Agree" "Disagree" "Disagree" "Disagree" ...
#>  $ Item_2 : chr  "Strongly Disagree" NA "Neutral" "Neutral" ...
#>  $ Item_3 : chr  NA "Agree" NA "Disagree" ...
#>  $ Item_4 : chr  "Neutral" "Disagree" "Disagree" "Neutral" ...
#>  $ Item_5 : chr  "Strongly Agree" "Strongly Agree" "Strongly Disagree" "Neutral" ...
#>  $ Item_6 : chr  "SA" "D" "N" "SD" ...
#>  $ Item_7 : chr  "A" "D" "N" "A" ...
#>  $ Item_8 : chr  "SA" "SD" "SD" "SA" ...
#>  $ Item_9 : chr  "N" "SD" "A" "SD" ...
#>  $ Item_10: chr  "D" "D" "N" "SA" ...

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
#>  $ Item_1 : chr  "Strongly Agree" "Disagree" "Disagree" "Disagree" ...
#>   ..- attr(*, "levels")= chr [1:5] "Strongly Disagree" "Disagree" "Neutral" "Agree" ...
#>  $ Item_2 : chr  "Strongly Disagree" NA "Neutral" "Neutral" ...
#>   ..- attr(*, "levels")= chr [1:5] "Strongly Disagree" "Disagree" "Neutral" "Agree" ...
#>  $ Item_3 : chr  NA "Agree" NA "Disagree" ...
#>   ..- attr(*, "levels")= chr [1:5] "Strongly Disagree" "Disagree" "Neutral" "Agree" ...
#>  $ Item_4 : chr  "Neutral" "Disagree" "Disagree" "Neutral" ...
#>   ..- attr(*, "levels")= chr [1:5] "Strongly Disagree" "Disagree" "Neutral" "Agree" ...
#>  $ Item_5 : chr  "Strongly Agree" "Strongly Agree" "Strongly Disagree" "Neutral" ...
#>   ..- attr(*, "levels")= chr [1:5] "Strongly Disagree" "Disagree" "Neutral" "Agree" ...
#>  $ Item_6 : chr  "SA" "D" "N" "SD" ...
#>   ..- attr(*, "levels")= chr [1:5] "SD" "D" "N" "A" ...
#>  $ Item_7 : chr  "A" "D" "N" "A" ...
#>   ..- attr(*, "levels")= chr [1:5] "SD" "D" "N" "A" ...
#>  $ Item_8 : chr  "SA" "SD" "SD" "SA" ...
#>   ..- attr(*, "levels")= chr [1:5] "SD" "D" "N" "A" ...
#>  $ Item_9 : chr  "N" "SD" "A" "SD" ...
#>   ..- attr(*, "levels")= chr [1:5] "SD" "D" "N" "A" ...
#>  $ Item_10: chr  "D" "D" "N" "SA" ...
#>   ..- attr(*, "levels")= chr [1:5] "SD" "D" "N" "A" ...

intdat <- likert2int(dat)
head(dat)
#>           Item_1            Item_2         Item_3            Item_4
#> 1 Strongly Agree Strongly Disagree           <NA>           Neutral
#> 2       Disagree              <NA>          Agree          Disagree
#> 3       Disagree           Neutral           <NA>          Disagree
#> 4       Disagree           Neutral       Disagree           Neutral
#> 5        Neutral    Strongly Agree          Agree Strongly Disagree
#> 6 Strongly Agree             Agree Strongly Agree    Strongly Agree
#>              Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> 1    Strongly Agree     SA      A     SA      N       D
#> 2    Strongly Agree      D      D     SD     SD       D
#> 3 Strongly Disagree      N      N     SD      A       N
#> 4           Neutral     SD      A     SA     SD      SA
#> 5    Strongly Agree     SA      N      N     SD       A
#> 6             Agree      A      D     SD     SD      SA
head(intdat)
#>   Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> 1      5      1     NA      3      5      5      4      5      3       2
#> 2      2     NA      4      2      5      2      2      1      1       2
#> 3      2      3     NA      2      1      3      3      1      4       3
#> 4      2      3      2      3      3      1      4      5      1       5
#> 5      3      5      4      1      5      5      3      3      1       4
#> 6      5      4      5      5      4      4      2      1      1       5

# }
```
