---
title: "Final_Week3_R_Assignment"
author: "Divya Purohith"
date: "November 7, 2018"
output: 
  html_document: 
    keep_md: yes
---

1. Please Do the following DataCamp activities.
The first two chapters of Writing Functions in R.
The whole course on "Introduction to the Tidyverse"

2. Read in the .csv file (call it rna_counts).

#read .csv file and call it rna_counts


```r
rna_counts <- read.csv("./eXpress_dm_counts.csv", row.names=1)
```

*DONE*: Write a function that calculates and can output mean expression (using the mean() function) for a given data column
. 

*DONE*: Write this function so that it can work on either raw counts (as provided) or transformed to log2 values of the counts, with a logical argument that allows the user to choose whether to do this on the raw scale or log2 scale (i.e. log2 transform, then calculate the mean). 
CHECK at end: Make sure all functions you write for this assignment have this argument. We often use log2 scale for RNAseq data. 
*DONE*: Demonstrate that it works on several of the columns of data.


```r
compute_mean <- function(columnData, isLog2 = FALSE) {
  #column_data = rna_counts[,columnIndex];
  mean_value = NA;
  if(isLog2 ==TRUE){
    log2_values = log2(columnData)
    log2_values <-  log2_values[is.finite(log2(columnData))] #remove -Inf. 
    mean_value = mean(log2_values, na.rm = TRUE)
  }else{
    mean_value <- mean(columnData, na.rm = TRUE);
  }
  return(mean_value)
}
```

#Example of some random numbers to see if the function worked

```r
#compute_mean(rna_counts[,1], isLog2 = FALSE);
compute_mean( rna_counts[,"F101_lg_female_hdhorn"], isLog2 = FALSE)
```

```
## [1] 1978.847
```

```r
compute_mean(rna_counts[,1], isLog2 = TRUE)
```

```
## [1] 9.033518
```

```r
compute_mean(rna_counts[,2], isLog2 = TRUE)
```

```
## [1] 9.045379
```

3. Now using the function you have written, 
Task1 DONE: use a loop to generate a vector of the mean expression value for each column (each sample). 
 
*Task2 DONE*: Make sure that your output vector have named elements (i.e. each element of the vector is associated with the names of the original columns of the data). 

*Task 3 DONE*: Confirm that this function is giving you the correct answer based on what you found in question 2. 

*Task2 analyse*: Do you notice any patterns for which individuals or tissues have the highest mean expression?


```r
compute_matrix_col_mean_loop = function(matrix, isLog2 = FALSE){
  outVector = NULL;
  if(isLog2==FALSE){
  for( col in 1 :ncol(matrix)) {
    key = colnames(matrix)[col]
    value = compute_mean(matrix[,col],isLog2)
    names(value) <- key
    outVector <- c(outVector, value)
  }
  return( outVector)
  }else{
    
    #TRANSFORMED VALUES
  }
}

meanFromLoop <- compute_matrix_col_mean_loop(rna_counts, isLog2 = TRUE)

#compute_matrix_mean_loop(rna_counts, isLog2 = TRUE)
```

4. Repeat this procedure (using the function you wrote, or slightly revised) that uses one of the apply family of functions to do the same as 3. Check which is faster (to compute not to write), and demonstrate how you did this.


```r
compute_matrix_col_mean = function(matrix, isLog2 = FALSE){
     apply(matrix,2, compute_mean, isLog2)
}
meanFromApply <- compute_matrix_col_mean(rna_counts, isLog2 = TRUE)

#check same
identical(meanFromLoop, meanFromApply)
```

```
## [1] FALSE
```

```r
#sort(matrixMean)
#compute_matrix_mean(rna_counts, isLog2 = TRUE)

#check which one is faster' -- loop is faster, apply is slower.
system.time(replicate(1000,compute_matrix_col_mean_loop(rna_counts, isLog2 = TRUE)))
```

```
##    user  system elapsed 
##       0       0       0
```

```r
system.time(replicate(1000,compute_matrix_col_mean(rna_counts, isLog2 = TRUE)))
```

```
##    user  system elapsed 
##   61.64    1.62   75.05
```

```r
#> system.time(replicate(1000,compute_matrix_mean_loop(rna_counts, isLog2 = TRUE)))
#   user  system elapsed 
# 19.23    0.00   19.44
#> system.time(replicate(1000,compute_matrix_mean(rna_counts, isLog2 = TRUE)))
#    user  system elapsed 
# 29.64    0.45   30.14 
```

5. What is a much easier way to do the operations we did in Q 3 and 4, (i.e. you don't need to write your own function) to calculate and output all of the column means? i.e. an Rish way of doing this?

#much way easier way to do means of all the columns

```r
div <- colMeans(rna_counts) 
head(div)
```

```
##  F101_lg_female_hdhorn F101_lg_female_thxhorn   F101_lg_female_wings 
##               1978.847               1983.250               1583.904 
##  F105_lg_female_hdhorn F105_lg_female_thxhorn   F105_lg_female_wings 
##               2105.712               1433.749               1869.962
```

```r
 # The disadvantage is, this one does not have isLog2 option.
```

6. It is common (say for a MAplot) to want the mean expression value of each given gene across all samples. Write a function to do this, and using one of the approaches from Q 3-5 generate and store these values in a variable.

```r
compute_matrix_row_mean = function(matrix, isLog2 = FALSE){
     apply(matrix,1, compute_mean, isLog2)
}

gene_mean_expression = compute_matrix_row_mean(rna_counts)
head(gene_mean_expression) #same gene, accross samples.
```

```
## FBpp0087248 FBpp0293785 FBpp0080383 FBpp0077879 FBpp0311746 FBpp0289081 
##    23.45455  3446.90909    79.54545   139.21818   145.09091  1485.90909
```

7. We are very interested in what is going on in the head horns between small males and large males. Using the type of tools you have written (feel free to modify as you need, but show the new functions) calculate the mean expression for the subset of columns for large and small male head horns. Note you are calculating means on a gene by gene basis, NOT sample by sample. Now calculate the mean difference (again gene by gene) between large male and small males (for head horns). i.e. first calculate the mean expression among individuals who are large males (head horns), ditto for the small males, and calculate their difference.

#grep function to extract sm_male_hdhorn

```r
sm_male_hdhorn <- rna_counts[,grep("sm_male_hdhorn",colnames(rna_counts))]
```

#grep function to extract lg_male_hdhorn

```r
lg_male_hdhorn <- rna_counts[,grep("lg_male_hdhorn",colnames(rna_counts))]
lg_male_hdhorn
```

```
##             M125_lg_male_hdhorn M160_lg_male_hdhorn M180_lg_male_hdhorn
## FBpp0087248                  14                  27                  21
## FBpp0293785                8273                2518                7039
## FBpp0080383                  48                  73                  47
## FBpp0077879                   5                   8                   6
## FBpp0311746                  67                 212                 125
## FBpp0289081                1420                 339                1746
## FBpp0311729                 873                1586                1066
## FBpp0085807                 108                  68                 115
## FBpp0081078                 237                 470                 183
## FBpp0312037                 948                1876                 821
## FBpp0302581                 680                 705                 760
## FBpp0084962                1613                1983                1993
## FBpp0311717                  45                  59                  54
## FBpp0301845                 115                 161                 115
## FBpp0070488                 785                  97                 578
## FBpp0070489                 787                  97                 581
## FBpp0070498                1100                 124                 843
## FBpp0079637                2554                1225                2256
## FBpp0307731                 683                 892                1038
## FBpp0072041                 362                 717                 320
## FBpp0085933                1674                 860                1779
## FBpp0310022                  91                  75                 115
## FBpp0310023                  65                  25                  74
## FBpp0073761                 167                 101                 201
## FBpp0075931                  41                  69                  83
## FBpp0293200                1040                 606                 944
## FBpp0076448                6846                1194                6785
## FBpp0076447                6293                 964                6218
## FBpp0290258                 688                 843                1018
## FBpp0087138                 192                 234                 196
## FBpp0300338                 346                 341                 522
## FBpp0307649                 320                 306                 481
## FBpp0289682                  82                 232                  59
## FBpp0306972                 439                 356                 544
## FBpp0306971                  23                  22                  20
## FBpp0290399                2478                1541                2525
## FBpp0301747                1001                2519                1092
## FBpp0070964                1722                1258                2120
## FBpp0307424                 543                 817                 731
## FBpp0071698                 968                 472                1026
## FBpp0070064                 520                1135                 599
## FBpp0309706                2954                1040                4244
## FBpp0307826                 383                 681                 740
## FBpp0289160                 419                 388                 505
## FBpp0301573                4676                5242                5599
## FBpp0307274                 374                 681                 641
## FBpp0306212                1723                2455                2145
## FBpp0078893                 135                 268                 168
## FBpp0297908                 571                1111                 657
## FBpp0084875                 532                 582                 668
## FBpp0072886                3753                4818                4411
## FBpp0304917                1756                2568                2189
## FBpp0084191                 475                 572                 440
## FBpp0075559                  97                  36                 109
## FBpp0300451                  81                  67                  88
## FBpp0300454                 654                 405                 805
## FBpp0300453                 846                 490                1009
## FBpp0301994                   5                   0                   3
## FBpp0081890                 911                1373                1072
## FBpp0309414                 737                 886                 896
## FBpp0099609                  24                  21                  30
## FBpp0085609                 214                 193                 192
## FBpp0271827                 419                 366                 648
## FBpp0292601                 230                 369                 349
## FBpp0304696                1203                 496                1094
## FBpp0087508                4127                2536                3449
## FBpp0087507                4056                2519                3378
## FBpp0306041                 900                1552                2210
## FBpp0088972                 869                 948                1263
## FBpp0301166                5992                1156                6511
## FBpp0304456                1663                   4                1828
## FBpp0072256                 955                1511                1225
## FBpp0305245                   8                 710                   9
## FBpp0309235                9311               10405                9600
## FBpp0072048                 266                 456                 379
## FBpp0309965                 484                 807                 654
## FBpp0072045                 118                 211                 132
## FBpp0290603                1418                1509                1579
## FBpp0079271                 181                 168                 321
## FBpp0292520                  37                  33                  53
## FBpp0110163                 701                1171                 865
## FBpp0076184                 596                1044                 723
## FBpp0306282                  26                  97                  98
## FBpp0311985                 282                 490                 370
## FBpp0289242                  30                  92                  45
## FBpp0085890                  87                 119                  74
## FBpp0085891                 749                1156                 801
## FBpp0290686                4822                4935                6816
## FBpp0290682                4737                4733                6470
## FBpp0311210                 810                2899                 980
## FBpp0110105                  19                  31                  19
## FBpp0305185                 394                 421                 458
## FBpp0086986                 159                 277                 198
## FBpp0079788                 996                2034                1184
## FBpp0077502                1577                1950                1832
## FBpp0082980                 334                 506                 376
## FBpp0081609                 305                 423                 422
## FBpp0312210                 316                 653                 403
## FBpp0082235                  33                  17                  22
## FBpp0072429                  65                 169                  98
## FBpp0072428                 316                 871                 439
## FBpp0305622                 633                 873                 943
## FBpp0072723                1427                1923                2161
## FBpp0080677                 309                 389                 361
## FBpp0309012                 310                 333                 403
## FBpp0086500                 343                 361                 433
## FBpp0085074                   2                   4                   1
## FBpp0085571                 710                 536                 954
## FBpp0078694                 360                 569                 477
## FBpp0311413                1385                1930                1872
## FBpp0083529                  11                  24                   1
## FBpp0310585                 135                 294                 237
## FBpp0311962                 106                 280                 208
## FBpp0291601                1108                 800                1467
## FBpp0111933                 123                  90                 110
## FBpp0306839                 661                 413                 547
## FBpp0077170                 686                 442                 570
## FBpp0306840                   3                   1                   1
## FBpp0304071                 338                 485                 411
## FBpp0300816                2248                1631                3153
## FBpp0289422                 600                 490                 652
## FBpp0099946                 967                1408                1321
## FBpp0308662                 533                 782                 703
## FBpp0112011                 316                 584                 410
## FBpp0085743                 493                 917                 696
## FBpp0073204                 259                 413                 286
## FBpp0290816                 116                 217                 178
## FBpp0305646                 165                 269                 241
## FBpp0290815                 294                 513                 434
## FBpp0111746                 656                1727                 769
## FBpp0271862                  57                 314                  19
## FBpp0086591                 264                 824                 273
## FBpp0073943                7742                7763                9169
## FBpp0309126                 469                 140                 654
## FBpp0082056                 219                 343                 208
## FBpp0089363                2420                4188                2803
## FBpp0084619                1124                 959                1037
## FBpp0310244                 383                 760                 487
## FBpp0073059                1885                  57                1016
## FBpp0293017                 346                 163                 940
## FBpp0072322                   5                   1                   2
## FBpp0072717                3948                2532                4795
## FBpp0300803                9059                7979               10676
## FBpp0070749                7700                8536                9300
## FBpp0297994                 148                 105                 181
## FBpp0099653                 358                 227                 576
## FBpp0305144                   2                   0                   0
## FBpp0086707                  97                 145                 177
## FBpp0071379                 900                1815                 939
## FBpp0076580                1591                 299                1973
## FBpp0073194                 221                 227                 193
## FBpp0291143                3399                6461                6979
## FBpp0070024                1846                3441                3822
## FBpp0305587                2577                2080                2873
## FBpp0072583                4684                3827                5213
## FBpp0304981                 863                1400                1215
## FBpp0083213                 867                1012                1071
## FBpp0088000                  59                  81                  83
## FBpp0087999                 212                 238                 304
## FBpp0291368                2794                2780                3177
## FBpp0078971                 859                1429                 975
## FBpp0308407                  81                 182                 100
## FBpp0084733                 341                 560                 481
## FBpp0305264                1685                2287                1868
## FBpp0300417                  32                  55                  22
## FBpp0297298                1114                1699                 846
## FBpp0081373                 201                 307                 240
## FBpp0082637                1825                2347                1913
## FBpp0079539                 172                 248                 175
## FBpp0078992                 284                 596                 426
## FBpp0311560                5066                5404                6362
## FBpp0079682                  16                  24                  21
## FBpp0087045                2865                3349                3477
## FBpp0310441                 420                 384                 572
## FBpp0083105                 576                 839                 630
## FBpp0070453                4988                3504                5036
## FBpp0089160                4873                2761                4970
## FBpp0309142                2255                4034                2558
## FBpp0303898                2541                3375                3780
## FBpp0290108                2636                2629                2535
## FBpp0306758                 411                 863                 412
## FBpp0083086                  24                  85                  34
## FBpp0304760                 632                1134                 894
## FBpp0074808                 428                 840                 582
## FBpp0309380                2994                3996                5326
## FBpp0074995                 753                 321                 705
## FBpp0089115                 565                 299                 539
## FBpp0086267                 101                 240                 146
## FBpp0289662                 224                 434                 246
## FBpp0081799                3720                2168                4650
## FBpp0312359                  36                   1                   1
## FBpp0306873                 445                 486                 460
## FBpp0306872                 483                 505                 491
## FBpp0292326                1856                2196                1677
## FBpp0307571                 332                 586                 314
## FBpp0088910                  79                  49                 147
## FBpp0311477                2130                2752                2697
## FBpp0083936                 101                 135                 181
## FBpp0081352                 191                 257                 230
## FBpp0072250                2201                1890                2876
## FBpp0297127                2133                1842                2790
## FBpp0085358                 203                  60                 283
## FBpp0288415                 114                  79                 111
## FBpp0306532                1033                1491                1310
## FBpp0099870                 700                3411                 727
## FBpp0290011                1017                1731                1265
## FBpp0300329                 937                1078                 958
## FBpp0293951                 421                 481                 443
## FBpp0309983                 764                 997                 769
## FBpp0075086                1263                1289                1583
## FBpp0303197                 291                 308                 437
## FBpp0292502                 383                 539                 485
## FBpp0074330                 400                 558                 555
## FBpp0077346                1195                1591                1247
## FBpp0309827                 237                 236                 163
## FBpp0079846                 192                 256                 230
## FBpp0081343                 597                 672                 687
## FBpp0071551                  56                  79                  67
## FBpp0077715                1370                2709                1772
## FBpp0082581                1032                2116                1195
## FBpp0304099                2007                1985                2202
## FBpp0073673                2319                2246                2638
## FBpp0301197                 192                 318                 569
## FBpp0304675                 465                 498                 495
## FBpp0304032                  41                 255                 108
## FBpp0308470                  50                  38                  77
## FBpp0306916                 872                1535                 998
## FBpp0306186                  10                  14                  10
## FBpp0306187                   8                  13                  30
## FBpp0310129                 493                 727                 624
## FBpp0073859                 548                 791                 688
## FBpp0087714                 701                 316                 979
## FBpp0074784                 160                 137                 203
## FBpp0082972                 858                 962                1337
## FBpp0071129                1648                2466                2167
## FBpp0293470                1331                1376                1331
## FBpp0087686                3384                2574                2828
## FBpp0304681                 370                  14                  63
## FBpp0074247                 256                 455                 272
## FBpp0309143                 788                1156                 944
## FBpp0073608                 460                1045                 579
## FBpp0307465                 734                1187                 870
## FBpp0309979                1634                2054                1791
## FBpp0289514                  13                  17                  25
## FBpp0081661                  49                 391                  64
## FBpp0304298                 740                 779                 997
## FBpp0080997                1031                1436                1432
## FBpp0303820                 169                 433                 225
## FBpp0084021                1233                1788                1556
## FBpp0291142                  37                  27                  38
## FBpp0292665                 450                 651                 526
## FBpp0292663                1527                2128                1863
## FBpp0075999                 527                 721                 674
## FBpp0112193                3070                4012                3338
## FBpp0309360                 294                 291                 321
## FBpp0305497                2537                3370                3339
## FBpp0292329                  57                  48                  66
## FBpp0085812                  23                  25                 160
## FBpp0081156                3841                4694                4874
## FBpp0086395                3201                2822                2805
## FBpp0305484                 240                 312                 209
## FBpp0290814                2755                4867                3015
## FBpp0111781                1807                2319                2324
## FBpp0084726                 899                1331                1076
## FBpp0087866                 165                 376                 186
## FBpp0290138                1688                2341                1893
## FBpp0290139                1571                2132                1746
## FBpp0307216                  97                 167                 113
## FBpp0079619                 397                 664                 357
## FBpp0307931                 397                 661                 356
## FBpp0071548                 998                 971                1034
## FBpp0303265                 171                  98                 178
## FBpp0304880                  59                  92                  70
## FBpp0087429                 942                 350                1147
## FBpp0087431                 963                 373                1142
## FBpp0289106                 795                 281                 884
## FBpp0087436                 876                 336                1033
## FBpp0074792                 619                 751                 820
## FBpp0308816                2052                   7                2417
## FBpp0297771                1411                1983                1927
## FBpp0308626                  18                  32                  57
## FBpp0293582                 157                  91                 558
## FBpp0083533                 439                 503                 497
## FBpp0312441                2668                4288                3020
## FBpp0308362                 292                 588                 308
## FBpp0081263                1454                 712                 975
## FBpp0073740                 377                 441                 567
## FBpp0308981                 167                 441                 211
## FBpp0309092                1131                2604                1551
## FBpp0074614                  65                 125                  74
## FBpp0071677                 224                 261                 516
## FBpp0292148                 579                 795                 640
## FBpp0302807                4085                5187                4771
## FBpp0074314                  13                  29                  11
## FBpp0080450                 230                 233                 202
## FBpp0305020                  76                  31                 113
## FBpp0297111                 241                   8                 278
## FBpp0076185                  83                 203                 143
## FBpp0085075                2794                4336                4760
## FBpp0304800                  56                  35                  75
## FBpp0311990                 144                  45                 265
## FBpp0301042                 107                 226                 216
## FBpp0085082                 276                 377                 278
## FBpp0081996                1734                1790                2575
## FBpp0304377               32085                2130               44876
## FBpp0081380                 357                 249                 319
## FBpp0293890                 460                 348                 451
## FBpp0075918                4813               16992                8050
## FBpp0300990                1676                2280                1913
## FBpp0300989                 774                 838                 758
## FBpp0086347                  58                  13                  15
## FBpp0304174                1432                2328                1700
## FBpp0088962                2891                4110                2714
## FBpp0072100                 319                 275                 262
## FBpp0309809                2366                 683                2506
## FBpp0110195                  98                  14                 143
## FBpp0100035                1370                1413                1512
## FBpp0305194                 449                 512                 546
## FBpp0081850                  69                  17                  73
## FBpp0070795                1019                1518                1245
## FBpp0086647                 713                 678                 636
## FBpp0293738                   1                   0                   0
## FBpp0303580                  30                   9                  39
## FBpp0086868                  42                  10                  47
## FBpp0312526                1062                1596                1490
## FBpp0291652                9320               21542               11092
## FBpp0074464                  92                 140                  94
## FBpp0071686                 995                1255                1101
## FBpp0307759               38896               42320               48244
## FBpp0291742                2811                2465                3349
## FBpp0081823               25752               27052               30577
## FBpp0086483                  83                 122                 147
## FBpp0087352                  88                 103                 116
## FBpp0073556                 308                   2                 516
## FBpp0086506                 386                 419                 498
## FBpp0086505                 377                 396                 471
## FBpp0310686                3470                2218                3870
## FBpp0085586               14351               12956               15364
## FBpp0308422                 557                 836                 636
## FBpp0072038                 358                  33                 388
## FBpp0071713                 141                 150                 132
## FBpp0292778                 133                 140                 139
## FBpp0309285                  53                  94                  67
## FBpp0309324                 176                 685                 319
## FBpp0080672               12089               20888               12661
## FBpp0085448                  74                 116                  98
## FBpp0087704                2381                2432                2692
## FBpp0111946                 168                 167                 270
## FBpp0078881                 217                 225                 273
## FBpp0083584                3981                1257                4664
## FBpp0070025                 730                1105                 868
## FBpp0086487               61083                 227               75557
## FBpp0085484               12439               12919               16296
## FBpp0309089                 716                2128                1038
## FBpp0312096                1483                1849                1847
## FBpp0306962                 519                 436                 540
## FBpp0079737                 789                 558                 806
## FBpp0306292                 739                1219                 798
## FBpp0085804                 116                  88                  82
## FBpp0088841                 141                 163                 174
## FBpp0082876                3802                2440                4067
## FBpp0083537                 665                 759                 809
## FBpp0292800                 811                1167                 863
## FBpp0305226                 418                 939                 515
## FBpp0089133                4805                4626                6557
## FBpp0290667                 356                 322                 408
## FBpp0297610                 588                 689                 739
## FBpp0304215                  50                 355                 119
## FBpp0078566                1205                1580                1413
## FBpp0073437                 452                 712                 588
## FBpp0084340                1221                1266                1561
## FBpp0309079                 517                 547                 653
## FBpp0113120                1171                1282                1507
## FBpp0293289                9676                3145               10828
## FBpp0290158                 263                 764                 430
## FBpp0310810                 308                 330                 360
## FBpp0073297                1411                3624                1607
## FBpp0085829                 630                 596                 812
## FBpp0088130                1075                2001                1168
## FBpp0080867                  89                 368                 101
## FBpp0292901                   8                  13                  16
## FBpp0083979                3919                3325                4472
## FBpp0307180                 133                 214                 141
## FBpp0081526                 317                 608                 359
## FBpp0081990                 312                 597                 341
## FBpp0087180                 113                 181                 181
## FBpp0301778                 108                  68                 106
## FBpp0070410                1021                 579                1104
## FBpp0289206                  76                   2                 111
## FBpp0072494                6509                7534                7959
## FBpp0303821                2006                2378                2515
## FBpp0073856                 964                1429                1104
## FBpp0079915                5921                6919                7086
## FBpp0305776                 691                 791                1222
## FBpp0300799                 238                 131                 277
## FBpp0080484                1981                2702                1945
## FBpp0077367                1513                3023                1887
## FBpp0304732                1776                3589                2227
## FBpp0075316                1647                1417                2463
## FBpp0289785                7608                6425                7497
## FBpp0289783                3587                2922                3531
## FBpp0078086                 306                 410                 244
## FBpp0303364                 191                 240                 186
## FBpp0303361                1382                1727                1589
## FBpp0303363                1391                1728                1601
## FBpp0075631                 601                 497                 785
## FBpp0075425                 443                 577                 627
## FBpp0087320                2114                2769                2739
## FBpp0076872                2280                1930                2802
## FBpp0305817                2492                3328                2562
## FBpp0305335                2413                2496                2548
## FBpp0085453                  86                 127                 134
## FBpp0077998                1305                1449                1669
## FBpp0292047                  66                  80                  61
## FBpp0074691                 659                1173                 780
## FBpp0306051                 122                 179                 118
## FBpp0087958                 460                 730                 619
## FBpp0303416                3909                4142                4446
## FBpp0300431                4314                4479                4893
## FBpp0300432                4216                4400                4789
## FBpp0089348                   2                   0                   3
## FBpp0309386                 621                1427                 781
## FBpp0290266                 672                1255                 790
## FBpp0073664                  25                  42                  23
## FBpp0297358                4220                5257                4509
## FBpp0084937                  24                  53                  45
## FBpp0072192                1326                1775                1180
## FBpp0271836                 623                1421                 821
## FBpp0087404                  18                  59                  24
## FBpp0293134                  59                 202                  82
## FBpp0079453                 981                1694                1308
## FBpp0072172                 361                 123                 228
## FBpp0302579                 184                 250                 257
## FBpp0070083                  13                   4                  21
## FBpp0308713                1192                1222                1200
## FBpp0073678                 298                 480                 443
## FBpp0304879               42272               43427               56269
## FBpp0305851                1003                2107                1541
## FBpp0290364                1084                 507                1472
## FBpp0305633                1559                2070                2070
## FBpp0081907                2188                2683                2948
## FBpp0081649                  39                  77                  23
## FBpp0310727                 803                1054                 827
## FBpp0082180                 144                 711                 185
## FBpp0290908                   1                   0                   2
## FBpp0290906                  48                 191                  53
## FBpp0086127                   4                  14                   2
## FBpp0305244                 161                 256                 178
## FBpp0070465                 115                 130                 164
## FBpp0311868                1456                 386                3936
## FBpp0071296                1796                1228                2368
## FBpp0305527                2197                1517                2925
## FBpp0289954                1898                1130                1977
## FBpp0084258                 999                1348                1128
## FBpp0308637                 728                 631                1058
## FBpp0077975                 160                 157                 204
## FBpp0084738                2024                1429                2289
## FBpp0079326                8821                7121                9856
## FBpp0086359                  48                  92                  66
## FBpp0293627                  44                  86                  62
## FBpp0308762                 246                 569                 229
## FBpp0289849                  39                  57                  57
## FBpp0083767                1386                 900                1530
## FBpp0077740                3745                2050                4582
## FBpp0290166                 590                 735                 419
## FBpp0298309                1865                2537                2032
## FBpp0082742                  36                  98                  51
## FBpp0080903                3241               16976                4291
## FBpp0311147               12264                3007               12356
## FBpp0307273                 545                 963                 704
## FBpp0311874                 767                  90                 791
## FBpp0290830                  44                  88                  63
## FBpp0309990                 131                 298                 198
## FBpp0304417               33087               21126               32398
## FBpp0111524                7013                2069                7029
## FBpp0075727                 332                 259                 330
## FBpp0305584                 368                 628                 456
## FBpp0075560                 354                 362                 499
## FBpp0292507                 604                1285                 628
## FBpp0084792                 934                2009                1003
## FBpp0084791                 272                 618                 299
## FBpp0310240               10328                7231               10769
## FBpp0311780                 107                 186                 142
## FBpp0081226                1974                 739                2325
## FBpp0291807                 721                1098                 763
## FBpp0297588                 728                1126                 770
## FBpp0291806                 459                 714                 509
## FBpp0080969                 491                 529                 592
## FBpp0087085                 334                 634                 333
## FBpp0303022                  65                   7                  35
## FBpp0302026                 255                 125                 391
## FBpp0289455                 257                 125                 395
## FBpp0079815                2044                2363                2343
## FBpp0083736                  39                  94                  24
## FBpp0076871                 154                 200                 193
## FBpp0081961                  88                 150                 126
## FBpp0086475                 534                 752                 591
## FBpp0077872                 537                 843                 719
## FBpp0082186                1927                2493                2318
## FBpp0075617                2491                3303                2613
## FBpp0292595                 848                 597                 997
## FBpp0297553                1748                3081                2793
## FBpp0074793                  38                  69                  40
## FBpp0080751                1838                3207                2811
## FBpp0080753                1694                3007                2693
## FBpp0308295                 433                 543                 443
## FBpp0307682                  54                  68                  70
## FBpp0074442                 510                 640                 588
## FBpp0310438                 619                 573                 616
## FBpp0075214                   5                  20                  12
## FBpp0312219                 504                 575                 576
## FBpp0311549                1137                1828                2092
## FBpp0085540                 608                 899                 708
## FBpp0089141               20714               19466               23354
## FBpp0110281                1784               17601                3223
## FBpp0303067                  13                  22                  15
## FBpp0303069                  13                  21                  15
## FBpp0084874                 624                 820                 752
## FBpp0070612                 603                 918                 682
## FBpp0072127                4762                4603                5395
## FBpp0311256                 842                1107                1100
## FBpp0111705                1404                2123                1530
## FBpp0085693                 245                 473                 330
## FBpp0076791                 343                 715                 404
## FBpp0289278                 193                  17                 258
## FBpp0292037                 144                  34                 241
## FBpp0271718                1218                1685                1498
## FBpp0312201                 202                 258                 250
## FBpp0074005                 973                1440                1275
## FBpp0087519                 448                 533                 652
## FBpp0080275                 121                 204                 166
## FBpp0307803                2156                3625                2161
## FBpp0293626                 878                1619                 843
## FBpp0308661                2644                3396                2975
## FBpp0079623                 411                 643                 320
## FBpp0082720                3199                 422                4043
## FBpp0297605                 107                 161                  95
## FBpp0086429                2237                3171                2549
## FBpp0081312                4556                 309                6613
## FBpp0311512                 238                 409                 299
## FBpp0307870                 313                 370                 384
## FBpp0088031                 392                 451                 459
## FBpp0301215                2930                 586                4206
## FBpp0087323                1259                1478                1656
## FBpp0072564                2786                2541                3240
## FBpp0302692                2778                2530                3231
## FBpp0086421                1623                1141                1639
## FBpp0309704                1071                 770                1507
## FBpp0310002                 261                 390                 369
## FBpp0302006                  21                  34                  16
## FBpp0304310                  45                 157                 108
## FBpp0076921                 102                 179                 143
## FBpp0290353                 875                 194                 247
## FBpp0078673                1096                 605                 757
## FBpp0076599                 323                 646                 417
## FBpp0083846                3791                3071                4211
## FBpp0305307                 594                1105                 771
## FBpp0110549                 575                 476                 690
## FBpp0077448                1126                1150                1306
## FBpp0306201                 534                 515                 638
## FBpp0306862                 304                 490                 396
## FBpp0080148                  92                  99                  74
## FBpp0289644                 784                 648                 896
## FBpp0292236                1029                1306                1109
## FBpp0310630                1086                1371                1171
## FBpp0292244                1087                1368                1171
## FBpp0290516                3117                 748                2954
## FBpp0289110                   5                  13                  14
## FBpp0301793                  19                 121                  54
## FBpp0293619                2647                1225                2574
## FBpp0293617                2308                1037                2335
## FBpp0290577                3304                1524                3334
## FBpp0082135                 109                  98                 179
## FBpp0071969                 211                 367                 243
## FBpp0071028                 897                1549                 939
## FBpp0079304                 306                 509                 513
## FBpp0291536                2578                4871                2650
## FBpp0289568                 635                 998                 546
## FBpp0070977                 781                 813                1043
## FBpp0309146                 187                5447                 496
## FBpp0073154                  11                  45                  10
## FBpp0297431                  92                 168                  79
## FBpp0111550                 346                 521                 392
## FBpp0297369                2834                2824                3096
## FBpp0292362                 222                 317                 158
## FBpp0076391                2012                1154                1207
## FBpp0076393                3622                2138                2213
## FBpp0076083                 443                 451                 450
## FBpp0085866                5461                4479                6276
## FBpp0309317                5585                2263                6501
## FBpp0085865                6731                5513                8019
## FBpp0085869                6897                3444                8419
## FBpp0310296                 346                 342                 376
## FBpp0308526                 377                 522                 346
## FBpp0308430                  54                 176                  62
## FBpp0079802                 432                1257                 495
## FBpp0301708                  77                 147                 102
## FBpp0291801                   0                  14                   4
## FBpp0291802                  64                 228                  92
## FBpp0081671                2040                4355                2578
## FBpp0311386                2005                5088                2571
## FBpp0304902                2452                2454                3211
## FBpp0310418                4444                5710                5987
## FBpp0307951                  20                  12                  28
## FBpp0077750                 175                 541                 207
## FBpp0304134                 512                 206                 588
## FBpp0296949                3783                3270                4142
## FBpp0075348                 406                 632                 486
## FBpp0088490                3208                3931                3422
## FBpp0088945                3273                4070                3521
## FBpp0071069                 319                 315                 387
## FBpp0303564               24374                5162               26217
## FBpp0072601                 532                 843                 570
## FBpp0304008                 342                1225                 380
## FBpp0307166                 665                1004                 750
## FBpp0305415                 507                 535                 630
## FBpp0305414                 543                 576                 663
## FBpp0306018                 189                 763                 190
## FBpp0290840                1152                1219                1346
## FBpp0304449                 890                1133                1151
## FBpp0290720                6513                2349                9167
## FBpp0303277                2948                1374                4279
## FBpp0085646                9030                3633               12823
## FBpp0303276                6417                2328                9029
## FBpp0076329                 165                 358                 131
## FBpp0081619                 267                 434                 252
## FBpp0086640                 651                 691                 589
## FBpp0289951                 800                1018                 722
## FBpp0073067                 187                 235                 211
## FBpp0302570                1841                2564                1810
## FBpp0306848                1167                1661                1305
## FBpp0087585                1613                2314                1799
## FBpp0297612               11131                3826                9500
## FBpp0113036               20828                7307               18878
## FBpp0079390                 721                 440                1015
## FBpp0074650                2653                3321                3129
## FBpp0086349                6780                1298                6501
## FBpp0293835                  23                 133                  77
## FBpp0081956                 204                 377                 318
## FBpp0082107                 220                 446                 331
## FBpp0087010                2152                1359                2144
## FBpp0070874                 326                 332                 457
## FBpp0075088                1083                1090                1436
## FBpp0307982              113375              128846               95084
## FBpp0072975                 965                1386                1251
## FBpp0087712                 225                 419                 361
## FBpp0074285                 649                 736                1138
## FBpp0291370                  47                 278                  57
## FBpp0288668               15288                 127               22097
## FBpp0087926                  77                 124                  76
## FBpp0083321                 245                 443                 320
## FBpp0290811                1159                1610                1322
## FBpp0310381                 444                 497                 537
## FBpp0310380                  69                  72                  89
## FBpp0292261                  16                  36                  38
## FBpp0305198                 204                 293                 222
## FBpp0297444                1190                2310                1532
## FBpp0290693                1050                1851                1291
## FBpp0297443                1188                2223                1497
## FBpp0085367                1836                3089                1705
## FBpp0100031                1097                1711                 995
## FBpp0312456                 178                 271                 246
## FBpp0084716                  26                  20                  22
## FBpp0078142                 148                 432                 136
## FBpp0311638                 240                 107                 247
## FBpp0081810                  61                  93                  88
## FBpp0112110                 446                 599                 435
## FBpp0076076                 305                 244                 435
## FBpp0306019                 395                1022                 533
## FBpp0305179                 559                 114                 559
## FBpp0292727                1054                 206                1034
## FBpp0292725                 983                 208                 979
## FBpp0289282                6787                7289                7782
## FBpp0072638                 520                 709                 679
## FBpp0082350                2425                 940                2926
## FBpp0080823                 100                  77                 107
## FBpp0077964                 394                3148                 533
## FBpp0291724                1787                1093                2043
## FBpp0291723                3067                1774                3432
## FBpp0081322                 469                 524                 484
## FBpp0079524                2793                1275                3807
## FBpp0073256                9309               40497               11522
## FBpp0081601               12320               19138               15082
## FBpp0075940                  91                 213                  97
## FBpp0079033                 202                 274                 276
## FBpp0111941                 465                 736                 616
## FBpp0292379                   0                  11                   7
## FBpp0080489                 149                 437                 212
## FBpp0086316                 312                 114                 346
## FBpp0083855                1989                2991                2887
## FBpp0308734                 124                 185                 193
## FBpp0085064                 556                 792                 840
## FBpp0271781                 870                1141                1102
## FBpp0079253                 135                 106                 162
## FBpp0079155                  16                  41                  60
## FBpp0300436                 579                 383                 349
## FBpp0307652                 726                1606                 701
## FBpp0074220                 544                 532                 600
## FBpp0074686                 361                 459                 428
## FBpp0074685                 352                 454                 423
## FBpp0304337                  90                  12                 118
## FBpp0089048                  89                   9                  81
## FBpp0304336                 204                  22                 234
## FBpp0304339                 119                  13                 155
## FBpp0073801                2379                  52                1935
## FBpp0083988                 827                1029                 937
## FBpp0077208                 386                 494                 488
## FBpp0289724                  11                  12                  35
## FBpp0290027                 837                1663                1050
## FBpp0083505                 848                2366                1152
## FBpp0076349                1116                1602                1288
## FBpp0072277                1426                8388                2409
## FBpp0271909                 220                 314                 291
## FBpp0292671                 633                 812                 705
## FBpp0293411                 434                 717                 551
## FBpp0088182                1270                1929                1565
## FBpp0300576                 143                 181                 145
## FBpp0075344                 664                1513                 879
## FBpp0311424                1692                2569                2108
## FBpp0086956                1189                2934                1349
## FBpp0078652                 116                  37                 194
## FBpp0099650                 112                  35                 189
## FBpp0084341                 247                 257                 565
## FBpp0086718                7927                8819               13118
## FBpp0310303                 825                1333                1145
## FBpp0290493                1592                1520                1933
## FBpp0076346                 938                 932                1142
## FBpp0310079                 164                 259                 185
## FBpp0310199                1490                1869                1977
## FBpp0309007                2268                2243                2138
## FBpp0290728                 771                 822                 803
## FBpp0074831                2451                4184                2701
## FBpp0311237                 853                1140                1230
## FBpp0307146                 146                 276                  79
## FBpp0311705                8776                5295               11625
## FBpp0081324                4099                3568                4786
## FBpp0083453                  67                 113                 116
## FBpp0071077                 556                 531                 613
## FBpp0071078                 125                 129                 134
## FBpp0112020                 154                 176                 150
## FBpp0081813                  30                 190                  18
## FBpp0081404                 563                 827                 654
## FBpp0291063                  25                  18                  28
## FBpp0291065                 625                 643                 828
## FBpp0310104                 681                 626                 944
## FBpp0291064                 239                 243                 355
## FBpp0071138                1969                2473                2287
## FBpp0081234                  62                  79                  69
## FBpp0305104                2693                2767                3435
## FBpp0075122                 124                 332                 122
## FBpp0290613                1923                1422                2126
## FBpp0088350                1337                 999                1734
## FBpp0079247                1178                1507                1529
## FBpp0292965                1108                1778                1371
## FBpp0298283                5565                5838                7661
## FBpp0081563                 206                 324                 215
## FBpp0308010                 252                1064                 418
## FBpp0072318                  97                 127                 101
## FBpp0311843                 186                 638                 212
## FBpp0308437                1259                 618                1305
## FBpp0312545                 419                 680                 470
## FBpp0290714                1060                2220                1467
## FBpp0290713                1636                3245                1821
## FBpp0307414                 125                 186                 122
## FBpp0099512                 153                  90                 280
## FBpp0082228                  60                 163                  38
## FBpp0271771                 401                1029                 725
## FBpp0292977                 409                1066                 745
## FBpp0312390                  12                  43                  10
## FBpp0306663                 476                 472                 586
## FBpp0304979                 681                 365                 965
## FBpp0088184                1341                  16                1088
## FBpp0308319               11964                9669               14391
## FBpp0290729                 109                 136                 106
## FBpp0303602                5058                7032                5181
## FBpp0075170                4896                6234                5481
## FBpp0075171                 725                 747                 669
## FBpp0309151                2554                4797                3911
## FBpp0305606                 466                 921                 511
## FBpp0311588                 214                 300                 233
## FBpp0312012                 372                 521                 461
## FBpp0082754                  65                  61                  49
## FBpp0077337                 143                 334                 178
## FBpp0301023                1541                1675                1723
## FBpp0080690                 504                 902                 458
## FBpp0085457                 567                 636                 718
## FBpp0080719                 584                 639                 735
## FBpp0291768                1667                1866                1278
## FBpp0086484                2041                2389                2681
## FBpp0309269                5072                6175                5209
## FBpp0306655                   3                   3                   5
## FBpp0308322                  13                  61                 103
## FBpp0078543                   5                  24                  52
## FBpp0070255                 602                 346                 476
## FBpp0088519                1241                1485                1289
## FBpp0079874                 869                1208                 996
## FBpp0309262                 433                 588                 570
## FBpp0089259                 360                 481                 386
## FBpp0303140                 127                  72                 111
## FBpp0303512                 480                1007                 545
## FBpp0086531                 600                 629                 702
## FBpp0076342                 235                 817                 280
## FBpp0076341                 224                 805                 276
## FBpp0074558                3919                5164                4470
## FBpp0112125                  53                  27                  59
## FBpp0304596                 368                 924                 649
## FBpp0082923                 284                 481                 378
## FBpp0310259                6049                1296                9134
## FBpp0080676                 673                1052                 863
## FBpp0305099                 205                 298                 247
## FBpp0304216                 211                 188                 292
## FBpp0089221                1066                 696                1148
## FBpp0311249                 686                 946                 908
## FBpp0088498                 180                 231                 228
## FBpp0307026                 472                 642                 597
## FBpp0083163                 185                 221                 187
## FBpp0307645                 231                 426                 224
## FBpp0079573                 556                 671                 639
## FBpp0300655               18421               19695               20719
## FBpp0084585                1955                3042                2458
## FBpp0289762                 828                 671                 843
## FBpp0289758                1085                 644                1007
## FBpp0304847                2780                3936                5361
## FBpp0297624                1638                1782                1641
## FBpp0307553                1542                1957                1566
## FBpp0073651                 252                 499                 398
## FBpp0077713                  26                  91                  25
## FBpp0076438                 852                 272                1217
## FBpp0309417                 489                1092                 777
## FBpp0071144                2922                1489                3548
## FBpp0111874                 178                  80                 189
## FBpp0075118                 598                 553                 687
## FBpp0310589                 112                  84                  65
## FBpp0077220                  66                2338                  11
## FBpp0304993              206535                 637              211373
## FBpp0306197                 816                1326                 964
## FBpp0085216                 252                 387                 307
## FBpp0075916                2540                2258                2041
## FBpp0310202                 764                 389                 746
## FBpp0084955                  11                  25                  17
## FBpp0084956                   4                  10                   4
## FBpp0081650                 483                 790                 465
## FBpp0076408                4937                7940                5057
## FBpp0072461                1310                 159                1420
## FBpp0084252                3359                5591                3733
## FBpp0074595                 302                 479                 324
## FBpp0304939                 134                 179                 109
## FBpp0307855                 292                 453                 365
## FBpp0307451                 385                 638                 398
## FBpp0302864                 456                 523                 579
## FBpp0072874                1430                1796                1643
## FBpp0304148                  54                  97                  88
## FBpp0294043                 132                  76                  76
## FBpp0309242                  64                 129                 112
## FBpp0079944                 367                 700                 480
## FBpp0075458                 926                 933                 917
## FBpp0099768                  44                  47                  47
## FBpp0071905                 900                1434                 986
## FBpp0305296                1786                2266                1598
## FBpp0301223                7769                8484                7686
## FBpp0078935                  45                  11                  40
## FBpp0307930                 266                 493                 335
## FBpp0072629                2737                1457                3580
## FBpp0290774                1835                2018                1739
## FBpp0291126                 298                 355                 310
## FBpp0083897                 876                1485                 902
## FBpp0310636                 898                1059                1181
## FBpp0309502                 123                 184                 152
## FBpp0290395                  44                  98                  59
## FBpp0290862                1534                2048                1845
## FBpp0078606                3838                2590                4996
## FBpp0071254                 119                 356                 114
## FBpp0079663               21859               17202               20443
## FBpp0301123                2317                3692                2662
## FBpp0309427                1718                2786                1915
## FBpp0309943                 915                 696                1016
## FBpp0082187                1799                1904                2289
## FBpp0088538                1082                 885                1063
## FBpp0312497                  49                  14                  46
## FBpp0076326                  74                  34                  11
## FBpp0291183                  21                  66                  28
## FBpp0077654                 550                 610                 506
## FBpp0073113                   4                   1                   1
## FBpp0310948                3303                2212                3113
## FBpp0307565                 128                  91                 237
## FBpp0310187                 363                 803                 399
## FBpp0290522                 334                 728                 347
## FBpp0071973                 865                 932                 704
## FBpp0301840                2667                 961                2933
## FBpp0292193                 252                 380                 317
## FBpp0304305                 237                 457                 256
## FBpp0310504                  50                   6                 107
## FBpp0312538                 118                 120                  92
## FBpp0308486                3881                3637                3693
## FBpp0073365                 882                1915                1994
## FBpp0307854                  70                 143                 102
## FBpp0082065                 700                1682                 870
## FBpp0303087                 203                 187                 183
## FBpp0076824                 270                 491                 289
## FBpp0078363                1089                1113                1113
## FBpp0302626               15895               17076               20733
## FBpp0291679                 759                 785                 879
## FBpp0071936                 261                 353                 299
## FBpp0307842                1272                1186                1323
## FBpp0077507                 578                 701                 737
## FBpp0086054                  27                  48                   7
## FBpp0088592                  23                  12                  20
## FBpp0086657                 131                 104                 127
## FBpp0303934                2103                3086                2286
## FBpp0309866                 103                 172                 145
## FBpp0082719                 645                1055                 818
## FBpp0077133                 239                 517                 295
## FBpp0110463                1907                1354                2436
## FBpp0087645                 450                 613                 554
## FBpp0292924                 746                2368                 929
## FBpp0302662                 540                 761                 655
## FBpp0311373                3177                5000                3561
## FBpp0075713                1456                 989                1425
## FBpp0304398                1411                 949                1373
## FBpp0290828                 576                 943                 500
## FBpp0084175                   4                   5                   4
## FBpp0084174                  24                  36                  15
## FBpp0077868                 391                 477                 453
## FBpp0073672                 362                 569                 419
## FBpp0085713                 362                 501                 470
## FBpp0112042                1445                1961                1676
## FBpp0072863                1429                1919                1672
## FBpp0303359                 164                 332                 164
## FBpp0088009                 479                1059                 512
## FBpp0307642                 709                 712                 700
## FBpp0297360                 560                1540                 613
## FBpp0081048                 249                 196                 278
## FBpp0070865                6983                 265                8867
## FBpp0077362                1428                1787                1686
## FBpp0070584                1658                1740                1855
## FBpp0290710                 869                1388                1085
## FBpp0305007                5971                8264                7291
## FBpp0083611                9581                8424               11441
## FBpp0080704                 611                1219                 614
## FBpp0309669                 238                 502                 323
## FBpp0111766                 147                  13                 115
## FBpp0084115                1377                2436                1577
## FBpp0071072                2343                3295                2332
## FBpp0308636                 428                 384                 509
## FBpp0297518                 263                 398                 259
## FBpp0290035                1147                1292                1373
## FBpp0309825                1707                1851                1963
## FBpp0309714                1357                1579                1359
## FBpp0088493                 186                 147                 179
## FBpp0075166                   7                  54                   5
## FBpp0076879                 191                 397                 201
## FBpp0084580                 133                 179                 163
## FBpp0070155                 214                 356                 246
## FBpp0089216                2279                1653                3036
## FBpp0087040                 251                 608                 388
## FBpp0089293                 102                  89                 188
## FBpp0079565                5529                5710                6397
## FBpp0071503                 273                 452                 310
## FBpp0072557                 388                 462                 484
## FBpp0288575                1008                 916                1162
## FBpp0072590                 563                 395                 493
## FBpp0290366                 593                 731                 666
## FBpp0074971                 475                 491                 557
## FBpp0084739                 112                 115                 107
## FBpp0305453                   4                   7                   8
## FBpp0305532                8307                1532               10529
## FBpp0082801                 653                 667                 827
## FBpp0070719                1229                1373                1763
## FBpp0292262                 220                5216                 375
## FBpp0305108                 165                 336                 225
## FBpp0084742                 640                 292                1027
## FBpp0084350                 221                 370                 315
## FBpp0081153              641058              317034              710696
## FBpp0086469                 747                 651                 978
## FBpp0311431                 730                1825                 811
## FBpp0081355                 643                1117                 780
## FBpp0081756                 404                1220                 516
## FBpp0074426                 767                1357                 987
## FBpp0088506                2568                5750                3225
## FBpp0081682                 422                1085                 659
## FBpp0087711                 178                 171                 165
## FBpp0309915                  11                  21                  21
## FBpp0070635                2166                1573                1985
## FBpp0075519                 317                 459                 268
## FBpp0082832                 609                 641                 616
## FBpp0308474                 259                 196                 192
## FBpp0300972                 146                 283                 181
## FBpp0088877                1236                 583                3352
## FBpp0110197                1177                 536                3135
## FBpp0073236                 295                 430                 358
## FBpp0307917                1011                1185                1197
## FBpp0309076                   4                   0                   0
## FBpp0081672               15560                2197               17428
## FBpp0309189                 154                 199                 152
## FBpp0071449                 323                 395                 320
## FBpp0297411                 105                1072                 181
## FBpp0071752                  95                  97                  96
## FBpp0309391                3643                  98                1320
## FBpp0307248                 867                1086                1058
## FBpp0071194                1023                 931                1622
## FBpp0078334                 181                 175                 140
## FBpp0078333                 534                 468                 477
## FBpp0309029                 843                 999                 911
## FBpp0086650                1187                  91                 299
## FBpp0071462                 166                 288                 211
## FBpp0309429                 324                 640                 452
## FBpp0306887                 395                  27                 279
## FBpp0081082                 174                 287                 226
## FBpp0304769                 446                 574                 506
## FBpp0288449                 273                 225                 233
## FBpp0298003                 746                 884                 654
## FBpp0085308                 605                 563                 839
## FBpp0072782                 269                 517                 322
## FBpp0300747                  17                  80                  43
## FBpp0306127                 119                 448                 187
## FBpp0309847                 488                 408                 380
## FBpp0084759                3079                1799                2654
## FBpp0310071                1034                1610                1153
## FBpp0307446                 288                 224                 364
## FBpp0290865                2447                1984                3410
## FBpp0307388                3447                3018                4764
## FBpp0084748                 546                 662                 578
## FBpp0084750                2488                2023                3463
## FBpp0309852                1879                2474                2455
## FBpp0087918                 370                 989                 487
## FBpp0312411                 380                 743                 420
## FBpp0088901                 640                 689                 880
## FBpp0079514                 297                 401                 425
## FBpp0291851                 125                 127                 105
## FBpp0304445                 124                 883                 213
## FBpp0077125                 168                 344                 197
## FBpp0084635                  33                  19                  44
## FBpp0079666                  36                  90                  35
## FBpp0305602                 498                1300                 561
## FBpp0292906                 147                 360                 170
## FBpp0307720                 387                 571                 418
## FBpp0271936                 683                3444                 859
## FBpp0300820                1065                5657                1473
## FBpp0085072                 148                 105                 147
## FBpp0084839                 150                 249                 185
## FBpp0087094                 479                 636                 549
## FBpp0292077                 819                1103                 777
## FBpp0308932                 462                 612                 433
## FBpp0086340                1347                 591                1121
## FBpp0297608                 758                 314                 619
## FBpp0310748                1629                1055                1768
## FBpp0301183                5569                2217                3748
## FBpp0081613                   3                  18                   9
## FBpp0083244                 257                 259                 270
## FBpp0078472                2294                1142                3162
## FBpp0078469                 568                 850                 691
## FBpp0306930                 549                 821                 669
## FBpp0297182                 356                 896                 512
## FBpp0310684                  92                 190                  44
## FBpp0312292                1565                3423                1547
## FBpp0086247                  20                   0                  18
## FBpp0309068                 281                  10                 242
## FBpp0080048                7456                8360                7479
## FBpp0310139                 584                  63                 545
## FBpp0082545                 127                 313                 160
## FBpp0290896                  16                  22                  14
## FBpp0087084                1361                1200                1579
## FBpp0291294                 338                 682                 346
## FBpp0070431                 416                 952                 494
## FBpp0310902                 473                 927                 445
## FBpp0289079                 298                 442                 403
## FBpp0074770                 291                  37                 309
## FBpp0076467                 819                 769                 859
## FBpp0087483                 273                 300                 416
## FBpp0073659                1020                1638                1124
## FBpp0304173                 613                 784                 676
## FBpp0308303                   4                   0                   2
## FBpp0310396                 707                 777                 758
## FBpp0072956                 284                 222                 324
## FBpp0089217                 649                 839                 904
## FBpp0080191                5030                2654                6356
## FBpp0293371                  55                  39                  89
## FBpp0078829                1268                3147                 977
## FBpp0078827                4939               13514                4100
## FBpp0311332                 361                 259                 378
## FBpp0311333                 280                 207                 298
## FBpp0297498                  29                  88                  77
## FBpp0291374                  44                 129                 106
## FBpp0290947                 111                 233                 170
## FBpp0306506                5691               64497               10620
## FBpp0309978                  11                   4                  22
## FBpp0291853                 441                 768                 438
## FBpp0289912                 305                 277                 252
## FBpp0077884                   9                   4                   5
## FBpp0077883                  47                  46                  46
## FBpp0308405                 678                3966                 917
## FBpp0073040                  47                 116                  61
## FBpp0309612                4460                3963                4399
## FBpp0311683                 246                 499                 237
## FBpp0310032                  84                 141                  77
## FBpp0077029                 477                 620                 603
## FBpp0074113                 437                 290                 535
## FBpp0077738                 242                 549                 227
## FBpp0292775                  60                 129                  57
## FBpp0112025               17576                5121               18298
## FBpp0089139               12313                8880               12634
## FBpp0297171                2783                 441                1985
## FBpp0290315                1373                1755                1221
## FBpp0076331                 332                 375                 370
## FBpp0300973                 297                 244                 330
## FBpp0310009                 409                 565                 449
## FBpp0082223                6067                7255                6855
## FBpp0307187                 257                 508                 244
## FBpp0307186                2181                3863                2281
## FBpp0304407                 304                 435                 298
## FBpp0072899                 304                 440                 299
## FBpp0307805                 168                 170                 300
## FBpp0084244                 457                 951                 537
## FBpp0112160                 223                 337                 254
## FBpp0308524                 459                 573                 571
## FBpp0072453                 557                 604                 571
## FBpp0073902                6002                6759                6781
## FBpp0300790                1573                2539                2022
## FBpp0304832                 244                 587                 300
## FBpp0294034                 248                 615                 308
## FBpp0307426                1086                4059                1845
## FBpp0081896                1188                4513                2046
## FBpp0084486                 133                 227                 149
## FBpp0080495                1759                1690                1989
## FBpp0307142                1697                1670                1922
## FBpp0080496                 554                 563                 648
## FBpp0079454                4532                5188                4891
## FBpp0075396                 239                 225                 273
## FBpp0077873                 119                 140                 155
## FBpp0075398                 637                1287                 878
## FBpp0087732                 100                 150                 106
## FBpp0306672                 972                1399                1344
## FBpp0306632                 902                3547                1290
## FBpp0083179                1385                5454                1881
## FBpp0081055                 425                 567                 569
## FBpp0290679                 570                 734                 759
## FBpp0309338                 655                1176                 874
## FBpp0087780                1579                3195                1490
## FBpp0111956                  91                 133                 108
## FBpp0308815                 821                 954                 980
## FBpp0304976                1369                1910                1686
## FBpp0111689                 859                1611                 846
## FBpp0305859                   7                  17                  30
## FBpp0074836                  12                  62                  53
## FBpp0289451                 354                 311                 369
## FBpp0083935                 187                 426                 222
## FBpp0075248                3631                4094                3741
## FBpp0305488                3710                4126                3807
## FBpp0076488                 626                1060                 671
## FBpp0111729                 616                1046                 656
## FBpp0076486                 579                 943                 626
## FBpp0297602                 613                1041                 644
## FBpp0076413                 541                 912                 600
## FBpp0292356                1089                 192                1158
## FBpp0301966                 296                  33                 345
## FBpp0080392                 270                 293                 344
## FBpp0086590                 690                 850                 864
## FBpp0310569                 237                 373                 236
## FBpp0083225                  36                 101                  68
## FBpp0082624                 107                 149                 140
## FBpp0311839                 693                1574                 972
## FBpp0078425                 200                 230                 273
## FBpp0309629                   2                   2                   2
## FBpp0307200                 462                 665                 707
## FBpp0082070                 107                 174                 179
## FBpp0084348                 795                 521                 857
## FBpp0110299                 120                  41                 107
## FBpp0074162                1930                 787                1949
## FBpp0297167                 381                 634                 500
## FBpp0078402                 151                 193                 173
## FBpp0080017                 579                 545                 721
## FBpp0082352                1676                1848                2232
## FBpp0085616                 689                 983                 795
## FBpp0070373                 377                 981                 416
## FBpp0309962                 307                 763                 406
## FBpp0297204                1583                2232                1761
## FBpp0311440                  37                 134                  25
## FBpp0289046                  38                 137                  25
## FBpp0070094                 122                  27                 175
## FBpp0306551                1554                2048                1672
## FBpp0311451                2401                3084                2579
## FBpp0306904                 162                1731                 325
## FBpp0292890                 167                1919                 348
## FBpp0079182                 266                 409                 293
## FBpp0081434                 354                 407                 416
## FBpp0310065                 440                1543                 499
## FBpp0290409                6979                  23                5890
## FBpp0305003                 290                 267                 253
## FBpp0290541                  89                 100                  94
## FBpp0110215                 235                 594                 282
## FBpp0110216                 231                 612                 273
## FBpp0289622                 205                 334                 252
## FBpp0311786               79774                  37              104226
## FBpp0288892                1178                1116                1442
## FBpp0080237                1644                3475                1745
## FBpp0309425                 465                1022                 518
## FBpp0079519                 396                1579                 553
## FBpp0073717                 594                 728                 709
## FBpp0080293                  39                 192                   8
## FBpp0082359                 140                 294                 213
## FBpp0074713                  22                  28                  41
## FBpp0294007                  84                 141                 105
## FBpp0075937                 918                1285                1076
## FBpp0302564                1499                1281                2875
## FBpp0310151                  64                 107                  77
## FBpp0307613                2409                1095                2458
## FBpp0307614                3704                1718                3881
## FBpp0305795                1774                1933                2209
## FBpp0305929                  10                  10                  28
## FBpp0303428                1917               11576                2437
## FBpp0077343                 822                1389                 646
## FBpp0293153                 356                 350                 388
## FBpp0308330                 155                 151                 146
## FBpp0074538                 353                 347                 387
## FBpp0074246                 155                 415                 207
## FBpp0289116                 349                 382                 468
## FBpp0078927                1021                1143                1516
## FBpp0111580                1140                 306                1139
## FBpp0111581                 519                 102                 432
## FBpp0311993                  77                 157                 108
## FBpp0301283                 384                 502                 478
## FBpp0310514                 809                1063                1049
## FBpp0083814                  38                  23                  48
## FBpp0309252                 846                1012                1084
## FBpp0309251                  19                  15                  17
## FBpp0085546                 372                 481                 421
## FBpp0070899                2675                4661                2982
## FBpp0308555                  35                  54                  31
## FBpp0073199               26409                2621               25308
## FBpp0085617                 204                 407                 236
## FBpp0083879                 314                 903                 368
## FBpp0076533                 204                 231                 264
## FBpp0072033                 673                 662                 733
## FBpp0300997                 434                  19                 512
## FBpp0290635                  87                   0                  79
## FBpp0290637                  21                   0                  25
## FBpp0307448                 106                 173                 113
## FBpp0305528                 107                  87                  76
## FBpp0309676                 741                 491                 673
## FBpp0289404                 315                 449                 518
## FBpp0082120                  33                  70                  12
## FBpp0075441                 558                 768                 780
## FBpp0081493                 405                 555                 431
## FBpp0311646                  30                 116                  43
## FBpp0297092                  41                 182                 134
## FBpp0079847                 165                 232                 205
## FBpp0078372                 613                 773                 758
## FBpp0311009                1354                2035                1391
## FBpp0312439                 166                 140                 138
## FBpp0085365                1169                  44                1096
## FBpp0082597                 477                  15                 424
## FBpp0291071                8472                2599                7406
## FBpp0084210                  31                  58                  41
## FBpp0073103                 446                 876                 615
## FBpp0076132                 278                 337                 257
## FBpp0075935                 829                1584                 854
## FBpp0075934                 603                1033                 540
## FBpp0310637               16203                2718               17730
## FBpp0079527                3053                8901                4079
## FBpp0076784                 580                 876                 653
## FBpp0113091                 568                 833                 594
## FBpp0290169                 210                 185                 234
## FBpp0081886                1094                2579                1370
## FBpp0087154                 197                 380                 362
## FBpp0300236                1575                 371                 943
## FBpp0083713                 812                1712                 920
## FBpp0084689                1167                1431                1547
## FBpp0081473                1289                2865                1431
## FBpp0305823                1222                2828                1353
## FBpp0074206                 314                 523                 466
## FBpp0074205                 296                 497                 436
## FBpp0083413                3195                4076                3953
## FBpp0311406                1647                  50                1675
## FBpp0111853                 109                 164                 141
## FBpp0071938                 615                 688                 784
## FBpp0071937                 730                 806                 910
## FBpp0085912                 370                 402                 624
## FBpp0311802                 376                 409                 634
## FBpp0089422                 979                1903                1241
## FBpp0110316                 824                1560                1036
## FBpp0293227                2300                3478                2780
## FBpp0290866                 144                  17                 541
## FBpp0292613                1206                1533                1435
## FBpp0292910                3047                2962                2901
## FBpp0301134                   0                   0                   0
## FBpp0288437                   6                   5                   4
## FBpp0312193                 118                 120                 117
## FBpp0076936                 137                 187                 105
## FBpp0305555                 255                 302                 234
## FBpp0301128                 188                 419                 172
## FBpp0301125                  72                 193                  67
## FBpp0300161                1608                2292                1884
## FBpp0084452                1801                2634                2096
## FBpp0306270                1468                2067                2080
## FBpp0307572                 197                 118                  97
## FBpp0303191                 631                1086                 740
## FBpp0077893                 314                 375                 410
## FBpp0083885                 161                  50                 136
## FBpp0290185                 287                  78                 232
## FBpp0080931                4367                 524                7024
## FBpp0078296                 326                 890                 403
## FBpp0311643                5533               11182                7520
## FBpp0076186                  39                  26                  26
## FBpp0110083                   7                   0                   5
## FBpp0297227                2201                1226                3040
## FBpp0311995                2253                3907                3091
## FBpp0312093                2798                4501                3235
## FBpp0073588                 832                1165                 950
## FBpp0304873                 175                 237                 236
## FBpp0302777                3963                6748                4374
## FBpp0070330               10086               12555               11658
## FBpp0289493                 662                3604                 706
## FBpp0290057                 126                 195                 158
## FBpp0311034                1193                1266                 469
## FBpp0071746                 359                 347                 136
## FBpp0305326                 470                1326                 806
## FBpp0311768                 768                  56                1133
## FBpp0311803                 185                 373                 197
## FBpp0271832                 396                 587                 436
## FBpp0089163                 747                1263                 786
## FBpp0305662                 540                 848                 628
## FBpp0305522                 164                 286                 220
## FBpp0290099                  18                  33                  23
## FBpp0073128                 113                 190                 132
## FBpp0070110                  71                  86                  57
## FBpp0290755               11751                5889               12450
## FBpp0078678               16814                8006               20344
## FBpp0078679               11460                5657               12192
## FBpp0078680               11471                5632               12136
## FBpp0290753               15650                9108               15212
## FBpp0078674               15755                9149               15351
## FBpp0311450               43735                3550               35027
## FBpp0083227                1547                1656                1522
## FBpp0070652                1646                2067                1799
## FBpp0087744                 377                1981                1714
## FBpp0083842                4345                 412                4090
## FBpp0310511                  84                  87                 118
## FBpp0074671                 990                1672                1321
## FBpp0297103                1255                1924                1714
## FBpp0307144                1267                1812                1752
## FBpp0082314                1390                1363                1596
## FBpp0080747                 307                 473                 414
## FBpp0307178                 824                 981                 946
## FBpp0079254                 118                 118                 126
## FBpp0079255                 450                 535                 568
## FBpp0078797                  10                  17                  12
## FBpp0312226               16008               13426               12985
## FBpp0076718                 295                 494                 457
## FBpp0311478                 818                1452                1051
## FBpp0070751                 594                1094                 825
## FBpp0301702                 642                1074                 676
## FBpp0071288                1097                 829                1060
## FBpp0084410                1092                1137                1458
## FBpp0290701                1735                2305                2059
## FBpp0304333                  49                  33                  55
## FBpp0304672                 136                 540                 157
## FBpp0312207                 617                 573                 694
## FBpp0078116                 217                 363                 251
## FBpp0081608                 944                1497                1126
## FBpp0303579                   0                   4                  10
## FBpp0070730                 121                 239                 138
## FBpp0086858                 101                   9                 119
## FBpp0100135                   7                  18                   8
## FBpp0304159                   2                   7                   2
## FBpp0079425                2333                3038                2528
## FBpp0072612                 805                 957                 699
## FBpp0083514                  56                 115                  71
## FBpp0083658                 664                 954                 848
## FBpp0300948                 446                1221                 449
## FBpp0309228                4057                3128                5517
## FBpp0077899                 112                 656                 160
## FBpp0309275                 154                  91                 239
## FBpp0099795                 527                 777                 744
## FBpp0088033                 111                 227                 124
## FBpp0292920                 589                 608                 590
## FBpp0079347                 513                1311                 532
## FBpp0309192                 132                 217                 251
## FBpp0307570                 552                 770                 604
## FBpp0306011                 859                1039                 781
## FBpp0080118                1630                1804                1737
## FBpp0310753                1751                 882                1149
## FBpp0081780                4187                6872                4995
## FBpp0311199                  66                  83                  82
## FBpp0303185                 114                 203                 161
## FBpp0297586                5076                9619                4591
## FBpp0290579                  75                  78                  36
## FBpp0293283                 416                 763                 796
## FBpp0083436                 639                 836                 758
## FBpp0083399                 330                 434                 379
## FBpp0075276                 671                 891                1001
## FBpp0072803                 525                 573                 595
## FBpp0293860                 330                 313                 432
## FBpp0083275                  52                   9                 142
## FBpp0310640                  96                  46                  72
## FBpp0305834                 215                  50                 290
## FBpp0072562                1939                1866                2073
## FBpp0078726                 156                 381                 114
## FBpp0293533                 634                 611                 539
## FBpp0293530                 525                 505                 467
## FBpp0072906                 661                1137                1087
## FBpp0111528                1241                2454                1229
## FBpp0304182                1128                1397                1463
## FBpp0074559                1169                 513                1808
## FBpp0309263                1340                 351                1501
## FBpp0302955                 254                1075                 299
## FBpp0304613                2746                3024                1972
## FBpp0077447                 514                 663                 602
## FBpp0085763                 394                 856                 425
## FBpp0090954                 441                 592                 547
## FBpp0306434                 302                 499                 378
## FBpp0089292                  33                 123                  35
## FBpp0310440                1126                 848                1024
## FBpp0303715                  32                  48                  14
## FBpp0309362                 770                1222                 772
## FBpp0305199                  69                  65                  56
## FBpp0076330                  47                 183                  59
## FBpp0288869                 851                1682                 990
## FBpp0311911                1529                1515                1746
## FBpp0307983                7175                8458                9166
## FBpp0302536                 947                1016                 877
## FBpp0306001                 249                 490                 255
## FBpp0307977                   1                   3                   0
## FBpp0305775                7669                3237                6766
## FBpp0297531                   2                   4                   0
## FBpp0301158                  33                  52                  28
## FBpp0089033                4993                8673                5720
## FBpp0304835                5020                8727                5748
## FBpp0073620                1423                 989                1697
## FBpp0089405                 126                 114                 121
## FBpp0088792                1641                1306                1536
## FBpp0292532                 392                 580                 373
## FBpp0307163                 742                1100                1027
## FBpp0080012                 853                 342                 853
## FBpp0312437                3267                1229                3740
## FBpp0078642               28871               53364               43171
## FBpp0085657                  39                 568                  54
## FBpp0293054                 390                 600                 464
## FBpp0311547                2030                4313                2699
## FBpp0302962                 147                 293                 218
## FBpp0076520                1661                2385                2154
## FBpp0073714                 468                 736                 814
## FBpp0070708                1145                6134                1970
## FBpp0070707                  65                 277                  95
## FBpp0310379                 503                 728                 571
## FBpp0086501                 858                 531                1097
## FBpp0074343                1106                 963                1136
## FBpp0074342                1176                1158                1207
## FBpp0311662                  10                  35                  15
## FBpp0291423                 158                 131                 284
## FBpp0087398                 392                 630                 460
## FBpp0305083                5654                  19                6040
## FBpp0310796                3312                  15                2800
## FBpp0080945                 873                 604                1565
## FBpp0072475                 463                 648                 567
## FBpp0070798               14952               57888               18428
## FBpp0085809                2159                3322                3013
## FBpp0084176                5087                4224                5439
## FBpp0305448                1725                1219                2011
## FBpp0306905                1576                2248                1707
## FBpp0112097                  84                  82                  47
## FBpp0311408                 356                   7                 310
## FBpp0297108                 116                  35                 144
## FBpp0311758                 295                 655                 433
## FBpp0305422                  54                  88                  71
## FBpp0075680                1207                2911                1470
## FBpp0087544                 326                1210                 316
## FBpp0290074                 621                2209                 578
## FBpp0087535                 698                 662                 699
## FBpp0081130                 207                 238                 173
## FBpp0289747                1581                1647                1843
## FBpp0309398                 505                 717                 718
## FBpp0309396                 172                 255                 224
## FBpp0288679                 617                 625                 781
## FBpp0302956                 382                 808                 501
## FBpp0271928                 102                 190                  96
## FBpp0073823                5247                5922                5263
## FBpp0079924                1034                2178                1725
## FBpp0087855                 166                 225                 189
## FBpp0077925                1010                3429                1086
## FBpp0080817                 550                 839                 763
## FBpp0311225                 556                 677                 699
## FBpp0071217                1475                1237                1698
## FBpp0289916                   0                  10                   3
## FBpp0070006                1101                1314                1362
## FBpp0070355                4597                4301                5501
## FBpp0304913                  27                   5                  82
## FBpp0083613                 124                 121                 117
## FBpp0306504                 555                 877                 619
## FBpp0077119                 624                1001                 736
## FBpp0311657                5186                6817                4643
## FBpp0309217                 273                 486                 314
## FBpp0075794                3767                4586                5187
## FBpp0084386                  24                   0                  23
## FBpp0074401                1306                3443                2396
## FBpp0293060                1086                1429                1419
## FBpp0297302                2933                3342                2542
## FBpp0289196               82382                 969               66429
## FBpp0307245               48122                 658               37124
## FBpp0289199               81566                 996               65642
## FBpp0290285                3677                 428                3298
## FBpp0070041                2946                3564                3621
## FBpp0082985                3141                3950                4562
## FBpp0293233                 527                 820                 631
## FBpp0070811                 487                 875                 614
## FBpp0077867                 317                 542                 407
## FBpp0301152               12255               12487               15684
## FBpp0307374                1151                1151                1247
## FBpp0088044                 525                 533                 542
## FBpp0310230                 872                1697                 962
## FBpp0297726                 898                 316                 644
## FBpp0312446               10899                3340                9909
## FBpp0074915               10977                3552               10058
## FBpp0074917               11469                3124               10217
## FBpp0086664                 378                 347                 419
## FBpp0308325                 450                 653                 538
## FBpp0099896                 574                1523                 775
## FBpp0291730                7459               19840               10595
## FBpp0291729                 176                 363                 202
## FBpp0111759                 302                 216                 378
## FBpp0311113                 451                 341                 537
## FBpp0072567                2888                 682                3671
## FBpp0072568                2769                 647                3433
## FBpp0072685                   0                   1                   1
## FBpp0083686                 425                 520                 533
## FBpp0290295                 635                 791                 719
## FBpp0309318                2173                4979                2245
## FBpp0088429                 425                1037                 419
## FBpp0291866                 217                1247                 284
## FBpp0290196                 715                 702                 973
## FBpp0290195                  66                  57                  80
## FBpp0073247                 104                   3                  24
## FBpp0304109                  25                   1                  10
## FBpp0073246                 160                   4                  37
## FBpp0079042                  17                  20                  17
## FBpp0303230                 293                 377                 240
## FBpp0082592                 784                1017                 978
## FBpp0308326                1754                1909                1360
## FBpp0298033                 548                 527                 519
## FBpp0293201                3419                5793                3387
## FBpp0303033                 232                 462                 315
## FBpp0071049                 531                 468                 664
## FBpp0112712                 944                1462                1034
## FBpp0112713                 626                 958                 680
## FBpp0304252                2570                2783                3389
## FBpp0076892                2656                2377                2942
## FBpp0073893                 648                1161                 597
## FBpp0304930                 321                 143                 383
## FBpp0304921                2242                 916                2023
## FBpp0304920                3517                1404                3343
## FBpp0304314                 298                 463                 343
## FBpp0310097                 233                 346                 283
## FBpp0075360                 652                1317                 602
## FBpp0080045                 262                 496                 382
## FBpp0082055                  75                 180                 103
## FBpp0080915                 233                 470                 321
## FBpp0086116                1855                 987                1745
## FBpp0310439                  38                  84                  45
## FBpp0303068                 498                 521                 586
## FBpp0297935                2155                2137                2367
## FBpp0307449                 116                 281                 149
## FBpp0082935                 134                 173                 151
## FBpp0112331                1295                1840                1617
## FBpp0302987                 188                  70                 184
## FBpp0310714                 516                 748                 618
## FBpp0083082               28383                1144               49013
## FBpp0076385               46500              224184               61522
## FBpp0271876                 180                 354                 240
## FBpp0312324                 291                 229                 418
## FBpp0292158                1082                 683                1736
## FBpp0302934                1227                1413                1546
## FBpp0307685                  85                 198                 129
## FBpp0300173                 217                 370                 227
## FBpp0088654                 209                 330                 209
## FBpp0077520                 159                 228                 197
## FBpp0075038                 545                 678                 658
## FBpp0077503                1566                2799                1707
## FBpp0305130                1424                1964                1612
## FBpp0298330                 598                 204                 666
## FBpp0298328                 284                  96                 304
## FBpp0082036                  77                 182                  99
## FBpp0082034                  23                  32                  22
## FBpp0303610                  80                 191                 101
## FBpp0082035                  14                  65                  28
## FBpp0297342                1762                2317                1894
## FBpp0297339                1851                2453                2011
## FBpp0087069                  36                 115                  84
## FBpp0293018                 204                 337                 160
## FBpp0311368                  86                  60                  52
## FBpp0083958                 185                 317                 128
## FBpp0300568                1246                1062                1464
## FBpp0071192                1650               10248                2290
## FBpp0087058                1220                8343                1634
## FBpp0298037                   1                   7                   5
## FBpp0309260                 483                 562                 545
## FBpp0072148                 134                 155                 192
## FBpp0072952                 724                 300                1240
## FBpp0079770                 425                 657                 574
## FBpp0305544                  98                  85                  57
## FBpp0087353                 910                1380                1139
## FBpp0305265                4935                8590                4539
## FBpp0293000                1461                2747                1400
## FBpp0307992                2722                2739                3302
## FBpp0292404                1921                2017                1777
## FBpp0293465                1075                1333                1294
## FBpp0307387                4908                2800                5049
## FBpp0087120                 213                 282                 229
## FBpp0297105                3042                1896                2737
## FBpp0297107                  17                  53                 179
## FBpp0304785                1188                1787                1549
## FBpp0080319                 783                 877                 685
## FBpp0304830                5113                5965                7285
## FBpp0085870                1255                2099                1578
## FBpp0302018                3246                 716                2592
## FBpp0312568                 889                1708                1006
## FBpp0073848                 198                 403                 252
## FBpp0305206                3880                 646                3399
## FBpp0310741                3960                2029                3432
## FBpp0292340                 198                 218                 174
## FBpp0078920                 430                 646                 486
## FBpp0307630                1221                 384                1502
## FBpp0310809                 469                 327                 604
## FBpp0292773                 256                 491                 185
## FBpp0307564                 377                 395                 303
## FBpp0305270                 514                 789                 586
## FBpp0072006                  94                  97                 132
## FBpp0079156                  94                  97                 160
## FBpp0080862                 257                 197                 304
## FBpp0070216                1200                1703                1784
## FBpp0309594                 375                 746                 549
## FBpp0076654                1560                1547                2025
## FBpp0304343               12592               16072               16301
## FBpp0308578                  63                 209                  86
## FBpp0293051               64641               14757               69439
## FBpp0078722                 287                 306                 186
## FBpp0078723                 223                 261                 141
## FBpp0290848                 786                 710                1089
## FBpp0290847                  80                 176                 107
## FBpp0081509                2216                3325                2634
## FBpp0111835                2099                 251                3180
## FBpp0078268               86847                8880              110370
## FBpp0073996                1407                2138                1240
## FBpp0311495                1852                1227                1751
## FBpp0073674                1225                1451                1391
## FBpp0087406                  30                  56                  90
## FBpp0305988                 338                 246                 274
## FBpp0071406                 314                 567                 490
## FBpp0312156                 209                 561                 320
## FBpp0291523                 193                 319                 228
## FBpp0300977                  64                  93                 112
## FBpp0084089                1285                3172                1463
## FBpp0085117                 167                 208                 219
## FBpp0111932                 300                 438                 425
## FBpp0077316                 581                 808                 848
## FBpp0288548                1645                2785                1170
## FBpp0288731                 819                1694                 600
## FBpp0078648                  93                 273                 145
## FBpp0297777                 107                 176                 139
## FBpp0303837                 162                 158                 151
## FBpp0076549                 940                 777                1050
## FBpp0300839                   2                   2                   0
## FBpp0099631                 162                 362                 190
## FBpp0306416                 337                 685                 398
## FBpp0099632                 309                 635                 374
## FBpp0301556                 330                 669                 391
## FBpp0075675                  92                 106                  83
## FBpp0081447                  47                 167                  78
## FBpp0082768                  19                  20                   7
## FBpp0300995                  21                  33                  41
## FBpp0078182                 842                 358                1011
## FBpp0288833                1592                 440                1253
## FBpp0080280               17671               24632               20301
## FBpp0077262                4660                2727                5403
## FBpp0290850                  21                  12                  12
## FBpp0300975                  47                  42                 129
## FBpp0308852                1916                3043                1870
## FBpp0306093                 621                2044                 654
## FBpp0086732                  48                 104                  53
## FBpp0300976                1072                 867                1067
## FBpp0110523                 699                 986                 797
## FBpp0110483                  98                 115                  91
## FBpp0073449                 403               10046                1131
## FBpp0071298                 207                 537                 233
## FBpp0075693                 934                1010                1396
## FBpp0305742                 108                 227                 165
## FBpp0075239                 759                1117                1087
## FBpp0085765                 143                  95                  67
## FBpp0308364                 229                 194                 224
## FBpp0071279                5293                9563                6629
## FBpp0301591                  15                  19                  32
## FBpp0088357                 286                 356                 199
## FBpp0309307                 400                 545                 477
## FBpp0074026                1443                 669                1556
## FBpp0088013                 310                 589                 397
## FBpp0290003                 317                 614                 404
## FBpp0297426                 679                 584                 686
## FBpp0071700                 213                 254                 260
## FBpp0070224                 371                 631                 418
## FBpp0074788                1649                1172                1904
## FBpp0311541                1942                1667                2406
## FBpp0082664                 996                1790                1046
## FBpp0302792                 523                1074                 590
## FBpp0305673                 702                 942                 965
## FBpp0071115                4562                4939                5066
## FBpp0078166                4267                4616                4735
## FBpp0309959                 165                 587                 246
## FBpp0070831                 283                1244                 534
## FBpp0071418                2058                2556                1976
## FBpp0300807                2054                2576                1983
## FBpp0075115                 583                1055                 623
## FBpp0079819                 941                1310                1020
## FBpp0082550                 837                1100                1132
## FBpp0297136                   8                  13                  15
## FBpp0305313                9178                2052                9190
## FBpp0086941                4473                1034                4658
## FBpp0311473                 954                1049                1564
## FBpp0079685                 237                 483                 159
## FBpp0071212               12740               12528               12740
## FBpp0297401                 816                3173                3591
## FBpp0309997                1237                1521                1361
## FBpp0085715                 129                1131                 215
## FBpp0293217                 179                1506                 274
## FBpp0074936                 435                 694                 479
## FBpp0084907                4502                6487                6369
## FBpp0088054                 127                 235                 198
## FBpp0083932                 396                 945                 564
## FBpp0070876                  68                 341                 119
## FBpp0291029                   3                   1                   0
## FBpp0077302                 236                 305                 246
## FBpp0309145                 204                 261                 220
## FBpp0078358                 384                 495                 308
## FBpp0088561                1798                3571                1952
## FBpp0293844                 750                 841                 694
## FBpp0293837                 410                 434                 339
## FBpp0308484                 116                 419                 148
## FBpp0309944                 618                 912                 602
## FBpp0077022                 145                 211                 191
## FBpp0077099                 204                 297                 225
## FBpp0088818                 370                 809                 618
## FBpp0084213                  97                 171                  94
## FBpp0311052                 849                1858                 908
## FBpp0291660                 779                 989                 934
## FBpp0309024                 870                1373                1102
## FBpp0289975                2329                 478                2060
## FBpp0082989               12313               20078               13401
## FBpp0080125                1037                 745                1026
## FBpp0080124                1057                 785                1080
## FBpp0084043                 112                 163                 135
## FBpp0075096                 766                1091                 913
## FBpp0074104                1012                1154                1376
## FBpp0082541                2197                2256                2514
## FBpp0079549                 592                1707                 991
## FBpp0072239                1003                1775                1123
## FBpp0290893               54966               63328               46309
## FBpp0304044                6853                6250                6853
## FBpp0289202                 189                 316                 228
## FBpp0307782                 118                 183                 186
## FBpp0305102                1201                2979                1445
## FBpp0085589                 276                 275                 350
## FBpp0310242                1404                1646                1803
## FBpp0071849                 275                 666                 466
## FBpp0071269                 787                1017                 963
## FBpp0078843                 831                 608                1197
## FBpp0307972                 640                 520                 809
## FBpp0082855                 206                 178                 260
## FBpp0076485                 733                1036                 522
## FBpp0300612               47398               89526               64081
## FBpp0074568                1527                2357                1467
## FBpp0309288                 766                1705                1114
## FBpp0082148                 449                 421                 467
## FBpp0303793                3080                3650                3298
## FBpp0086939                 651                 931                 868
## FBpp0082813                1048                 899                1057
## FBpp0084852                6174              278992               14601
## FBpp0084729                 360                 502                 519
## FBpp0075186                 427                 451                 651
## FBpp0306920                 765                1032                 990
## FBpp0303568                2143                 787                2026
## FBpp0303566                2813                1053                2666
## FBpp0303567                2362                 882                2204
## FBpp0303572                 380                 133                 352
## FBpp0079844                  42                  50                  49
## FBpp0304186                 237                 365                 248
## FBpp0079699                 553                 546                 716
## FBpp0311506                 211                  88                 193
## FBpp0303939                 857                 399                 850
## FBpp0297063                2163                 249                2698
## FBpp0099825                1661                 215                2009
## FBpp0309030                7115                8314                7227
## FBpp0289825                 745                 425                 615
## FBpp0310679                 141                 179                 156
## FBpp0290777                 163                 299                 248
## FBpp0306931                3611                3268                2606
## FBpp0306932                4031                3391                2852
## FBpp0111703                9796                9378               10595
## FBpp0079621                 543                3119                 960
## FBpp0079620                 506                2875                 882
## FBpp0099403                 324                 533                 338
## FBpp0305725                 351                 310                 307
## FBpp0294008                 431                1116                 429
## FBpp0073262                 789                1376                 991
## FBpp0312024                1091                1547                1352
## FBpp0305599                1085                1519                1333
## FBpp0308632                  80                 172                  93
## FBpp0310095                1012                 718                1125
## FBpp0081255                  88                 736                 157
## FBpp0311801                 257                 337                 276
## FBpp0088552                 461                 350                 469
## FBpp0088549                 134                  97                 117
## FBpp0088550                 422                 336                 458
## FBpp0088551                  61                  48                  61
## FBpp0311398               20162                2296               18714
## FBpp0297162                4143                 471                3597
## FBpp0289985                 133                 490                 197
## FBpp0289573                 282                 194                 196
## FBpp0074230                1686                2631                1789
## FBpp0074365                2162                2686                2394
## FBpp0311684                  71                 159                 123
## FBpp0074145                 177                 161                 265
## FBpp0074144                 167                 153                 257
## FBpp0078189                 136                 252                 129
## FBpp0309136                3312                4659                4532
## FBpp0312067                 329                 422                 409
## FBpp0311350               15921               17255               18862
## FBpp0071676                 483                1026                 881
## FBpp0077724                 289                 607                 475
## FBpp0082198                3096                3359                4439
## FBpp0304357                2145                 717                2114
## FBpp0076244                3390                2669                3661
## FBpp0290592                  38                  65                  28
## FBpp0302852               14333                1306               11116
## FBpp0084974                1041                 487                1283
## FBpp0307415                1045                1638                1103
## FBpp0082030                  21                  28                  24
## FBpp0300571                 656                 153                 568
## FBpp0310367                 201                   4                 748
## FBpp0289826                 231                 299                 318
## FBpp0082642                 538                 933                 646
## FBpp0305826                 144                  92                 117
## FBpp0086416                1123                1861                1274
## FBpp0099994                 117                 164                 156
## FBpp0303192                  84                 233                 108
## FBpp0078161                 319                 700                 304
## FBpp0304320                  28                  40                  59
## FBpp0083033                2066                2352                1986
## FBpp0073935                  50                  60                  26
## FBpp0304105                 133                 259                 122
## FBpp0309250                 110                 193                  83
## FBpp0310006                1803                1523                1752
## FBpp0071624                1195                2246                1486
## FBpp0301213                1460                1646                1906
## FBpp0084688                 624                 961                 691
## FBpp0291670                 418                 714                 479
## FBpp0289873                 226                 392                 278
## FBpp0306680                  35                   0                   1
## FBpp0080559                 263                 296                 259
## FBpp0080560                 369                 321                 340
## FBpp0082061                 124                 159                 118
## FBpp0300539                 319                1427                 417
## FBpp0099784                 372                1473                 457
## FBpp0070202                 352                 457                 486
## FBpp0073530                 139                 323                 183
## FBpp0075120                  69                 303                 140
## FBpp0305621                 141                 258                 170
## FBpp0113023                  40                  20                  47
## FBpp0088471                1087                1241                1268
## FBpp0304171                5456                4189                5496
## FBpp0089164               15536                3354               12140
## FBpp0309972                 613                 766                 697
## FBpp0297436                 109                  41                 140
## FBpp0075677                1043                1762                1260
## FBpp0305943                 792                 763                 881
## FBpp0077012                3551                1943                3913
## FBpp0305095               25957               30212               32878
## FBpp0072781                1286                1626                1350
## FBpp0081600                 446                 519                 574
## FBpp0081310                 443                 536                 634
## FBpp0311267                 709                1505                 858
## FBpp0073976                 495                 803                 678
## FBpp0112329                5960                5393                5801
## FBpp0074702              237599                  36              233554
## FBpp0310405                  15                   6                  13
## FBpp0308016                3293                9661                4092
## FBpp0088478                  51                 157                  38
## FBpp0304476                 928                2532                1070
## FBpp0292427                4554               14299                5659
## FBpp0301586                 370                 274                 368
## FBpp0071277                 754                1343                 838
## FBpp0078450                 462                 983                 679
## FBpp0308997                  77                 394                 105
## FBpp0083007                  64                 308                  85
## FBpp0078404                 588                 760                 665
## FBpp0089179                 955                1177                1212
## FBpp0305722                  95                 134                 137
## FBpp0297580                1842                1610                1722
## FBpp0086016                 244                 407                 279
## FBpp0309558                 259                 191                 281
## FBpp0089192                 673                 638                1006
## FBpp0085524                1121                1105                1496
## FBpp0308657                 259                 262                 283
## FBpp0308658                 231                 230                 257
## FBpp0085694                4441               10041                4879
## FBpp0309999                   0                   0                   3
## FBpp0087756                  96                 286                 215
## FBpp0072177               81605               45586               74888
## FBpp0085720              384466              246404              384398
## FBpp0309138                 432                 749                 586
## FBpp0309137                1095                1877                1398
## FBpp0087673                   0                  33                  82
## FBpp0312008                2744                2030                3243
## FBpp0077073                 128                 139                 173
## FBpp0077074                 133                 162                 179
## FBpp0305363                 606                 804                 763
## FBpp0111556                  39                   3                  58
## FBpp0293539                 434                 529                 522
## FBpp0077625                 321                 571                 473
## FBpp0074847                 139                 217                 178
## FBpp0075498                3074                5899                3562
## FBpp0311978               10953               10472               13235
## FBpp0099868                 237                 386                 238
## FBpp0077896                 526                 415                 558
## FBpp0290552                 674                1073                1452
## FBpp0304691                 888                1403                1932
## FBpp0292227                 822                1306                1810
## FBpp0309795                 261                 569                 366
## FBpp0073135                 239                 347                 174
## FBpp0110561                 948                 505                 692
## FBpp0082798                 460                 466                 572
## FBpp0076375                 919                 652                 667
## FBpp0303522                 205                 277                 185
## FBpp0079757                1393                1623                1640
## FBpp0084926                  74                  73                  87
## FBpp0084678                3563                 181                4370
## FBpp0291398                2686                5360                3420
## FBpp0302943                 180                   7                 264
## FBpp0291373                 247                  11                 375
## FBpp0099801                 433                  33                 428
## FBpp0302948                 493                  32                 710
## FBpp0301784                 401                  26                 419
## FBpp0075952                1081                1202                1170
## FBpp0075785                 808                2509                 743
## FBpp0075071                  85                 149                 133
## FBpp0087646                 488                 752                 646
## FBpp0081840                 412                 513                 487
## FBpp0072977                 493                 923                 710
## FBpp0309742                  83                  14                 100
## FBpp0310070               18664               18920               27340
## FBpp0293997                 534                 854                 587
## FBpp0086427                  37                  37                  59
## FBpp0111676                  64                  51                  84
## FBpp0071536                  53                 172                  32
## FBpp0290444                 600                 513                 653
## FBpp0074453                 205                 159                 238
## FBpp0081280                 581                1192                 631
## FBpp0307212                 843                1647                 919
## FBpp0076461                 281                 349                 403
## FBpp0073874                 197                 353                 324
## FBpp0310509                 375                 500                 437
## FBpp0307038                  62                1975                  85
## FBpp0301792                   1                   0                   1
## FBpp0086503                 954                1890                1749
## FBpp0086669                6164                 613                4624
## FBpp0300754                1537                2553                1608
## FBpp0288750                 595                 177                1076
## FBpp0309393                 219                  58                 403
## FBpp0309394                 618                 184                1139
## FBpp0292044                 210                 226                 186
## FBpp0311694                 776                1659                 841
## FBpp0300207                  70                 103                 113
## FBpp0306617                1171                8462                2092
## FBpp0081412                1197                1614                1232
## FBpp0311120                  88                  51                  40
## FBpp0089099                  87                  48                  40
## FBpp0306427                1124                1409                1496
## FBpp0297861                 713                 434                 451
## FBpp0312556                 228                  10                 202
## FBpp0072175                1908                  49                1754
## FBpp0072176                2117                  57                1948
## FBpp0309841                2833                2507                2714
## FBpp0081010                 928                2810                1162
## FBpp0305300                1056                1954                 915
## FBpp0081031                1577                2999                1443
## FBpp0297159                1158                2143                1002
## FBpp0311540                 310                 100                 188
## FBpp0074076                 791                 972                 837
## FBpp0077206                2352                2969                3027
## FBpp0309042                 368                 602                 423
## FBpp0302561               13422                 514               16752
## FBpp0077193                 176                 311                 204
## FBpp0084662                1783                2278                2258
## FBpp0071256                 572                 663                 570
## FBpp0306949                 511                 588                 513
## FBpp0087514                1997                3841                2212
## FBpp0081004                  29                  82                  31
## FBpp0305755                 412                 423                 428
## FBpp0075275                1297                1269                1234
## FBpp0082569                 817                 659                1123
## FBpp0079305                4654                4213                5878
## FBpp0305784                 647                1017                 796
## FBpp0075267                 464                 670                 490
## FBpp0112507                 869                1215                 614
## FBpp0073359                 170                 286                 145
## FBpp0306285                  64                  91                  84
## FBpp0309690                 925                1223                1264
## FBpp0076990                 986                 754                1425
## FBpp0310005                  10                  27                   9
## FBpp0077624                 151                 208                 220
## FBpp0288476                2283                3657                2512
## FBpp0086114                 891                1084                1048
## FBpp0297907                7306                9871                8631
## FBpp0305594                 892                1177                 960
## FBpp0311378                 135                 207                 185
## FBpp0070677                2279                1747                2125
## FBpp0301227                1107                1463                1327
## FBpp0304278                 134                 216                 229
## FBpp0304734                9874                6999                9592
## FBpp0304802                9382                 693                8815
## FBpp0311905                 330                 584                 426
## FBpp0111655                 594                1138                 749
## FBpp0308545                 704                1494                1064
## FBpp0310067                 307                 330                 183
## FBpp0076328              130512                 670              158271
## FBpp0293240                 485                 474                 488
## FBpp0070997                 169                 289                 246
## FBpp0099836                 534                 184                 910
## FBpp0076823                 205                   4                 471
## FBpp0308810                  16                  31                  24
## FBpp0086727                 268                 449                 353
## FBpp0071142                 573                 623                 697
## FBpp0293260                 858                1217                1285
## FBpp0308603                1650                1402                2065
## FBpp0080918                 109                 132                 120
## FBpp0087286                1378                   0                 252
## FBpp0072675                 267                 305                 335
## FBpp0310232                3490                6659                4884
## FBpp0078482                  14                   0                  24
## FBpp0082676                  43                   2                 109
## FBpp0078379                 551                 935                 630
## FBpp0110392                1161                 892                1067
## FBpp0311091                1272                 597                1012
## FBpp0084999                5915                2974                4754
## FBpp0084998               12130                6541               10705
## FBpp0304477               17591                9989               16115
## FBpp0071744                 978                1442                1013
## FBpp0302924                1557                2621                1846
## FBpp0292491                1687                2828                1951
## FBpp0083800                 130                 670                 129
## FBpp0077033                 256                 621                 319
## FBpp0303527                 161                 407                 233
## FBpp0083704                  21                  19                  27
## FBpp0306035                 783                 139                 731
## FBpp0074937                5275                6119                5729
## FBpp0311608                 570                 279                 784
## FBpp0080679                 637                1200                 772
## FBpp0088169                 567                 516                 517
## FBpp0304017                  91                 157                 102
## FBpp0100104                 698                1032                 724
## FBpp0301108                  84                  92                 110
## FBpp0071994                 112                 185                 116
## FBpp0084575                1216                1856                1243
## FBpp0310702                 270                 626                 403
## FBpp0071136                 479                 518                 669
## FBpp0290388                 829                 823                 766
## FBpp0088377                  65                  91                  78
## FBpp0303668                 235                 346                 262
## FBpp0084248                 270                 306                 153
## FBpp0088518                 285                 498                 380
## FBpp0073508                  22                  46                  35
## FBpp0290924                1535                2375                1608
## FBpp0111603                  95                 147                 110
## FBpp0111604                 317                 552                 404
## FBpp0289067                1051                1753                1536
## FBpp0291053                4153                6804                4855
## FBpp0291052                4426                7044                5195
## FBpp0306149                 190                 287                 285
## FBpp0099909                 240                 393                 368
## FBpp0293077                 455                 610                 585
## FBpp0312375                1438                2500                1787
## FBpp0076470                1324                1897                1621
## FBpp0306699                 687                1015                 803
## FBpp0306701                 283                 395                 327
## FBpp0310302                1372                1971                1200
## FBpp0080638                 744                1003                 705
## FBpp0301780                  50                  70                 125
## FBpp0080165                   3                   1                   1
## FBpp0074009                1790                 605                1877
## FBpp0305995               82374               10746               90041
## FBpp0070363                 888                 438                 802
## FBpp0070364                3635                2287                4145
## FBpp0070362                 910                 479                 841
## FBpp0084495                 611                 632                 873
## FBpp0084494                 533                 515                 715
## FBpp0087762                7352                4834                6935
## FBpp0292240                   3                  96                  14
## FBpp0113017                3435                  14                4657
## FBpp0075761                 755                 848                1044
## FBpp0074733                 303                 443                 327
## FBpp0291675                 257                 376                 274
## FBpp0311748                  18                   2                  57
## FBpp0084776                 109                 153                 102
## FBpp0077252                  28                  58                  44
## FBpp0075533                 165                 261                 146
## FBpp0088417                 327                 308                 407
## FBpp0297991                  42                 117                  41
## FBpp0291528                  72                 197                  92
## FBpp0074845                 472                2629                 630
## FBpp0309048                1082                1819                1312
## FBpp0077042                1598                1460                1484
## FBpp0308589                1711                1611                1580
## FBpp0078689                3689                4525                4634
## FBpp0289611                1207                1669                1707
## FBpp0291030                1460                 720                1324
## FBpp0072908                5800                7559                7526
## FBpp0074122                  38                  44                  44
## FBpp0099826                 125                  82                 110
## FBpp0070402                 342                 552                 377
## FBpp0310622                 825                 876                 710
## FBpp0111786                 157                 202                 171
## FBpp0303394                 313                 349                 305
## FBpp0080257                1788                2179                2395
## FBpp0078326                 624                1118                 740
## FBpp0302575                 672                 809                 703
## FBpp0110132                 518                 572                 490
## FBpp0070862                 157                 224                 181
## FBpp0077107                 484                 967                 598
## FBpp0290722                 345                 517                 415
## FBpp0291457                 937                1099                 913
## FBpp0077147                1231                2016                1062
## FBpp0292780                 215                 186                 306
## FBpp0292781                 336                 336                 519
## FBpp0289186                   1                   0                   0
## FBpp0296944               20298               10645               19428
## FBpp0296943                7535                9874                6752
## FBpp0292316                6244                4162                6133
## FBpp0306777                 171                 200                 194
## FBpp0079946                 772                 409                 925
## FBpp0080387                 312                 340                 371
## FBpp0085851                1185                 112                1319
## FBpp0086057                 172                 329                 286
## FBpp0306074                 766                 778                 982
## FBpp0306682                 970                 941                1199
## FBpp0307737                1737                 793                1927
## FBpp0081476                1199                 559                1386
## FBpp0291789                2428                1109                2931
## FBpp0303374                1225                 627                 979
## FBpp0079912                 109                  99                 181
## FBpp0271815               21870                2953               16094
## FBpp0087047                  56                  83                  58
## FBpp0311524                 305                 501                 325
## FBpp0308533                  62                  38                  45
## FBpp0294031                  61                  35                  43
## FBpp0087979                2706                4788                3604
## FBpp0082579                 211                 350                 268
## FBpp0081733                 319                 553                 340
## FBpp0086104                1398                2335                1217
## FBpp0304462                2576                6058                2694
## FBpp0297663                8925                8601               16270
## FBpp0082593                 134                 175                 196
## FBpp0078006                 144                 198                 215
## FBpp0073648                 131                 402                 134
## FBpp0306724                 322                 477                 408
## FBpp0304888                3157                4149                3170
## FBpp0304885                4226                5538                4343
## FBpp0084731                  69                 263                  99
## FBpp0084730                1514                6529                2353
## FBpp0087135                4355                4426                4449
## FBpp0111944                  48                  81                  64
## FBpp0303775                2387                2337                3318
## FBpp0311670                 474                 861                 449
## FBpp0309872                 407                 795                 472
## FBpp0290012                 455                1282                 704
## FBpp0300406                2147                1754                1930
## FBpp0311421                1024                1008                1543
## FBpp0073005                 229                 258                 375
## FBpp0073310                 750                 983                 821
## FBpp0291479                 545                 525                 526
## FBpp0074460                  66                  67                  92
## FBpp0291264                1973                2719                2072
## FBpp0087366                 306                 531                 313
## FBpp0076338                 333                 470                 359
## FBpp0294000                 804                 963                 731
## FBpp0080114                2590                2307                3134
## FBpp0301163                 165                 113                 133
## FBpp0084441                 585                 761                 605
## FBpp0071672                 807                1803                1128
## FBpp0300801                1223                1904                1564
## FBpp0309600                 600                 971                 622
## FBpp0073474                2297                3665                2293
## FBpp0074919                 230                 349                 246
## FBpp0074742                1343                 950                1273
## FBpp0074708                2877                1617                2889
## FBpp0305392                 373                2100                 573
## FBpp0081453                3408               20559                5572
## FBpp0077427                 117                 208                  81
## FBpp0084254                 445                 643                 631
## FBpp0309851                  39                  37                  93
## FBpp0305604                  44                  38                  67
## FBpp0311260                2519                2563                2963
## FBpp0305786                2318                 748                1339
## FBpp0305785               21284                6685               13516
## FBpp0081507                 423                 904                 692
## FBpp0083099                  19                  25                  24
## FBpp0303771                 460                 713                 564
## FBpp0303930                2358                1094                2079
## FBpp0303931                2784                1273                2453
## FBpp0291072                2744                1259                2402
## FBpp0111686                  32                  83                  24
## FBpp0304507                2175                3558                2641
## FBpp0077174                1077                1533                1345
## FBpp0303628               16140               22141               20296
## FBpp0076121                 441                 470                 494
## FBpp0086049                   7                  26                   8
## FBpp0087272                 423                 427                 510
## FBpp0306755                 638                1026                 695
## FBpp0307450                  34                  68                  35
## FBpp0304330                 441                 824                 477
## FBpp0084807                   0                   2                   0
## FBpp0080287                1404                1603                1615
## FBpp0073391                  33                  62                  54
## FBpp0085157                 282                 229                 234
## FBpp0290176                 660                 736                 761
## FBpp0099488                 297                 340                 334
## FBpp0081725                 480                 836                 705
## FBpp0077788                 231                 499                 255
## FBpp0309407                1019                2757                1688
## FBpp0289098                  14                  58                  15
## FBpp0086771                 293                 944                 409
## FBpp0077121                  48                 221                  43
## FBpp0075456                4448                2109                5327
## FBpp0075859                 684                 318                 828
## FBpp0304907                 347                 186                 337
## FBpp0084264                 601                1328                 670
## FBpp0312580                 744                 901                 576
## FBpp0075790                 830                1215                 956
## FBpp0075791                 806                1192                 945
## FBpp0304695                 189                 241                 207
## FBpp0087277                 473                 634                 519
## FBpp0307452                 617                1006                 845
## FBpp0087692                 141                 112                 138
## FBpp0304736                 104                3539                 134
## FBpp0087221                 417               14609                 634
## FBpp0085514                 209                 391                 266
## FBpp0089175                  37                  69                  47
## FBpp0073488                 238                 226                 254
## FBpp0312046                 580                 452                 567
## FBpp0305549                2287                1531                4313
## FBpp0082205               56886                 517               59005
## FBpp0304812                 794                 480                 949
## FBpp0077617                 308                 344                 444
## FBpp0082213                  90                  40                 160
## FBpp0304650                  49                  80                  77
## FBpp0086665                2567                 295                2850
## FBpp0077909                 367                 442                 395
## FBpp0304054                  38                  35                  45
## FBpp0304053                 200                 248                 238
## FBpp0304049                 177                 225                 217
## FBpp0288814                  32                  28                  38
## FBpp0072673                 432                 543                 498
## FBpp0084099                1353                3101                1476
## FBpp0083452                2009                1964                2309
## FBpp0288456                1878                3067                2052
## FBpp0311069                 894                1172                 998
## FBpp0082540                 202                 241                 164
## FBpp0304632                  28                   4                  10
## FBpp0304623                 136                  37                 131
## FBpp0072825                 131                 131                 170
## FBpp0112050                 633                 668                 805
## FBpp0305280                   2                   1                   2
## FBpp0086084                 257                 127                 174
## FBpp0084165                 557                 668                 573
## FBpp0293620                1849                1427                2229
## FBpp0292197                  56                 235                  71
## FBpp0290456                 458                 560                 439
## FBpp0310597                2461                 363                3760
## FBpp0308219                 203                 716                 284
## FBpp0073267                 363                 362                 394
## FBpp0311871                4126                4307                4767
## FBpp0305239                 154                 282                 236
## FBpp0309465                 554                3248                 891
## FBpp0311517                 331                 316                 326
## FBpp0086373                1906                5024                2780
## FBpp0304067                1651                1743                1553
## FBpp0305757                2569                5124                2595
## FBpp0303193                1233                2530                1222
## FBpp0087986                 268                1151                 353
## FBpp0082270                4169                5279                4836
## FBpp0077204                 722                 635                 715
## FBpp0290837                 937                1858                 942
## FBpp0087719                1097                1113                1077
## FBpp0087720                1205                1220                1163
## FBpp0087718                  73                  45                  64
## FBpp0293114                  28                  33                  39
## FBpp0087721                 890                 941                 877
## FBpp0300426                1280                1268                1231
## FBpp0307227                 173                 204                 212
## FBpp0078598                 994                1547                1200
## FBpp0308990                5148                6781                6048
## FBpp0074616                 185                 209                 168
## FBpp0075523                 719                 848                 674
## FBpp0309601                2937                 908                2868
## FBpp0086020                 643                2247                 751
## FBpp0309015                  76                 168                 160
## FBpp0070160                 281                 346                 334
## FBpp0077962                2061                3383                2543
## FBpp0073100                1040                1061                1132
## FBpp0099998                2297                2691                2498
## FBpp0074400                 566                 690                 552
## FBpp0076201                 222                 252                 301
## FBpp0087368                 648                 298                 633
## FBpp0100067                  98                  45                  88
## FBpp0086958                 423                 699                 574
## FBpp0073438                1934                2104                2170
## FBpp0312146                1031                1094                1125
## FBpp0070411                  60                 182                  59
## FBpp0307743                3548                2354                2783
## FBpp0072703                 774                1159                1002
## FBpp0311281                6970               14946               10279
## FBpp0082383                 312                 556                 349
## FBpp0077054               15704                 719               14784
## FBpp0310632                 171                 259                 225
## FBpp0075632                2115                1335                3009
## FBpp0290049                  90                  97                 100
## FBpp0083070                 687                 590                 802
## FBpp0083673                 437                 571                 545
## FBpp0312104               10055               10876               12341
## FBpp0307598                 119                 207                 148
## FBpp0080003                 884                1319                 950
## FBpp0087481                 235                 378                 292
## FBpp0301728                2314                2610                2686
## FBpp0082847                   5                   0                   6
## FBpp0307411                 269                 480                 619
## FBpp0072077                 928                1278                1113
## FBpp0297347               11328               15822               13651
## FBpp0089315                 389                 631                 489
## FBpp0305278                 174                 265                 210
## FBpp0306125                 203                 505                 250
## FBpp0085707                 465                1301                 583
## FBpp0310802                 506                 567                 507
## FBpp0079631                 837                 524                 670
## FBpp0298308                1224                7189                1665
## FBpp0087774                 468                 620                 618
## FBpp0080016                 523                 987                 731
## FBpp0075169                 157                 181                 165
## FBpp0304267                 489                 525                 479
## FBpp0297286                 408                 415                 372
## FBpp0304213                 139                 267                 228
## FBpp0075700                  91                  96                 128
## FBpp0309730                 741                1151                1019
## FBpp0080701                 272                 541                 395
## FBpp0306553                2958                  58                2873
## FBpp0088734                 120                   4                 131
## FBpp0304855                 146                  97                 188
## FBpp0304854                 303                 192                 370
## FBpp0310015                 562                 363                 654
## FBpp0072022                 437                 544                 440
## FBpp0081737                 169                 241                 153
## FBpp0072540                  33                  32                  28
## FBpp0303106                 202                 299                 136
## FBpp0308767                 199                 297                 136
## FBpp0309884                  48                  83                  56
## FBpp0305671                1698                2933                1611
## FBpp0290965                   3                   0                   0
## FBpp0076094                  17                  44                 111
## FBpp0312515                 382                 509                 475
## FBpp0080022                 451                 465                 548
## FBpp0071222                 225                 340                 281
## FBpp0292592                 397                 762                 427
## FBpp0080049                1485                3216                1283
## FBpp0082865                 549                 552                 613
## FBpp0070172                1702                2648                2649
## FBpp0304575                  66                 169                  48
## FBpp0297275                 115                 166                 111
## FBpp0077623                  88                 118                 118
## FBpp0088430                 305                 524                 433
## FBpp0086848                 598                 582                 680
## FBpp0071583                 284                 340                 367
## FBpp0291372               10272               13833               11113
## FBpp0070998                 163                 594                 220
## FBpp0294038                2782                3536                2653
## FBpp0073906                 332                 490                 461
## FBpp0305942                2692                2608                3215
## FBpp0073890                 345                 263                 373
## FBpp0293606                  98                 169                 137
## FBpp0083026                 114                 198                 162
## FBpp0071470                2548                1910                3156
## FBpp0303745                 258                 242                 293
## FBpp0081331                 173                 211                 176
## FBpp0081330                3646                5685                4620
## FBpp0289146                 974                 868                1157
## FBpp0303963                 218                 313                 202
## FBpp0079220                 673                 908                 682
## FBpp0311688                5308                1238                5168
## FBpp0305479                1946                1912                1950
## FBpp0071376                 408                 475                 550
## FBpp0070644                 147                 148                 107
## FBpp0085953                1134                 509                1253
## FBpp0113073                 390                 221                 462
## FBpp0070907                 292                 333                 322
## FBpp0083898                8844                4236                8303
## FBpp0083899                6950                2234                6528
## FBpp0072839               13480               30361               18561
## FBpp0311204                9235                3226                5801
## FBpp0290407                 956                 950                 970
## FBpp0077678                 163                 224                 195
## FBpp0306037                1149                1185                1044
## FBpp0304606                 471                1408                 714
## FBpp0305783                1091                 744                1215
## FBpp0310086                1827                1722                2077
## FBpp0300565                1605                1544                1844
## FBpp0306444                1634                2331                1865
## FBpp0074177                  26                   0                  12
## FBpp0311744                 720                1718                 899
## FBpp0113092                 314                 613                 350
## FBpp0305470                  62                  85                  87
## FBpp0074709                 451                 780                 580
## FBpp0309305                 442                 381                 355
## FBpp0083174                 399                 709                 521
## FBpp0304992                 179                  50                 115
## FBpp0078139                1189                1049                1416
## FBpp0303007                9051                1915                8250
## FBpp0305433                  95                 164                 158
## FBpp0271920                 101                 167                 160
## FBpp0085394                 295                1644                 255
## FBpp0289731                 324                1824                 277
## FBpp0290630                2161                2461                3298
## FBpp0073538                 670                1510                 759
## FBpp0078414                 526                 730                 638
## FBpp0309926                3406                 887                2404
## FBpp0089002                 176                 125                 200
## FBpp0100043                3412                6737                4414
## FBpp0082727                 737                 848                 766
## FBpp0070415                 246                 388                 275
## FBpp0080943                1232                1613                1649
## FBpp0085658                2236                3237                2622
## FBpp0075646                 667                1121                 814
## FBpp0077132                 523                 527                 566
## FBpp0082468                 347                 568                 357
## FBpp0311863                 256                 431                 311
## FBpp0302568                1038                1593                1301
## FBpp0071811                1013                1150                1208
## FBpp0082539                1050                4855                2072
## FBpp0111711                1853                5603                2528
## FBpp0297439                1433                2877                1886
## FBpp0297438                1875                4909                2571
## FBpp0307469                 351                 401                 366
## FBpp0088886                2235                2605                2408
## FBpp0079093                 223                 241                 270
## FBpp0073835                 674                 870                 843
## FBpp0306600                 182                 795                 162
## FBpp0112215                 210                 986                 193
## FBpp0310241                 625                 702                 790
## FBpp0290716                  10                   7                   5
## FBpp0310184                 641                 743                 801
## FBpp0111975                 209                 873                 323
## FBpp0087487                  83                 357                 129
## FBpp0070637                1060                1191                1082
## FBpp0301756                 574                 116                 608
## FBpp0289650                 299                  57                 270
## FBpp0309333                1582                 366                1497
## FBpp0311785                 386                 564                1050
## FBpp0304678                 191                 140                 254
## FBpp0304676                 221                 210                 199
## FBpp0110232                1793                1589                1905
## FBpp0078457                 795                 587                 875
## FBpp0078454                2131                1763                2217
## FBpp0111520                 169                 211                 185
## FBpp0111521                 194                 239                 232
## FBpp0309018                  62                 315                  76
## FBpp0077714                 673                 975                 541
## FBpp0290558                2178                2433                1806
## FBpp0084110                 824                 634                 931
## FBpp0085121                 720                 649                 832
## FBpp0071235                1026                1263                1028
## FBpp0311987               40830               42008               50085
## FBpp0112471                3016                3926                3003
## FBpp0302837                  81                  49                 100
## FBpp0083397                2053                1930                2459
## FBpp0110120                1046                1034                1315
## FBpp0090943                 883                1740                1050
## FBpp0311184                 254                 491                 480
## FBpp0079294                 372                 686                 500
## FBpp0312522                 166                 244                 177
## FBpp0082044                1971                2015                2922
## FBpp0305767                 398                 363                 493
## FBpp0289118                 384                 354                 478
## FBpp0305931                1858                1317                1626
## FBpp0311859                 546                1240                 789
## FBpp0304807                 113                 358                 149
## FBpp0293948                2154                2057                2480
## FBpp0309949                 718                 251                 935
## FBpp0112375                  51                 313                  65
## FBpp0310726                 115                 133                 153
## FBpp0303765                 358                 441                 421
## FBpp0297250                 291                 201                 274
## FBpp0078399               10746               16201               14571
## FBpp0308728                 730                 783                 823
## FBpp0308653                2028                1924                2129
## FBpp0303635                 664                 559                 656
## FBpp0079437                7400                7818               10340
## FBpp0085372                 216                 117                 279
## FBpp0306395                1829                2467                1842
## FBpp0083022                 308                 327                 359
## FBpp0293159                 497                 771                 680
## FBpp0293161                 499                 774                 677
## FBpp0311289                 229                 334                 287
## FBpp0293162                 233                 335                 290
## FBpp0071848                 904                1196                1083
## FBpp0070298                2006                2214                2576
## FBpp0083948                4565                6068                4085
## FBpp0306961                 901                1068                 807
## FBpp0077841                 328                 728                 245
## FBpp0301955                1266                1349                1612
## FBpp0307022                 360                 496                 363
## FBpp0296971                 418                 271                 462
## FBpp0075027                1674                3281                2169
## FBpp0308602                 122                 101                 101
## FBpp0078405                 552                 637                 693
## FBpp0077492                 467                 477                 573
## FBpp0076451                 767                1493                 979
## FBpp0078565                 213                 388                 257
## FBpp0080819                  25                   3                  44
## FBpp0303977                  40                  45                  71
## FBpp0310135                 637                 516                 623
## FBpp0086328                5414                9469                5725
## FBpp0073767                 133                 165                 160
## FBpp0070113                 258                 463                 344
## FBpp0087515                1190                1752                1361
## FBpp0082982                 247                 168                 164
## FBpp0310166                3794                5420                5014
## FBpp0290024                 167                   5                 305
## FBpp0300430                 188                   5                 364
## FBpp0300414                2025                  12                3908
## FBpp0070482                 290                 327                 388
## FBpp0080448                 500                 466                 604
## FBpp0075644                  89                  40                  78
## FBpp0289778                 120                  69                 134
## FBpp0293231                  48                 140                  34
## FBpp0080456                  18                   0                  25
## FBpp0077890                 193                 188                 178
## FBpp0304246                  96                 132                 116
## FBpp0074730                  78                 112                  80
## FBpp0075538                1009                 726                1182
## FBpp0305182               96674               92762              104120
## FBpp0084016                4956                9798                4671
## FBpp0071571                6309                3286                4057
## FBpp0071570                4874                2465                3102
## FBpp0304765                 618                 669                 779
## FBpp0110241                2160                2247                2534
## FBpp0091112                  16                  26                  21
## FBpp0084812                 803                1233                1118
## FBpp0312237                3224                2658                3345
## FBpp0311354                 923                 465                 864
## FBpp0309268                 136                  68                 125
## FBpp0071629                 407                1374                 563
## FBpp0079568                1975                3387                1975
## FBpp0307849                 806                1355                 825
## FBpp0297078                5295                4750                4623
## FBpp0081074                2919                2340                2210
## FBpp0083512                  18                  11                  26
## FBpp0310305                 291                 227                 331
## FBpp0310304                1032                 766                1110
## FBpp0075238                1530                2645                2090
## FBpp0308840                 158                 254                 184
## FBpp0078012                 173                 268                 205
## FBpp0304959                1474                2934                1856
## FBpp0075212                1231                2681                1698
## FBpp0080595                 385                 257                 356
## FBpp0297152                 691                1208                 606
## FBpp0071168                  66                 183                  82
## FBpp0310508                 284                 271                 433
## FBpp0082590                   8                 216                  18
## FBpp0088146                 741                 151                 436
## FBpp0291325               24337               16303               25526
## FBpp0087222                  12                   5                  16
## FBpp0087733                1099                1208                1154
## FBpp0293020                   7                   0                   6
## FBpp0099713                 675                 374                 895
## FBpp0111811                 175                 107                 286
## FBpp0305723                 154                 222                 183
## FBpp0086223                 253                 998                 280
## FBpp0072931                 404                 647                 476
## FBpp0070054                 427                 537                 602
## FBpp0309712                 598                1119                 646
## FBpp0307880                 544                  64                  76
## FBpp0307664                 382                  49                  69
## FBpp0307797                1565                2093                1475
## FBpp0311695                 237                 363                 215
## FBpp0081828                 963                1171                1006
## FBpp0077085                1565                2087                1473
## FBpp0082685                2527                7936                3255
## FBpp0082686                2439               11465                3461
## FBpp0298359                1379                 275                1279
## FBpp0302585                6397                4211                6258
## FBpp0311357               10948                 661               11631
## FBpp0088628                 445                 756                 430
## FBpp0088627                 998                1508                1052
## FBpp0081980                 432                 810                 590
## FBpp0074972                1547                2288                1939
## FBpp0305446                 242                 243                 295
## FBpp0083501                 782                1287                1060
## FBpp0305410                 376                 495                 294
## FBpp0111757                 470                 615                 369
## FBpp0305411                 467                 614                 369
## FBpp0086704                 345                 449                 325
## FBpp0078120                 523                 983                 470
## FBpp0305015                 294                 387                 258
## FBpp0305016                1107                1657                1196
## FBpp0289642                 758                1180                 866
## FBpp0296980                1917                2249                1907
## FBpp0072618                 174                 176                 174
## FBpp0072620                  36                  46                  34
## FBpp0306964                 279                 214                 237
## FBpp0074348                3717                3380                4012
## FBpp0085137                 261                 386                 329
## FBpp0072323                 812                 992                 924
## FBpp0078315                 224                 374                 342
## FBpp0301208                 222                 219                 238
## FBpp0090963                  46                 115                  41
## FBpp0308828                 253                 136                 232
## FBpp0077396                 268                 430                 294
## FBpp0304567                 383                 596                 460
## FBpp0089095                  80                 103                  91
## FBpp0310306                 512                1563                 729
## FBpp0310459                 418                 850                 421
## FBpp0300317                2955                5388                2426
## FBpp0084623               11582               21783               14643
## FBpp0075318                1177                1321                1559
## FBpp0075319                1071                1164                1389
## FBpp0291544                  29                  84                  55
## FBpp0304250                 143                 278                 210
## FBpp0084240                3053                3156                3346
## FBpp0307707                 202                 499                 245
## FBpp0311006                 129                 280                 141
## FBpp0083850                 690                 803                1001
## FBpp0293271                 704                 815                1016
## FBpp0071945                 966                1003                 882
## FBpp0086203                 387                 494                 363
## FBpp0311640               15225               16041               19228
## FBpp0310889                1566                2076                1411
## FBpp0308775                 225                 548                 312
## FBpp0308774                 341                 844                 473
## FBpp0082287                  26                  10                  23
## FBpp0289769                  82                 125                  79
## FBpp0078887                 260                 491                 295
## FBpp0303137                 454                 465                 397
## FBpp0292214                2366                5547                2819
## FBpp0076112                 778                 919                 811
## FBpp0303404                 809                 850                1200
## FBpp0083632                 839                 824                1185
## FBpp0301099                7669                9338                8676
## FBpp0301949                6937                8338                7766
## FBpp0306969                 808                1336                1007
## FBpp0303328                8037                1633                8247
## FBpp0298338                  43                  43                  47
## FBpp0086101                 325                 598                 381
## FBpp0304370                  84                 111                 108
## FBpp0303933                   6                  16                   5
## FBpp0289092                 251                 116                 189
## FBpp0079045                 341                 417                 499
## FBpp0086969                 410                 481                 507
## FBpp0088857                1909                1799                1682
## FBpp0088858                1998                1917                1733
## FBpp0291582                 881                1141                1119
## FBpp0112323                 202                 734                 228
## FBpp0305229                 406                 409                 376
## FBpp0307247                 156                 280                 208
## FBpp0308789                 927                1099                1093
## FBpp0086980                 429                 385                 577
## FBpp0085855                   2                   0                   4
## FBpp0087198                  11                   4                   0
## FBpp0301154                1326                1566                1402
## FBpp0309937                 722                1467                 721
## FBpp0087722                1142                1824                1397
## FBpp0311236                   3                   8                   8
## FBpp0076252                  16                  18                  17
## FBpp0297429                 481                 763                 415
## FBpp0078636                 844                1108                1240
## FBpp0290430                  75                  74                  90
## FBpp0290429                 175                 240                 209
## FBpp0077028               10189                 488               14044
## FBpp0111469                 923                1023                1106
## FBpp0074481                1019                 835                1109
## FBpp0085545                 771                 926                 778
## FBpp0086181                 611                1108                 782
## FBpp0086182                 596                1089                 775
## FBpp0071810                  31                   3                  20
## FBpp0072893                 324                1462                 519
## FBpp0309120                   2                  34                   7
## FBpp0072117                 878                1271                 995
## FBpp0079493                  33                  41                  21
## FBpp0087864                 900                1571                1024
## FBpp0307775                1370                2073                1806
## FBpp0089026                1132                1736                1500
## FBpp0089025                1474                2330                1930
## FBpp0311942                 805                1471                 979
## FBpp0074414                   2                   3                   3
## FBpp0072593                 980                1173                1013
## FBpp0310538                 350                 762                 657
## FBpp0078832                  19                 146                  21
## FBpp0308780                 900                1137                1050
## FBpp0084349                 920                1567                1213
## FBpp0079586                2338                4008                2491
## FBpp0309518                 357                 326                 385
## FBpp0075612               15959               15701               18664
## FBpp0079885                 431                 521                 501
## FBpp0079886                 911                1026                1008
## FBpp0301059                   3                   8                   4
## FBpp0311367                   6                   7                  38
## FBpp0088565                6313                6249                6761
## FBpp0076237                1292                3279                1800
## FBpp0081216                1179                 987                1283
## FBpp0303947                2522                3272                3388
## FBpp0293574                 268                 121                 228
## FBpp0083595                 750                 439                 624
## FBpp0311704                3360                2539                4157
## FBpp0081863                 331                 658                 431
## FBpp0088926               39716               79471               44090
## FBpp0088452               29105               56773               31597
## FBpp0084561                 127                 123                 159
## FBpp0289119                   2                   1                   1
## FBpp0305476                 267                 506                 283
## FBpp0271885                 236                 455                 254
## FBpp0293167                 440                 886                 520
## FBpp0310058               18358               17211               22339
## FBpp0070317                 340                 637                 375
## FBpp0290995                 479                 589                 601
## FBpp0271898                 990                2544                 903
## FBpp0087065                2272                 168                3185
## FBpp0304111                 121                 328                 287
## FBpp0304752                1126                1124                1463
## FBpp0304980                 995                 978                1264
## FBpp0072423                 214                1349                 249
## FBpp0304874                   0                   8                   2
## FBpp0306941                 619                 611                 695
## FBpp0089380                 254                 314                 268
## FBpp0088152                2315                2080                2187
## FBpp0088153                2337                2247                2360
## FBpp0310844                  34                  56                  40
## FBpp0078478                 460                 940                 511
## FBpp0079693                   1                   0                   1
## FBpp0310307                3462                7735                4520
## FBpp0297611                 298                 525                 442
## FBpp0076153               22850               21319               29006
## FBpp0309532               24603               23280               31584
## FBpp0290743                  18                 224                  34
## FBpp0077354                  63                  57                  56
## FBpp0310842                5106                2105                5345
## FBpp0084788                 143                 221                 190
## FBpp0113077                 459                 849                 439
## FBpp0309776                 118                 219                 107
## FBpp0087006                 366                 252                 262
## FBpp0086387                 355                 443                 352
## FBpp0079732                1035                1133                 981
## FBpp0087251                 450                 836                 611
## FBpp0085961                 952                1312                1313
## FBpp0082547                 517                 586                 649
## FBpp0086663                 713                1015                 324
## FBpp0297292                 314                 535                 374
## FBpp0311797                 705                1488                 948
## FBpp0073557                 564                1120                 726
## FBpp0305432                1662                1328                2149
## FBpp0305431                1419                1247                1658
## FBpp0072930                 131                 294                 182
## FBpp0073400                 628                2208                1422
## FBpp0289693                1033                 278                1090
## FBpp0080887                  15                  11                  13
## FBpp0086205                 170                 221                 260
## FBpp0083195                1479                1436                1585
## FBpp0080886                  23                  29                  28
## FBpp0308361                 654                 547                 855
## FBpp0086435                 848                1345                 810
## FBpp0099793                1657                1788                2512
## FBpp0082794                   0                   0                   0
## FBpp0307859                2900                3843                3043
## FBpp0303864                 835                1251                 809
## FBpp0077122                 453                 778                 524
## FBpp0112349                2473                2031                2129
## FBpp0307592                 701                 517                 857
## FBpp0112984                1399                1089                1688
## FBpp0311763                1717                3756                2250
## FBpp0290588                7073                6807                6172
## FBpp0303266               17341               16912               15067
## FBpp0291327                2524                2201                2293
## FBpp0084241                1045                2139                1164
## FBpp0080755                  94                  78                 172
## FBpp0308859                 181                 136                 371
## FBpp0073803                 120                   5                 209
## FBpp0300986                 180                  44                 152
## FBpp0290880                  64                  53                  80
## FBpp0088686                 834                1050                1163
## FBpp0088685                  39                  56                  38
## FBpp0289803                  38                  37                  33
## FBpp0304195                2721                3952                3111
## FBpp0110174                 785                1306                1039
## FBpp0076837               62779                2184               57699
## FBpp0076647                1427                3011                1582
## FBpp0300836                1919                1792                2396
## FBpp0100080                2633                3866                3004
## FBpp0079218                2829                 134                2871
## FBpp0078361                 115                 276                 125
## FBpp0077144                1336                1277                1693
## FBpp0082438                 350                 385                 328
## FBpp0085716                1751                1539                1681
## FBpp0074609                1378                1245                1538
## FBpp0071921                 872                  10                 652
## FBpp0082225                7199                  21                4984
## FBpp0088909                   2                   3                   0
## FBpp0310568                 275                 417                 347
## FBpp0303669                 574                 675                 575
## FBpp0086186                 501                 463                 557
## FBpp0303034                 911                 963                 866
## FBpp0087376                 403                 598                 544
## FBpp0087375                 412                 606                 550
## FBpp0298350                  15                   6                  14
## FBpp0075317                 856                1336                 975
## FBpp0072660                1195                1557                1622
## FBpp0297922                   6                  13                   4
## FBpp0077335                  24                  72                  33
## FBpp0302002                 122                 134                 137
## FBpp0271693                 363                 292                 337
## FBpp0312512                1129                 985                1051
## FBpp0303188                 250                 499                 460
## FBpp0074251                 696                1442                 832
## FBpp0089088               19246               10566               16910
## FBpp0077103                 218                 333                 273
## FBpp0310656                2029                 782                1786
## FBpp0086603                4827                5031                5422
## FBpp0080564                1045                1186                1349
## FBpp0303613                1073                 493                1275
## FBpp0090953                 976                1648                 941
## FBpp0306723                 953                1359                 887
## FBpp0090952                2592                3743                2691
## FBpp0083248                2768                6029                3637
## FBpp0072854                 116                 193                  98
## FBpp0303172                2583                1839                2706
## FBpp0288784               17112                9709               17671
## FBpp0306142                  53                  18                  47
## FBpp0311343                2585                 225                2183
## FBpp0309942                 270                 928                 382
## FBpp0297229                 245                 466                 335
## FBpp0078463                4087               10051                3817
## FBpp0293292                 144                 374                 265
## FBpp0304791                 195                 367                 490
## FBpp0306684                1255                 266                1582
## FBpp0293338                  98                1158                  71
## FBpp0305213                 136                1546                 126
## FBpp0304655                 117                 128                 118
## FBpp0309618                 564                 721                 543
## FBpp0305074                 682                 860                 651
## FBpp0309171                 911                1056                 863
## FBpp0290490                 778                1284                 895
## FBpp0070953                1257                1875                1558
## FBpp0077314                1922                1855                2796
## FBpp0079090                 136                 257                 322
## FBpp0290511                 125                 251                 143
## FBpp0085224                 340                 626                 346
## FBpp0077934                 868                1980                1017
## FBpp0306799                1207                2048                1379
## FBpp0078343                 268                 516                 366
## FBpp0082596               14021                 658               14572
## FBpp0087770                  24                  39                  43
## FBpp0291059                 620                2016                 770
## FBpp0289521                 787                1503                 919
## FBpp0300815                 256                 290                 274
## FBpp0111666                 594                 690                 705
## FBpp0304349                1452                 705                1662
## FBpp0311618                 120                  61                 134
## FBpp0289959                  16                   9                  11
## FBpp0291628                  53                  72                  81
## FBpp0293235                1404                2976                1120
## FBpp0304386                 217                 285                 277
## FBpp0290275                 112                 145                 112
## FBpp0089344                  54                  79                  50
## FBpp0309392                4180                2432                3911
## FBpp0110110                2350                  13                2198
## FBpp0289288                 215                 375                 237
## FBpp0311906                1502                1795                2215
## FBpp0112293                 109                 601                 126
## FBpp0112292                 411                2738                 543
## FBpp0310349                 160                1035                 174
## FBpp0076099                 198                 264                 167
## FBpp0307929                  72                  88                  81
## FBpp0071478                1896                1321                2340
## FBpp0292511                3629                4983                3791
## FBpp0303451                 662                1514                 647
## FBpp0073084                 148                  81                 165
## FBpp0073627                 242                 474                 287
## FBpp0077836                5444                5157                6627
## FBpp0304201                3729                3434                4453
## FBpp0074107                 530                 715                 633
## FBpp0310364                 113                 353                 105
## FBpp0080011                  23                  33                  56
## FBpp0076343                 899                 987                 768
## FBpp0305262                  87                 124                 105
## FBpp0075645                  69                 244                 149
## FBpp0290798                 233                 377                 299
## FBpp0293213                 618                1113                 660
## FBpp0078625                  37                  52                  40
## FBpp0312080               35615               35861               49096
## FBpp0303072                 217                 538                 250
## FBpp0306751                 218                 552                 261
## FBpp0085281                   2                   8                   6
## FBpp0293600                 164                 330                 316
## FBpp0099843                4649                6132                4762
## FBpp0310025                1541                1217                1678
## FBpp0305548                 106                 140                  98
## FBpp0071593                1207                1758                1500
## FBpp0111700                 933                 528                 781
## FBpp0306893                 591                 857                 597
## FBpp0071469                   5                  17                  23
## FBpp0081027                 518                 752                 604
## FBpp0086629                 747                2173                1195
## FBpp0303879                 302                 201                 285
## FBpp0303082                 250                 344                 304
## FBpp0072072                2683                2199                3024
## FBpp0079589                 173                 488                 521
## FBpp0307453                 155                 240                 205
## FBpp0077339                  14                  21                  19
## FBpp0293270                 255                 460                 287
## FBpp0304988                3339                7560                3063
## FBpp0301574                  17                  41                  12
## FBpp0081989                 407                 657                 449
## FBpp0304694                 475                 840                 762
## FBpp0305426                 326                 637                 442
## FBpp0307147                 182                 564                 181
## FBpp0082591                 424                 600                 570
## FBpp0080826                 501                 907                 619
## FBpp0310843                6099                 325               11999
## FBpp0303944                 774                 783                 949
## FBpp0083249                  36                  35                  22
## FBpp0303595                   9                   3                   4
## FBpp0294023                 529                 622                 559
## FBpp0294024                 578                 681                 630
## FBpp0077911                 536                 889                 756
## FBpp0080700                 481                 588                 519
## FBpp0306780                 519                 916                 704
## FBpp0312224               37397               42325               42091
## FBpp0099770                1846                2436                1887
## FBpp0071316                 144                 433                 185
## FBpp0304066                8380                7862                9900
## FBpp0078756                4406                4253                5080
## FBpp0312149                1756                1926                2204
## FBpp0083656               24869                 335               30931
## FBpp0271854                 166                 218                 248
## FBpp0084161                7443               10957                8021
## FBpp0301218                   3                  20                  14
## FBpp0083373                1221                1303                1467
## FBpp0083975                 237                 325                 244
## FBpp0078663                2232                3759                2119
## FBpp0082329                 123                  51                 213
## FBpp0086786                 935                 685                1107
## FBpp0070517                  11                   9                  20
## FBpp0291631                 161                  36                 146
## FBpp0085353                 600                 353                 596
## FBpp0085351                 615                 368                 587
## FBpp0079614                 218                 288                 269
## FBpp0078265                2175                6232                2784
## FBpp0075168                4582                1998                5269
## FBpp0081602                 174                 163                 321
## FBpp0074055                1186                1812                1293
## FBpp0075707               41998                  70               46100
## FBpp0303809                2732                4079                2858
## FBpp0086322                 495                 937                 591
## FBpp0074017               25887               35433               33716
## FBpp0307742                 816                 804                 815
## FBpp0307741                 215                 236                 220
## FBpp0310028               11051                6799                8884
## FBpp0075837                   6                  84                  26
## FBpp0289271                 165                 181                 187
## FBpp0071530                1338                2297                1412
## FBpp0305284                4089                3541                3678
## FBpp0087196                1090                 949                 981
## FBpp0087534                3554                1436                5098
## FBpp0081209                 308                 393                 395
## FBpp0073009                 113                 456                 158
## FBpp0308705                1112                 633                1056
## FBpp0083503                 298                 449                 368
## FBpp0080449               13762                 726               12648
## FBpp0078376                3790                3741                3453
## FBpp0072135                2726                3693                3057
## FBpp0082599                  12                  34                  12
## FBpp0073293                 161                 189                 288
## FBpp0291316                 328                 440                 374
## FBpp0080889                1027                 901                1119
## FBpp0305150                  59                 108                  68
## FBpp0309296                   6                  12                   9
## FBpp0309298                   4                   9                  14
## FBpp0297427                 902                1033                 747
## FBpp0076695                 299                 468                 357
## FBpp0304366                9064               14268                8914
## FBpp0302767                4108                3314                4533
## FBpp0073562                3463                2807                3872
## FBpp0303477                 139                 186                 169
## FBpp0303476                  93                 113                 124
## FBpp0088396                9317                9825                9139
## FBpp0110350                 106                  72                  75
## FBpp0110346                 190                 118                 143
## FBpp0297695                  74                  42                  65
## FBpp0110394                  71                  39                  59
## FBpp0088091                 197                 120                 148
## FBpp0288481                  55                 146                  98
## FBpp0305289                 873                1574                1130
## FBpp0071505                  15                  56                  35
## FBpp0071507                  25                  60                  57
## FBpp0290663                6185                6233                8054
## FBpp0305367                  10                  10                   3
## FBpp0306714                 921                1386                 970
## FBpp0306430                 544                 566                 708
## FBpp0071509                1688                3596                2312
## FBpp0302593                 228                 261                 223
## FBpp0077167                4473                2504                2821
## FBpp0113027                  68                 234                  72
## FBpp0304593                8493               23548               11146
## FBpp0306438                 335                 766                 451
## FBpp0308267                3121                5478                3546
## FBpp0304638                3709                5087                2697
## FBpp0304642                3751                4362                2670
## FBpp0309276                 534                 849                 631
## FBpp0306903                2605                5081                3914
## FBpp0294020                1365                1841                1549
## FBpp0087126                 208                 306                 235
## FBpp0075323                   4                   4                   5
## FBpp0301732                1260                2224                1664
## FBpp0292147                 405                 649                 764
## FBpp0307963                 143                 226                  59
## FBpp0303370                1153                1421                1377
## FBpp0111712                 333                 637                 637
## FBpp0070791                 139                 196                 182
## FBpp0086110                1415                2559                1659
## FBpp0085449                 125                  82                 147
## FBpp0311537                4012                4432                4586
## FBpp0074717                 230                  55                  55
## FBpp0310174                 456                 608                 460
## FBpp0298370               11649                  24               17394
## FBpp0310631                 286                 880                 395
## FBpp0087347                 317                 300                 300
## FBpp0072151                2337                4629                2497
## FBpp0080521                 456                1167                 567
## FBpp0312423                  10                  19                   6
## FBpp0079964                  16                  32                   9
## FBpp0073643                  12                  14                   8
## FBpp0306704                  19                  56                  22
## FBpp0309257                 492                2767                 750
## FBpp0099426                 291                 546                 290
## FBpp0303989                 175                 197                 172
## FBpp0303990                 231                 261                 231
## FBpp0309988                 248                 332                 377
## FBpp0304594                 807                1104                 945
## FBpp0079702                  88                 133                  85
## FBpp0079183                1121                 109                1339
## FBpp0086904                 903                1258                1245
## FBpp0311996                  28                  20                  20
## FBpp0080774                  11                   4                  25
## FBpp0072016                  37                  17                 104
## FBpp0079780                1052                1430                1360
## FBpp0083678                 197                  46                 214
## FBpp0073016                 713                 878                 900
## FBpp0297442                2300                2906                2831
## FBpp0305318                2889                4147                4202
## FBpp0303571                 418                 406                 406
## FBpp0113010                 132                 103                 197
## FBpp0307736                 524                 827                 468
## FBpp0307747                  89                 223                 119
## FBpp0311844                 823                 895                 897
## FBpp0083168                 213                 295                 277
## FBpp0312381                 507                 657                 624
## FBpp0306611                3629                5173                3678
## FBpp0310682                1358                1137                1747
## FBpp0076363                1520                1636                1839
## FBpp0311505                  35                  41                  61
## FBpp0311562                 894                 970                 960
## FBpp0304270                 761                 853                 729
## FBpp0305302                5740                5906                6634
## FBpp0305303                3215                3381                3463
## FBpp0080408                 256                 562                 373
## FBpp0080256                 286                 352                 271
## FBpp0297243                 581                1012                 829
## FBpp0082250                  23                  18                  23
## FBpp0306911                2339                2097                3417
## FBpp0312031                2021                2023                2434
## FBpp0303895                   6                   1                  13
## FBpp0100187               33050               43742               31862
## FBpp0087031                   3                  22                   8
## FBpp0309567                 839                1309                1055
## FBpp0076183                  11                 237                  40
## FBpp0298026                1404                  30                1395
## FBpp0082767                 377                 419                 417
## FBpp0071259                 624                 629                 578
## FBpp0306412                1392                1716                1584
## FBpp0086643                  75                  28                 228
## FBpp0306740                4979                4167                5227
## FBpp0073750                 210                 391                 331
## FBpp0308667                 349                 327                 274
## FBpp0303465                 651                1007                 836
## FBpp0077210                 784                1070                1117
## FBpp0306846                  88                 256                  25
## FBpp0310039                 942                1293                1104
## FBpp0081592                2049                2369                2613
## FBpp0075923                  14                  23                  16
## FBpp0082068                  63                   7                 108
## FBpp0082600                   6                  33                  21
## FBpp0079399                  56                 611                  93
## FBpp0309930                1955                2330                2228
## FBpp0078385                   1                   9                  17
## FBpp0084162                 314                 482                 418
## FBpp0303849                 563                 743                 694
## FBpp0310826                1008                1321                1386
## FBpp0088080                2022                1764                2405
## FBpp0310275                   4                 112                   3
## FBpp0074825                7125                4609                7902
## FBpp0071732                   3                  20                  23
## FBpp0086096                 398                 538                 563
## FBpp0083371               17499               16300               23084
## FBpp0070594                   0                   3                   0
## FBpp0078319                 340                 822                 398
## FBpp0312482                  47                   8                  67
## FBpp0087236                  19                  14                  68
## FBpp0080532                 459                 851                 550
## FBpp0306146                  11                  25                  16
## FBpp0071275                 615                 628                 782
## FBpp0076732                  46                   1                 115
## FBpp0071892                3619                3743                4084
## FBpp0082867                1627                1738                1837
## FBpp0081659                4209                3121                4615
## FBpp0304538                6742                9381                9389
## FBpp0307011                 699                 859                 912
## FBpp0082996                 566                 754                 712
## FBpp0089177                  52                  28                  64
## FBpp0071940                1666                1938                1848
## FBpp0303003                1252                 963                1047
## FBpp0075833                 320                 359                 335
## FBpp0300826                 308                 497                 605
## FBpp0303827                2717                2899                3253
## FBpp0306442                4251                5886                5213
## FBpp0308324               14126               14679               17455
## FBpp0085917               83075               77991               91413
## FBpp0305067                2077                2817                2735
## FBpp0306023                 489                 467                 616
## FBpp0305158                  84                 203                 108
## FBpp0307150                   9                   5                  10
## FBpp0110455                   2                  20                   2
## FBpp0305520                2064                2192                2889
## FBpp0076459                  96                  89                 109
## FBpp0085065                1317                1991                1749
## FBpp0082121                1356                2242                1597
## FBpp0073110                  16                  20                   3
## FBpp0088021                2911                2700                3536
## FBpp0086098               11382                9136               11748
## FBpp0084782                 133                 205                 184
## FBpp0084842                 562                1010                 696
## FBpp0082984                 588                 642                 849
## FBpp0087647                 254                 482                 297
## FBpp0099972                2432                2665                4141
## FBpp0080622                 493                 773                 560
## FBpp0087092                 213                 398                 217
## FBpp0075466                 160                  98                 204
## FBpp0290948                 294                 293                 380
## FBpp0087941                  82                 338                 204
## FBpp0070359                   0                   0                   0
## FBpp0078894                 181                 229                 181
## FBpp0303176                  78                 216                  10
## FBpp0305308                 686                1457                 801
## FBpp0071145                 120                 132                 119
## FBpp0084989                 423                 524                 614
## FBpp0070143               15660               15259               18776
## FBpp0073458                 838                1062                1109
## FBpp0085430                 491                 975                 643
## FBpp0084172                 303                 340                 265
## FBpp0112128                  95                 151                 163
## FBpp0293236                6802                7122                6902
## FBpp0305603                3959                5933                5222
## FBpp0291704                2099                3408                3412
## FBpp0072083                2089                1849                2400
## FBpp0307562                 582                 878                 733
## FBpp0298351                2051                2536                2207
## FBpp0070875                   1                   0                   6
## FBpp0072518                 248                 496                 330
## FBpp0084050                 308                 308                 443
## FBpp0080659                6067                5681                6668
## FBpp0307999                  27                  84                  30
## FBpp0304005                1146                2374                1396
## FBpp0078929                 498                 532                 539
## FBpp0079470                 260                 282                 341
## FBpp0074092                  16                  18                  32
## FBpp0310415                   0                   2                   2
## FBpp0305562               18717               17885               22071
## FBpp0305250                  21                 174                  19
## FBpp0081390                  25                  63                  21
## FBpp0076723                 360                 523                 422
## FBpp0089108                  55                  31                   9
## FBpp0085122                 129                 139                 159
## FBpp0099934               25848                  81               21385
## FBpp0071897                 593                 565                 773
## FBpp0072129                  93                 189                 113
## FBpp0304566                1779                2141                2299
## FBpp0309606                2782                3451                3657
## FBpp0079233                3061                4702                3620
## FBpp0086627                6949                5510                7887
## FBpp0301157                1476                1799                1919
## FBpp0307729                 293                 583                 383
## FBpp0100180               35770               43338               52859
## FBpp0071847                 710                 979                 892
## FBpp0111920                3508                4892                5051
## FBpp0070860                 339                 509                 401
## FBpp0087118                 526                 893                 465
## FBpp0079495                9105                7863               10964
## FBpp0087367                1930                2387                2243
## FBpp0304236                 511                1297                 598
## FBpp0077247                   5                 176                   1
## FBpp0084959               15930               16711               20933
## FBpp0087013                 443                 752                 638
## FBpp0307127                 966                1432                1293
## FBpp0289480                 557                 748                 701
## FBpp0309705                1470                1242                2122
## FBpp0083124                 336                 368                 400
## FBpp0304919                   7                   1                  13
## FBpp0080715                   9                   2                  13
## FBpp0310904                  54                 129                  85
## FBpp0290000                 634                1091                 723
## FBpp0112156                  35                  96                  57
## FBpp0289181                 311                 369                 421
## FBpp0084950                4933                5595                5562
## FBpp0070817                  43                  29                  49
## FBpp0086701               28951               26962               36971
## FBpp0072144                1089                 950                1323
## FBpp0310558                1265                 615                1367
## FBpp0074662                4398                5088                4679
## FBpp0311889                 166                 314                 214
## FBpp0071223                 186                 298                 213
## FBpp0305858               18111               17818               20855
## FBpp0079979                1458                1506                1530
## FBpp0311474                3131                4800                3706
## FBpp0087346                 392                 384                 576
## FBpp0310165                 618                 705                 601
## FBpp0306592                 115                 253                 127
## FBpp0309710                  12                   6                  15
## FBpp0306426                 784                1480                 728
## FBpp0087463                 239                 236                 301
## FBpp0303030                 823                1151                 972
## FBpp0087870                2311                4273                3230
## FBpp0305141                2507                3031                3079
## FBpp0291643                 444                 932                 494
## FBpp0078061                  14                   4                  18
## FBpp0310943                 422                 373                 552
## FBpp0307559                  32                   0                  33
## FBpp0085875                 251                 333                 339
## FBpp0074191                 299                 379                 325
## FBpp0310331                 300                 374                 405
## FBpp0083687                  44                  29                  35
## FBpp0076686                 309                 352                 369
## FBpp0307934                 239                 204                 279
## FBpp0099679                 577                 765                 723
## FBpp0311983                2309                2861                2849
## FBpp0297140                 910                1238                1165
## FBpp0070703                  99                 118                 261
## FBpp0074863                  57                  59                 131
## FBpp0308792                 344                 236                 378
## FBpp0070295                 222                 258                 225
## FBpp0083415                1022                1489                1203
## FBpp0077308                5050                4901                6407
## FBpp0075970                 118                 181                 161
## FBpp0084774                4894                2533                5896
## FBpp0079375                 883                1798                1036
## FBpp0311396                 545                 803                 585
## FBpp0311481                8888                8770               10745
## FBpp0310769                  77                  29                 129
## FBpp0304189                 940                2025                1278
## FBpp0072941                 513                 716                 618
## FBpp0085718                   1                   2                   0
## FBpp0083843                2395                2291                2906
## FBpp0079324                1415                2174                1835
## FBpp0298273                  64                  80                  58
## FBpp0073316                 541                 821                 611
## FBpp0081245                 346                 402                 409
## FBpp0312000                2449                2700                2105
## FBpp0302563                 308                 427                 334
## FBpp0305750                 417                 543                 459
## FBpp0077735                 740                 795                 742
## FBpp0076704                   1                  58                   0
## FBpp0085690                1468                2718                1695
## FBpp0113033                  38                  23                  29
## FBpp0072463                 200                 419                 261
## FBpp0082735                2160                4083                2706
## FBpp0079219                 449                 782                 510
## FBpp0300512               15104               14862               19180
## FBpp0087734                 413                 668                 461
## FBpp0070716                 516                 836                 622
## FBpp0311458               16559               17666               22318
## FBpp0306837               21472               21000               24760
## FBpp0081350                 122                 260                 129
## FBpp0077145                 214                 308                 187
## FBpp0308582                 256                 734                 345
## FBpp0083962                 186                 267                 244
## FBpp0086468                2284                2525                2865
## FBpp0089414                  79                 202                 142
## FBpp0072691                1145                1623                1411
## FBpp0075349                1313                2742                1771
## FBpp0071497                1137                1819                1465
## FBpp0078246              240440                9115              253082
## FBpp0290448                 692                1222                 922
## FBpp0080120                  54                  74                  72
## FBpp0081371                3170                6132                3523
## FBpp0309989                 781                1147                 932
## FBpp0081372                   2                   0                   1
## FBpp0292398                1359                3289                1661
## FBpp0306090                  24                   6                  28
## FBpp0309009                   5                  11                   9
## FBpp0292508                 178                 267                 254
## FBpp0078447                1596                2304                1845
## FBpp0083799                1500                 314                1685
## FBpp0305395                 269                  36                 500
## FBpp0070476                 172                 249                 206
## FBpp0082137                 526                 707                 621
## FBpp0071748                2090                1814                2351
## FBpp0099824               15122                8598               16284
## FBpp0305495                 240                 282                 301
## FBpp0073430                6199                 207                1045
## FBpp0311691                3300                3529                4600
## FBpp0312315                 166                 342                 259
## FBpp0071846               16854               19222               20789
## FBpp0290696                 652                 898                 780
## FBpp0084505                   0                   0                   3
## FBpp0310390                1564                2410                2295
## FBpp0311384                 642                1040                 862
## FBpp0080628                  68                 172                  89
## FBpp0070949                 428                 373                 478
## FBpp0074909                2699                 190                5084
## FBpp0085902                1095                1257                1284
## FBpp0077538                 275                 392                 368
## FBpp0291744                  20                   3                  21
## FBpp0311526                 240                 237                 173
## FBpp0083928                1768                2335                2271
## FBpp0305267                2168                2921                2983
## FBpp0071451                2106                2356                2499
## FBpp0291478                1056                1195                1104
## FBpp0308496                3592                5089                3943
## FBpp0080648                 940                 950                1350
## FBpp0302815                  94                 345                  88
## FBpp0087969                2238                2115                2800
## FBpp0079041                1088                1055                1235
## FBpp0309738                2460                3716                3191
## FBpp0312179                 140                 187                 196
## FBpp0300667                1356                1922                1795
## FBpp0301600                2599                3767                3487
## FBpp0309234                 968                2461                3101
## FBpp0292605                   5                   6                   9
## FBpp0309201               19864               21538               27245
## FBpp0309477                5020                4426                5419
## FBpp0310207                 722                1027                 929
## FBpp0306002                 573                 905                 684
## FBpp0077806                 173                 191                 188
## FBpp0083451                 999                1077                1263
## FBpp0311966                  11                   7                  22
## FBpp0086465                1254                2334                1443
## FBpp0306603                3399                2739                4774
## FBpp0303937                 161                 751                 301
## FBpp0076098                 387                 638                 441
## FBpp0072349                  14                  21                   3
## FBpp0304055                   1                   2                   7
## FBpp0305777                 306                 478                 348
## FBpp0087957                  58                  55                  90
## FBpp0311728                  36                   2                  19
## FBpp0085166               30316               32694               40989
## FBpp0084911                1186                3268                1535
## FBpp0081867                  49                  55                  61
## FBpp0305959                 366                 657                 351
## FBpp0311265                  43                  35                  39
## FBpp0085260                 703                1460                 951
## FBpp0075034                1709                2354                1640
## FBpp0113056                 345                 432                 337
## FBpp0079492                   6                   0                   0
## FBpp0305334                1891                2200                2456
## FBpp0081860                  52                  63                  66
## FBpp0310433                1294                2076                1661
## FBpp0084528                 383                 367                 572
## FBpp0075382                2373                2469                3296
## FBpp0312205                5151                2621                5706
## FBpp0306036                 503                1143                 588
## FBpp0087859                 842                1232                 972
## FBpp0308731                1323                1685                1359
## FBpp0076545                 169                 129                 200
## FBpp0304646                1595                1403                1634
## FBpp0304645                3455                3144                3779
## FBpp0305717                7896               14859               10951
## FBpp0086269               16705               17166               20804
## FBpp0304381                  12                  45                  31
## FBpp0312542                   6                  20                   6
## FBpp0307760                 369                1201                 474
## FBpp0289972                 619                 790                 688
## FBpp0086235                   0                   3                   1
## FBpp0305836                 798                 862                 941
## FBpp0311779                8191                7990                9053
## FBpp0081879                1235                1456                1750
## FBpp0072146                 231                 466                 336
## FBpp0074693                   0                   1                   0
## FBpp0087241                3926                3433                3858
## FBpp0302861                 113                 782                 202
## FBpp0077933                   9                   0                   0
## FBpp0306868                1084                1325                1593
## FBpp0301709                   3                   0                  10
## FBpp0311461                4466                4256                5439
## FBpp0112463                 221                 337                 283
## FBpp0085703                3885                4804                5022
## FBpp0291141                   1                   9                   0
## FBpp0112271                   0                   0                   1
## FBpp0078701                   0                   0                   0
## FBpp0080335                 467                1416                 765
## FBpp0087227                 758                1106                1002
## FBpp0071825                  56                 346                  58
## FBpp0297872                   0                  13                   0
## FBpp0309036                 341                 511                 497
## FBpp0086841                 301                 555                 383
## FBpp0111805                 382                 559                 644
## FBpp0083502                 619                 957                 790
## FBpp0075854                 897                1627                1347
## FBpp0110410                 642                1175                 759
## FBpp0072029                  51                  26                  62
## FBpp0077839                 336                 531                 368
## FBpp0071981                  52                   8                  35
## FBpp0086795                 519                 724                 640
## FBpp0073974                1287                1107                1718
## FBpp0311639                1670                  20                2211
## FBpp0070306                 300                 403                 321
## FBpp0292879                  64                  27                 109
## FBpp0078416               14623               14757               18013
## FBpp0312189                 516                 692                 857
## FBpp0074513                 994                1340                1262
## FBpp0306392                 431                 562                 528
## FBpp0083972                 480                 408                 521
## FBpp0292258                 243                 203                 250
## FBpp0079892                1451                3117                1877
## FBpp0113050                   2                   0                   2
## FBpp0073074                   3                   8                   3
## FBpp0300658                  17                  29                  32
## FBpp0086334                   7                 107                  11
## FBpp0085560                  20                  66                  52
## FBpp0099382                  11                  11                   6
## FBpp0290229                  63                  94                  87
## FBpp0072687               25607               25052               33467
## FBpp0073805                 185                 279                 220
## FBpp0076655                   1                   2                   8
## FBpp0308558                 839                 888                 957
## FBpp0083411                 406                 557                 478
## FBpp0071703                2773                1982                3257
## FBpp0071794               14064               21989               15939
## FBpp0082953                 250                 353                 260
## FBpp0076134                2780                4508                3589
## FBpp0300391                 990                1038                1347
## FBpp0082877                1416                1734                1814
## FBpp0073058                 180                 289                 219
## FBpp0079550                 960                1054                1247
## FBpp0306730                  30                 110                  30
## FBpp0309239                  76                 144                  99
## FBpp0086067                1998                2110                2086
## FBpp0075676                 356                 470                 440
## FBpp0113028                   0                   7                   4
## FBpp0290083                3409                4369                3662
## FBpp0310192                1032                1157                1294
## FBpp0081744                 249                 396                 250
## FBpp0076643                2236                3806                3131
## FBpp0091107                  46                  56                  61
## FBpp0301542                   0                   0                   0
## FBpp0086666                  72                  19                  81
## FBpp0306203                 142                  75                 263
## FBpp0311405                 258                 153                 613
## FBpp0292380                  97                 232                 121
## FBpp0308546                 102                 239                 122
## FBpp0311452               18776               18129               22564
## FBpp0088502                 123                 121                 145
## FBpp0305517                 121                 120                 141
## FBpp0071476                 806                 811                 969
## FBpp0076142                 686                1086                 929
## FBpp0304016                 236                 296                 313
## FBpp0110438                 576                 808                 696
## FBpp0079352                 227                 231                 413
## FBpp0073459                 336                 448                 328
## FBpp0310669                   8                   4                   3
## FBpp0308423                   8                  23                  30
## FBpp0309555                 345                 148                 380
## FBpp0304126                   4                   9                   8
## FBpp0291626                 258                 525                 407
## FBpp0310050                   0                   1                   0
## FBpp0071587                2069                2189                2748
## FBpp0085636                1427                 801                1613
## FBpp0110260                 502                 697                 888
## FBpp0100141                2176                1805                2604
## FBpp0303962                 840                1170                 870
## FBpp0083028                 130                 184                 130
## FBpp0072458                 144                  79                 223
## FBpp0082758                 524                 985                 746
## FBpp0312460                 350                 205                 309
## FBpp0297522               16197               15528               19254
## FBpp0070058                9768               20653                9873
## FBpp0086844                1669                2272                2095
## FBpp0304101                2940                2142                3328
## FBpp0081205                  43                   6                  42
## FBpp0081302                 343                 799                 399
## FBpp0311550                1966                1547                2038
## FBpp0077189                  39                 206                  44
## FBpp0300893                2251                1542                2602
## FBpp0082462                 528                 707                 727
## FBpp0088522                1323                1607                1766
## FBpp0075043                 694                1260                1024
## FBpp0070924                 253                 339                 261
## FBpp0112365                  58                 101                  98
## FBpp0077004                 422                1046                 532
## FBpp0070244                  86                 260                 126
## FBpp0312506                 271                 558                 380
## FBpp0312112                  13                  31                  12
## FBpp0311959                1046                1728                1293
## FBpp0070326                  28                  11                  29
## FBpp0087232                 304                 199                 363
## FBpp0305840                 670                1158                 779
## FBpp0081068                 773                1146                1357
## FBpp0290355                   0                   0                   1
## FBpp0082932                   0                   5                   3
## FBpp0311991                 640                 880                 829
## FBpp0311484                  29                  15                  63
## FBpp0290699                  10                   7                   3
## FBpp0303153                   2                   0                   0
## FBpp0073983                1227                1843                1512
## FBpp0070102                  71                 881                 226
## FBpp0075184                   8                   2                  18
## FBpp0082549                 365                 425                 484
## FBpp0075717                   8                   0                  11
## FBpp0084456                   2                   2                   4
## FBpp0084457                   1                   1                   0
## FBpp0082125                   2                   3                   2
## FBpp0289083                 411                 526                 467
## FBpp0308369                 724                 478                 816
## FBpp0077763                 600                 893                 679
## FBpp0289361                  61                  96                  65
## FBpp0082328                  19                  15                  14
## FBpp0312508                 233                 360                 349
## FBpp0309829                 606                 731                 903
## FBpp0077149                 512                 425                 665
## FBpp0070855                 114                  14                  70
## FBpp0083084                 101                  14                 115
## FBpp0099895                 443                 693                 912
## FBpp0085725                 828                1129                1007
## FBpp0301972                  16                  29                  39
## FBpp0082745                   3                   2                   1
## FBpp0110263                   5                   0                   3
## FBpp0084017                 365                 479                 442
## FBpp0081276                   0                   0                   0
## FBpp0072096                 546                 600                 757
## FBpp0305981                   0                   2                   0
## FBpp0087511                1212                1281                1394
## FBpp0304014                  98                   0                 244
## FBpp0076458                 135                 283                 211
## FBpp0311371                1301                2160                1723
## FBpp0298340                   4                   1                   1
## FBpp0307926                   1                   8                   0
## FBpp0082782                   6                   5                  10
## FBpp0305852                2956                4113                4219
## FBpp0079577                 257                 530                 309
## FBpp0080024                 323                 408                 452
## FBpp0072021                 574                 719                 636
## FBpp0099646                 248                 398                 157
## FBpp0310878                 266                 388                 435
## FBpp0297621                 262                 824                 596
## FBpp0304354                 350                 400                 362
## FBpp0290593                  12                  52                  34
## FBpp0302908                   0                   0                   0
## FBpp0304385                 244                 478                 380
## FBpp0297132               11241               23732               12195
## FBpp0082231                 121                 117                 108
## FBpp0309066                 637                1097                 804
## FBpp0289888                 378                 630                 577
## FBpp0081442                 105                 136                 126
## FBpp0081444                 830                1496                1345
## FBpp0311070                   9                   2                   7
## FBpp0075261                  92                 150                 104
## FBpp0306219                 933                1283                1151
## FBpp0309195                 655                 747                 622
## FBpp0306864                   0                   0                   0
## FBpp0302818                 390                 425                 416
## FBpp0288705                 417                 368                 510
## FBpp0304290                 457                 936                 560
## FBpp0081533                 272                 347                 309
## FBpp0312451                 629                 631                 685
## FBpp0088368                3186                2921                3419
## FBpp0303612                  14                   4                  30
## FBpp0088872                  20                   5                   8
## FBpp0312215                  22                   5                   8
## FBpp0311394                7254                6996                8917
## FBpp0084626                 509                 809                 443
## FBpp0071535                  93                 143                 126
## FBpp0290642                 668                 512                 733
## FBpp0080789                   4                   6                   5
## FBpp0293880                   0                   1                   5
## FBpp0086751                  98                  65                 183
## FBpp0075284                1514                 109                1649
## FBpp0310472                 423                 566                 413
## FBpp0074694                   0                   1                   0
## FBpp0086582                 226                1706                 335
## FBpp0083006                   1                   0                   0
## FBpp0308760                   3                   0                   6
## FBpp0083549                  13                  15                   4
## FBpp0303960                   0                   0                   1
## FBpp0305746                 211                  27                 243
## FBpp0073104                 388                  56                 420
## FBpp0078655                8893                7354               10968
## FBpp0303791                  48                  49                 122
## FBpp0309221                  58                  67                 140
## FBpp0099899                7386               11281                9606
## FBpp0086024                 410                 480                 430
## FBpp0292215                3314                3505                4007
## FBpp0080390                 620                 474                 853
## FBpp0304878                 487                 547                 650
## FBpp0074807                  24                   1                  15
## FBpp0304388                2837                3523                2402
## FBpp0071818                 333                 442                 487
## FBpp0309313                   1                   1                   5
## FBpp0309311                   0                   1                   1
## FBpp0112504                1061                1782                1249
## FBpp0301986                   5                  13                   1
## FBpp0309678                 932                1025                1050
## FBpp0311114                 129                 252                 151
## FBpp0312005                 353                 611                 346
## FBpp0311872                 119                 171                 221
## FBpp0289675                 199                 630                 301
## FBpp0300789                   0                   1                   1
## FBpp0078449                2315                2536                2791
## FBpp0070930                  46                  30                  37
## FBpp0077885                 341                 478                 387
## FBpp0077676                 159                 277                 195
## FBpp0079575                 241                 313                 325
## FBpp0074843                1177                1716                1747
## FBpp0305515                1606                2599                1918
## FBpp0311613               26140               41439               29978
## FBpp0084117                1303                  29                 830
## FBpp0081123                 113                 232                 156
## FBpp0311123                  12                  29                   8
## FBpp0110094                   1                   7                   1
## FBpp0082645                1052                1136                1192
## FBpp0312034                 572                 738                1199
## FBpp0082770                  31                 534                  40
## FBpp0072184                 560                 692                 637
## FBpp0077892                   4                   1                   3
## FBpp0085097                   4                  70                   6
## FBpp0087052                 134                 190                 168
## FBpp0288689                   0                   1                   0
## FBpp0288671                 438                 429                 453
## FBpp0078422                 414                 601                 561
## FBpp0079060                  26                 167                   9
## FBpp0293494                  14                  20                  29
## FBpp0309279                   0                   0                   1
## FBpp0306674                   0                   0                   0
## FBpp0304588                  24                 140                  44
## FBpp0312566                  19                 124                  39
## FBpp0308353                   0                   6                   0
## FBpp0078260                   7                   4                   2
## FBpp0087858                   7                   0                   5
## FBpp0306658                  19                  63                  28
## FBpp0099812                 955                1205                1063
## FBpp0072187                 313                 433                 402
## FBpp0297308                   1                   5                   0
## FBpp0297309                   0                   8                   0
## FBpp0297314                   0                   2                   0
## FBpp0083817                   0                   2                   0
## FBpp0307641                1635                2303                2232
## FBpp0088678                 635                1060                 939
## FBpp0088679                 488                 929                 779
## FBpp0099722                 237                1399                 381
## FBpp0071813                 566                 849                 628
## FBpp0088602                   7                  14                  19
## FBpp0312318                  15                  21                  22
## FBpp0307010                  22                  53                  51
## FBpp0083131                 275                 180                 310
## FBpp0304368                1051                 900                 896
## FBpp0303494                 891                2273                1462
## FBpp0292100                   3                   3                   3
## FBpp0305596                  13                   9                   8
## FBpp0290238                   7                  11                   4
## FBpp0309933                1024                1150                1085
## FBpp0084714                 328                  84                 540
## FBpp0088656                 864                2082                1221
## FBpp0088085                  27                  78                  57
## FBpp0080722                1216                 928                1270
## FBpp0075139                 281                 301                 353
## FBpp0303088                  47                  32                  20
## FBpp0303373                   0                   2                   0
## FBpp0306805                   2                   0                   2
## FBpp0076651                  25                   8                  40
## FBpp0099504                  80                 140                  83
## FBpp0292305                   0                  12                   0
## FBpp0305227                 153                 171                 222
## FBpp0305170                  16                  13                   4
## FBpp0078810                 542                 955                 655
## FBpp0077511                 912                1538                1118
## FBpp0075764               18335               19712               23855
## FBpp0305677                 931                 969                1459
## FBpp0306622                1060                 980                1288
## FBpp0076621                   3                   1                   4
## FBpp0075684                 252                 493                 319
## FBpp0083861                3374                3674                4544
## FBpp0087985                  89                  48                  89
## FBpp0290861                   2                   0                   8
## FBpp0074949                 529                 986                 658
## FBpp0070826                   3                   2                   3
## FBpp0075561                 697                1834                 859
## FBpp0311917                2870                2882                3203
## FBpp0084027                   6                  23                  11
## FBpp0086314                 529                 583                 675
## FBpp0303596                  30                  67                  41
## FBpp0071087                 101                  97                 101
## FBpp0290569                 806                1077                 990
## FBpp0076792                 431                 225                 486
## FBpp0306232               55959               48820               68491
## FBpp0082190                  16                 195                  37
## FBpp0311276                4708                6068                5376
## FBpp0087259                   8                   4                  13
## FBpp0303108                 424                1399                 617
## FBpp0306279                5633                4436                5228
## FBpp0292882                1851                1903                2244
## FBpp0304342                 219                 266                 210
## FBpp0293781                 209                 260                 210
## FBpp0073384                   5                  10                   5
## FBpp0303666                5136                5562                5506
## FBpp0072477                1207                1068                1481
## FBpp0070979                  16                  20                  16
## FBpp0308286                  15                 135                  43
## FBpp0297316                   4                   8                   2
## FBpp0308011                 277                 403                 363
## FBpp0311226                1357                2061                2034
## FBpp0085489                1460                2215                2198
## FBpp0306720                   1                   2                   1
## FBpp0304563                2626                2769                3171
## FBpp0075749                  19                  18                  48
## FBpp0311514                1660                 715                2230
## FBpp0079833                 423                 319                 521
## FBpp0070543                 555                 768                 755
## FBpp0303999                 229                 322                 312
## FBpp0112315                 922                 435                 959
## FBpp0112316                 717                 336                 743
## FBpp0086393                 331                 386                 377
## FBpp0309967                 309                 506                 396
## FBpp0303780               11065               11374               13548
## FBpp0087004                 585                 507                 895
## FBpp0303294                 730                 941                 787
## FBpp0082543                 270                 619                 310
## FBpp0080694                  98                 131                  95
## FBpp0292321                 351                 422                 430
## FBpp0310359                   2                   3                   2
## FBpp0310587                5172                   2                9870
## FBpp0073802                5333                   2               10230
## FBpp0076155                2183                1523                2689
## FBpp0082297                 674                 662                 902
## FBpp0291576                 425                 984                 580
## FBpp0088620                 276                 353                 366
## FBpp0090982                 285                 407                 309
## FBpp0292059                  56                  61                  32
## FBpp0304400                2124                1667                2541
## FBpp0075715                1216                1070                1448
## FBpp0304082                 560                 761                 712
## FBpp0070457                1323                 940                1427
## FBpp0301687                  66                  64                  95
## FBpp0309589                  38                  37                  51
## FBpp0075245                   1                   2                   0
## FBpp0308313                  37                  60                  37
## FBpp0305844                 499                 737                 649
## FBpp0086266               11820               19551                9531
## FBpp0080698                3740                   3                4616
## FBpp0086361                   0                   0                   0
## FBpp0075431                1224                1387                1464
## FBpp0290046                1393                2100                1532
## FBpp0080721                 207                 229                 278
## FBpp0310397                 650                 987                 836
## FBpp0112197                 120                 380                 137
## FBpp0077230                  54                 182                  60
## FBpp0304503                 460                 642                 532
## FBpp0074664                2301                1984                2567
## FBpp0311460               28048               32226               37510
## FBpp0112333                 590                 497                 607
## FBpp0078240                4499                1818                8811
## FBpp0084247                 264                 902                 232
## FBpp0310721                 220                 409                 347
## FBpp0072460                  44                 275                  74
## FBpp0078383                 969                1557                1134
## FBpp0305707                  37                  22                  20
## FBpp0081882                2285                3313                2614
## FBpp0074318                 303                 697                 401
## FBpp0086002                4227                7738                5245
## FBpp0081480                 234                 190                 235
## FBpp0081481                 261                 221                 273
## FBpp0304773                   2                   4                   0
## FBpp0088692                 678                 869                 666
## FBpp0081814                  22                 150                  22
## FBpp0271799                  13                  29                  36
## FBpp0289380                   0                   0                   2
## FBpp0303146                1772                1941                1716
## FBpp0310411                3422                4729                3803
## FBpp0302674                 173                 429                 224
## FBpp0079258                 134                 382                 164
## FBpp0290319                   5                   1                   6
## FBpp0290318                   1                   0                   0
## FBpp0072224                 295                 351                 310
## FBpp0077674                 108                 170                 140
## FBpp0081044                   0                   0                   0
## FBpp0077326                 250                 598                 254
## FBpp0083696                  92                  29                  69
## FBpp0081159                 941                1602                1145
## FBpp0305155                   0                   6                   2
## FBpp0309347                   4                   6                  17
## FBpp0303214                1296                1774                1712
## FBpp0311414                 285                   7                 327
## FBpp0306948                 296                   3                 296
## FBpp0087086                5574                6631                6603
## FBpp0073989                2015                2407                2400
## FBpp0080687                2847                4721                3384
## FBpp0110412               16518               16087               16181
## FBpp0074687                 265                 493                 361
## FBpp0073324                1130                   3                2080
## FBpp0304395                   1                   0                   0
## FBpp0081148                 306                 603                 303
## FBpp0293332                  11                   6                  15
## FBpp0293331                   0                   1                   3
## FBpp0073900                 134                 208                 159
## FBpp0304862                  25                  42                  22
## FBpp0081317                 544                1781                 690
## FBpp0072881                 142                  18                 179
## FBpp0110109                  37                   9                  45
## FBpp0304165                 890                1039                1045
## FBpp0112438                1584                1859                1834
## FBpp0305969                   2                   0                  13
## FBpp0071255                  30                  35                  56
## FBpp0311603                1516                2472                1768
## FBpp0110179                 251                 643                 229
## FBpp0290333                1755                 757                1807
## FBpp0070333                 480                 799                 626
## FBpp0311994                 518                 617                 705
## FBpp0297282                8744                2114               11100
## FBpp0081544                5079                2381                5640
## FBpp0085373                1520                1748                1948
## FBpp0301800                  67                 107                  68
## FBpp0074161                2230                2525                2626
## FBpp0311178                1497                1136                1663
## FBpp0088190                1635                1619                1707
## FBpp0305835                2403                2363                2585
## FBpp0082998                 204                 507                 249
## FBpp0083126                 461                 751                 538
## FBpp0070302                 191                 301                 221
## FBpp0083630                1223                1333                1694
## FBpp0310876                 657                 471                 803
## FBpp0074707                 419                 737                 502
## FBpp0076001                 366                 473                 456
## FBpp0085195                4462                6278                6779
## FBpp0292371                 202                 533                 404
## FBpp0112608                 317                 736                 492
## FBpp0081754                1033                1726                1090
## FBpp0073173                1338                2343                1516
## FBpp0307198                 635                 564                 772
## FBpp0087335                  76                 135                  99
## FBpp0083757                  34                   4                  32
## FBpp0304397                   0                   0                   0
## FBpp0087499                 506                 647                 625
## FBpp0078995                   3                   0                   1
## FBpp0070814                 397                 653                 530
## FBpp0306039               25352               21965               32808
## FBpp0311507                  46                 113                 467
## FBpp0110166               52892                 110               55253
## FBpp0075136               80158                 168               82713
## FBpp0081910                 114                 215                 137
## FBpp0075148                 503                 633                 533
## FBpp0080255                 266                 532                 466
## FBpp0071631                1240                2202                1858
## FBpp0085226                 235                  33                  19
## FBpp0070654                  19                  39                  23
## FBpp0307128                 112                 109                 224
## FBpp0077416                 176                 154                 301
## FBpp0308331                  14                   5                  19
## FBpp0079203                 101                 243                 159
## FBpp0088528                  25                  26                  20
## FBpp0289706                 492                1208                 582
## FBpp0302795                 493                1207                 594
## FBpp0086535                 684                1325                 842
## FBpp0306251                 150                 187                 183
## FBpp0082803                 368                 619                 425
## FBpp0303867                  64                 125                 101
## FBpp0099725                 840                1263                 973
## FBpp0099726                5321               10294                6610
## FBpp0075080                 376                 587                 684
## FBpp0070994                 594                1084                 908
## FBpp0074760                   9                   9                  12
## FBpp0288916                1461                 597                1591
## FBpp0288915                1187                 479                1245
## FBpp0087146                7050                 921                7025
## FBpp0305186                 328                 594                 794
## FBpp0077109                  68                   6                 189
## FBpp0308232                1475                2356                1945
## FBpp0084930                 438                 791                 593
## FBpp0288730                 402                  39                 265
## FBpp0080710                 334                  32                 239
## FBpp0082729                 369                 669                 504
## FBpp0082570                 409                 590                 468
## FBpp0079429                  87                 216                  61
## FBpp0302535                  69                  77                 102
## FBpp0302530                 664                1246                1190
## FBpp0072802               29866               29844               39586
## FBpp0304595                 204                 211                 615
## FBpp0090944                 215                 361                 272
## FBpp0290912                 818                1552                1123
## FBpp0311627                 211                 305                 516
## FBpp0082175                 230                  24                 537
## FBpp0082849                  94                 288                 119
## FBpp0309932                 366                 351                 459
## FBpp0079111                1079                1437                1049
## FBpp0082154                   6                  20                  15
## FBpp0088969                   9                   3                   3
## FBpp0309815                 972                1652                1140
## FBpp0086244                 543                1298                 748
## FBpp0086845                 228                 299                 308
## FBpp0076861                 295                 461                 281
## FBpp0084471               93771                2515               89728
## FBpp0084120                1582                1511                1984
## FBpp0087676                 495                 651                 510
## FBpp0078602                2740                3286                2852
## FBpp0082737                 627                 816                 754
## FBpp0306422                 568                1297                 864
## FBpp0311282                2020                2853                2193
## FBpp0302782                  52                  90                  68
## FBpp0305990                1878                2695                2059
## FBpp0077706                   0                   0                   2
## FBpp0083507                 612                 570                 659
## FBpp0072097               19265               20864               21764
## FBpp0312110                3211                5333                2905
## FBpp0079472                1652                2739                2280
## FBpp0082850                 163                 468                 292
## FBpp0071350                 101                 217                  88
## FBpp0111884                  10                  20                   2
## FBpp0081800                 525                 829                 618
## FBpp0293284                   4                   1                   1
## FBpp0312478                2648                3471                3415
## FBpp0305462                3815                2228                3879
## FBpp0293064                 424                 498                 810
## FBpp0083832                   0                   3                   1
## FBpp0311887                1708                1945                1963
## FBpp0074123                 341                 135                 660
## FBpp0081872                   7                  12                   6
## FBpp0293004                 171                 115                 206
## FBpp0083072                  85                  39                  94
## FBpp0076608                 413                1318                 466
## FBpp0081834                2622                2542                3010
## FBpp0310043                   3                  11                   3
## FBpp0076868                 254                 580                 320
## FBpp0291019                1655                3272                2040
## FBpp0071681                 337                 243                 242
## FBpp0311873                 111                 215                 128
## FBpp0081548                 429                 713                 485
## FBpp0072570                1146                2130                1600
## FBpp0072569                1213                2252                1781
## FBpp0113013                 313                 361                 389
## FBpp0309017                 183                 270                 248
## FBpp0289214                  71                  17                  84
## FBpp0312030                 130                1220                 167
## FBpp0309363                 177                 124                 125
## FBpp0075445                  31                 197                  12
## FBpp0307732                1871                 711                1977
## FBpp0311933                  91                  52                 124
## FBpp0304934                 838                1128                 698
## FBpp0071609                 430                 405                 546
## FBpp0077043                 411                 481                 560
## FBpp0081459                 448                 621                 564
## FBpp0309566                2560                1837                2940
## FBpp0308311                  52                 146                  78
## FBpp0087182                  80                 825                  93
## FBpp0309344                  28                   4                  27
## FBpp0080872                 955                1215                1142
## FBpp0308273                 115                 249                 174
## FBpp0306890                 237                 361                 365
## FBpp0302969                  22                  11                  36
## FBpp0070037                3132                4833                3715
## FBpp0308926                 835                1310                 964
## FBpp0309765                 847                1096                 749
## FBpp0079574                 188                 262                 139
## FBpp0304756                   5                   7                   2
## FBpp0311535                 105                  65                 372
## FBpp0307576                   7                   8                   8
## FBpp0091111                6415               11520                6796
## FBpp0071600                 560                 861                 757
## FBpp0311533                1403                1597                1833
## FBpp0111906                 539                 830                 498
## FBpp0311555                 473                 547                 520
## FBpp0071516                 673                1084                 966
## FBpp0309685                 675                1076                 949
## FBpp0110478                 354                 470                 418
## FBpp0088027                 304                 671                 475
## FBpp0079641                 243                1346                 347
## FBpp0086820                 227                 408                 324
## FBpp0075400                  37                  65                  64
## FBpp0311466                  18                   5                   3
## FBpp0312108                 311                  25                 164
## FBpp0084329                1360                 523                 921
## FBpp0085562                 369                 464                 430
## FBpp0075508                1132                3485                1785
## FBpp0307394                  20                   0                   0
## FBpp0303919                 825                 561                 835
## FBpp0087583                5596                8311                6045
## FBpp0311922                  67                 452                  46
## FBpp0300656                  94                 813                  64
## FBpp0083740                 773                1002                 876
## FBpp0072334                 109                 181                  94
## FBpp0297101                   0                   1                   1
## FBpp0307770               19049               19249               20332
## FBpp0088505               15158               15353               16238
## FBpp0302735                 217                 206                 275
## FBpp0307181                1660                1750                2098
## FBpp0303609                  10                  14                   5
## FBpp0305564                   9                  20                   7
## FBpp0080121                 275                 348                 309
## FBpp0070935                 158                 205                 145
## FBpp0307367                 247                 216                 159
## FBpp0082127                  14                  51                  37
## FBpp0085553                   1                   0                   2
## FBpp0303516                1063                1309                1091
## FBpp0303319                  15                  11                  23
## FBpp0100089                  23                  13                  34
## FBpp0077357                 102                 223                  92
## FBpp0077129                1573                1671                1971
## FBpp0305743                  69                 151                  65
## FBpp0073088                 428                 978                 459
## FBpp0085119                 445                 493                 580
## FBpp0073098                 765                 906                 866
## FBpp0087607                 815                1009                1032
## FBpp0081504                 344                 121                 502
## FBpp0073791                1061                1883                 812
## FBpp0099814                 116                 107                 223
## FBpp0080432                  23                  17                  51
## FBpp0073421                  19                 122                 111
## FBpp0305845                  26                  67                  64
## FBpp0304259                  20                  13                  13
## FBpp0081552                1147                2102                1472
## FBpp0311716                  20                  25                  13
## FBpp0302766                1450                1547                1845
## FBpp0087437                1291                1566                1993
## FBpp0083546                  56                 107                  76
## FBpp0306398                 197                  28                  71
## FBpp0077055                  33                   7                  96
## FBpp0293275                 973                1223                1118
## FBpp0311274                 177                 303                 282
## FBpp0310687                2615                3692                3382
## FBpp0080320                 959                1597                1104
## FBpp0078250                 366                   0                  33
## FBpp0293601                  29                  47                  33
## FBpp0309831                  33                  45                  54
## FBpp0075534                 104                 146                  98
## FBpp0309283                1475                 749                2091
## FBpp0307590                 255                 316                 120
## FBpp0074660                 680                 702                 702
## FBpp0074661                 256                 284                 284
## FBpp0074729                 131                 156                 142
## FBpp0288420                   1                   1                   0
## FBpp0070864                 381                 304                 436
## FBpp0082571                8060                4864                9800
## FBpp0084778                 271                 541                 290
## FBpp0087973                 353                 436                 347
## FBpp0087206                  49                 119                  83
## FBpp0075697                 918                1116                 838
## FBpp0305575                   6                   7                  11
## FBpp0083160                 482                 567                 558
## FBpp0083923                 335                 711                 454
## FBpp0307666                 934                3611                2293
## FBpp0079267                 242                 353                 271
## FBpp0112047                 349                 252                 344
## FBpp0309820                 261                 410                 320
## FBpp0297937                   7                   5                   2
## FBpp0309326                   0                   0                   1
## FBpp0297362                2236                2669                2783
## FBpp0303045                   2                  10                   7
## FBpp0073149                 300                 700                 420
## FBpp0306543                  25                  32                  28
## FBpp0084012                 244                 318                 308
## FBpp0088171                   1                   0                   0
## FBpp0073196                 320                 326                 367
## FBpp0077011                 426                 489                 478
## FBpp0311629                 149                 201                 140
## FBpp0081451                 129                 223                 111
## FBpp0081582                 631                 839                 862
## FBpp0077688                  40                   1                   8
## FBpp0086271                  13                   8                  13
## FBpp0073148                 471                 888                 807
## FBpp0086702                  63                  81                 106
## FBpp0309238                 755                 590                 569
## FBpp0292109                  42                  46                  71
## FBpp0075864                  13                  19                  14
## FBpp0075938                 741                 850                 836
## FBpp0081187                  19                  16                  25
## FBpp0303075                  61                 106                  78
## FBpp0311631                  11                  28                  19
## FBpp0309034                1208                1532                1589
## FBpp0292351                 192                 531                 259
## FBpp0084155                  52                 364                  58
## FBpp0113041                4296                5740                4880
## FBpp0310080                 323                 506                 417
## FBpp0085481                  10                  46                  45
## FBpp0304799                 122                  79                 106
## FBpp0082242                  16                  10                  25
## FBpp0072848                1687                2472                2023
## FBpp0291346                1051                1720                1228
## FBpp0087939                3478                2318                3675
## FBpp0088018                1218                2204                1431
## FBpp0083610                  15                  23                  31
## FBpp0074525                 808                1098                 921
## FBpp0305702                5676               14116                6845
## FBpp0297081                6201               15917                7459
## FBpp0073982                 323                 624                 381
## FBpp0303585                 403                 549                 563
## FBpp0087055                 102                 164                  78
## FBpp0290487                  82                 153                  91
## FBpp0072112                1067                1304                1302
## FBpp0288766                2380                2131                2939
## FBpp0309306                 733                 862                 931
## FBpp0293874                 172                 303                 185
## FBpp0073828                1289                1480                1331
## FBpp0271901                 606                 823                 591
## FBpp0309625                  92                 329                 119
## FBpp0289952                 660                1113                 689
## FBpp0083665                 381                 564                 444
## FBpp0080855                 137                 187                 118
## FBpp0079584                 111                 207                 169
## FBpp0078388                 348                 358                 437
## FBpp0306644                1817                2992                2130
## FBpp0297504                  91                  54                 103
## FBpp0072641                  72                  46                 140
## FBpp0083238                1270                1395                1364
## FBpp0073120                 181                 442                 252
## FBpp0311482                 616                 592                 849
## FBpp0082326                  32                  77                  32
## FBpp0309739                1455                7045                2369
## FBpp0078431                1282                1772                1531
## FBpp0084545                 356                 964                 461
## FBpp0311161                3034                1030                3410
## FBpp0079801                 173                 272                 212
## FBpp0312025                   0                   0                   2
## FBpp0305994                 196                 228                 147
## FBpp0301711                 146                 167                  99
## FBpp0307618                 989                1094                2057
## FBpp0304260              237053              254213              288473
## FBpp0082510                 597                 745                 754
## FBpp0311982                4783                4156                5015
## FBpp0293863                 189                 320                 175
## FBpp0303400                 429                 679                 624
## FBpp0079576                 594                1492                1044
## FBpp0305024                 619                 830                 640
## FBpp0288974                9308                 130               13616
## FBpp0309996                8580                 102               10802
## FBpp0312428                9031                 106               11322
## FBpp0084196                  20                  56                  76
## FBpp0079616                 731                 813                 758
## FBpp0311678                   7                  10                  17
## FBpp0086115                 457                 720                 547
## FBpp0072060                4135                1074                5027
## FBpp0309617                  54                  56                 129
## FBpp0083005                 159                 305                 157
## FBpp0308258                 260                 498                 354
## FBpp0301991                 370                 819                 462
## FBpp0081593                 307                 351                 393
## FBpp0086597                  89                 252                 112
## FBpp0309051                  59                 128                 130
## FBpp0289447                  26                  69                  62
## FBpp0307794                   5                   4                  10
## FBpp0296928                 493                 759                 602
## FBpp0082392                 415                 616                 480
## FBpp0291505                  11                   4                   9
## FBpp0086207                1439                2478                1914
## FBpp0076789                  39                 131                  64
## FBpp0078124                 198                 179                 205
## FBpp0310069                1204                1431                1413
## FBpp0290546                1164                1396                1368
## FBpp0310068                  16                   6                  10
## FBpp0081288                 132                  83                 120
## FBpp0308778                 264                 198                 216
## FBpp0073440                 215                 354                 159
## FBpp0078054                 409                 655                 451
## FBpp0072145                 592                 838                 946
## FBpp0289505                   8                   7                  14
## FBpp0082973                 540                1183                 650
## FBpp0289444                1376                2678                1685
## FBpp0293211                 573                 539                 639
## FBpp0099977                 333                 461                 419
## FBpp0291491                2770                2431                3347
## FBpp0089109                3807                2399                4412
## FBpp0311799                7287               14724               11802
## FBpp0312573                 967                 492                 989
## FBpp0088139                 125                 356                 124
## FBpp0304824                 879                1052                1029
## FBpp0289452                  48                  36                  71
## FBpp0086323                 192                 394                 242
## FBpp0310042                  10                   5                   2
## FBpp0111303                 428                 615                 520
## FBpp0087479               17533               20026               18886
## FBpp0072720                 478                 806                 742
## FBpp0073017                1191                1926                1504
## FBpp0291732                  55                 272                 148
## FBpp0288739                  41                  69                  49
## FBpp0087636                1040                 341                1679
## FBpp0089247                   1                   1                   6
## FBpp0307408                  56                 131                 160
## FBpp0297354                  49                 126                 145
## FBpp0083969                 189                 241                 164
## FBpp0311543                 613                 673                 624
## FBpp0086984                 143                 535                 159
## FBpp0081401                7676               10216                8257
## FBpp0083134                 410                 572                 588
## FBpp0083769                 157                 293                 196
## FBpp0302571                 539                 890                 726
## FBpp0076705                  11                  18                  20
## FBpp0076656                  48                  33                  62
## FBpp0088287                2959                1553                2982
## FBpp0310817                 375                  11                 254
## FBpp0303487                 123                 311                 212
## FBpp0074381                1894                1453                1874
## FBpp0309384                4211                4840                4007
## FBpp0304108                1024                1258                1280
## FBpp0288791                  56                 123                  89
## FBpp0288779                1945                3267                2469
## FBpp0291674                  42                  65                  52
## FBpp0309447                 307                 391                 285
## FBpp0087062                1256                 919                1296
## FBpp0099892                6715                1339                7717
## FBpp0306952                  18                   2                  58
## FBpp0080778                 120                 451                 201
## FBpp0304364                  24                  28                  18
## FBpp0111294                 438                 725                 874
## FBpp0075930                1167                2399                1457
## FBpp0080264                  10                   9                  39
## FBpp0081704                1013                1371                1168
## FBpp0309074                 734                 734                 926
## FBpp0086261                1251                1291                1642
## FBpp0088881                3949                5752                5959
## FBpp0300167                2245                 371                2980
## FBpp0300169                2309                 379                3066
## FBpp0290365                   1                   0                   1
## FBpp0070723                 384                 554                 333
## FBpp0083818                 907                 258                1035
## FBpp0082111                1336                1936                1925
## FBpp0079171                 706                3387                1089
## FBpp0308244                  18                  73                  21
## FBpp0311888                 484                 673                 561
## FBpp0304265               12378               12449               14552
## FBpp0305258                 256                1165                 331
## FBpp0298346                 377                 390                 517
## FBpp0077081                1434                1634                1503
## FBpp0077537                  57                  54                  44
## FBpp0086741                  72                  85                 207
## FBpp0290377                 347                  25                 113
## FBpp0290380                 279                  24                 101
## FBpp0089007                  50                   3                  15
## FBpp0075250                 277                 516                 359
## FBpp0085264                 677                 473                 707
## FBpp0309483                 589                 828                 660
## FBpp0074213                 949                1352                1091
## FBpp0085466                 107                  86                 179
## FBpp0099923                4575                2658                4464
## FBpp0085204                 375                1189                 450
## FBpp0070748                 251                 314                 295
## FBpp0081096                 478                 733                 541
## FBpp0074261                 724                 791                 927
## FBpp0305700                 523                 580                 531
## FBpp0311816                4733                5659                5131
## FBpp0307190                   0                   0                   0
## FBpp0085255                1162                1044                1211
## FBpp0112117                1509                4466                1639
## FBpp0086965                 218                 328                 295
## FBpp0291553                  32                 112                  22
## FBpp0309175                 533                 656                 584
## FBpp0083076                 749                 765                 935
## FBpp0072116                 298                 421                 371
## FBpp0075581                8887                3575                7986
## FBpp0307389                 138                 337                 190
## FBpp0088990                1318                3026                2245
## FBpp0304361                 233                 262                 255
## FBpp0073235                 300                 252                 383
## FBpp0084466                 405                 663                 528
## FBpp0077277                 897                1226                1009
## FBpp0312199                3386                2311                3481
## FBpp0086767                  97                 706                 119
## FBpp0303181                   1                   5                   2
## FBpp0298306                  47                  39                  24
## FBpp0088517                4569                2921                5046
## FBpp0289815               18204               16738               17776
## FBpp0311531                 845                1647                1040
## FBpp0076833                4176                5018                4542
## FBpp0311387                4353                5809                4272
## FBpp0293147                  16                  19                  19
## FBpp0293145                   9                   8                   5
## FBpp0293149                  97                 157                 122
## FBpp0304061                 193                 391                 226
## FBpp0075395                 400                 779                 370
## FBpp0075104                3098                3236                4379
## FBpp0309389                   2                   5                  12
## FBpp0309390                   3                  10                  25
## FBpp0086674                 455                 808                 501
## FBpp0072035                5878               10061                6785
## FBpp0082655                1228                1593                1426
## FBpp0309448                5719                4268                6707
## FBpp0087524                 450                 590                 457
## FBpp0304214                 749                1622                1034
## FBpp0070262               26373                2397               28221
## FBpp0297643                1145                 981                1695
## FBpp0291922                1110                 870                1611
## FBpp0291923                1228                 552                1779
## FBpp0291924                 709                 486                1055
## FBpp0305374                 375                 365                 390
## FBpp0305376                 182                 194                 205
## FBpp0086994                1263                3001                 544
## FBpp0087716                 745                2927                 756
## FBpp0308386                 696                2721                 693
## FBpp0271912                 849                 878                1182
## FBpp0074146                 173                  75                 112
## FBpp0110565                 126                  37                 317
## FBpp0293109                  65                 103                  69
## FBpp0304253                 952                1361                1184
## FBpp0306915                3418                4746                2915
## FBpp0112205                  32                 119                  48
## FBpp0073797                  19                  38                  12
## FBpp0075119                1748                1969                1997
## FBpp0079091                 701                 976                 904
## FBpp0100186               39952               44625               61613
## FBpp0072004                 268                 502                 281
## FBpp0078710                 290                 328                 255
## FBpp0311825               13193               11894               15054
## FBpp0073572                 729                1041                 786
## FBpp0304749                 678                1066                 844
## FBpp0298366                 460                 596                 822
## FBpp0076359               17154               18696               23136
## FBpp0075013                 173                 232                 185
## FBpp0075012                2226                3122                2606
## FBpp0293864                 650                 427                 758
## FBpp0310529                 142                 301                 206
## FBpp0311454                7536                8574                9630
## FBpp0110337                 831                 495                1282
## FBpp0290817                 727                 486                1144
## FBpp0085317                 233                 352                 290
## FBpp0308451                1173                2032                1426
## FBpp0303631                 459                 127                 618
## FBpp0079073                1682                1208                1634
## FBpp0071427                5474                9236                6588
## FBpp0288543                 188                 288                 221
## FBpp0304443                 436                 601                 515
## FBpp0075260                 477                 658                 579
## FBpp0071178                 546                 655                 604
## FBpp0087399                  19                  46                  45
## FBpp0290496                  49                 111                  81
## FBpp0305501                1665                1895                1736
## FBpp0082883                2922                  13                2635
## FBpp0308306                 440                 563                 458
## FBpp0072834                  10                  19                  33
## FBpp0072833                  34                  61                 115
## FBpp0312410                 143                 267                 150
## FBpp0310321                 926                1030                1031
## FBpp0070129                 210                 241                 196
## FBpp0099820                 289                 441                 305
## FBpp0085775                1096                1586                1426
## FBpp0309737                 240                 779                 432
## FBpp0307647                1593                2389                1801
## FBpp0088269                1227                1301                1388
## FBpp0312218                 229                 334                 274
## FBpp0311936                  65                1554                  89
## FBpp0303776                 344                 561                 360
## FBpp0077652                1054                1973                1342
## FBpp0070760                5892                3249                6134
## FBpp0304264                  62                 153                  93
## FBpp0308772                1398                2697                1683
## FBpp0291114                 269                  20                 317
## FBpp0291113                 274                  20                 322
## FBpp0291632               47043               43278               50657
## FBpp0075718                 463                 632                 592
## FBpp0081545                2395                1197                2808
## FBpp0071461                5099                7052                6251
## FBpp0071459                 788                1006                 926
## FBpp0292237                  80                 155                  39
## FBpp0073029                 929                 778                 475
## FBpp0073200                   7                   4                   3
## FBpp0303833                 333                 481                 424
##             M257_lg_male_hdhorn
## FBpp0087248                  23
## FBpp0293785                9060
## FBpp0080383                  96
## FBpp0077879                   6
## FBpp0311746                 176
## FBpp0289081                1133
## FBpp0311729                4062
## FBpp0085807                  43
## FBpp0081078                 349
## FBpp0312037                1246
## FBpp0302581                 558
## FBpp0084962                1889
## FBpp0311717                  35
## FBpp0301845                 156
## FBpp0070488                 432
## FBpp0070489                 431
## FBpp0070498                 571
## FBpp0079637                1524
## FBpp0307731                1271
## FBpp0072041                 646
## FBpp0085933                 842
## FBpp0310022                 101
## FBpp0310023                  38
## FBpp0073761                 142
## FBpp0075931                  21
## FBpp0293200                1227
## FBpp0076448                2556
## FBpp0076447                2266
## FBpp0290258                 764
## FBpp0087138                 242
## FBpp0300338                 615
## FBpp0307649                 556
## FBpp0289682                 167
## FBpp0306972                 413
## FBpp0306971                  27
## FBpp0290399                1723
## FBpp0301747                2787
## FBpp0070964                1472
## FBpp0307424                 978
## FBpp0071698                 210
## FBpp0070064                 927
## FBpp0309706                2447
## FBpp0307826                 532
## FBpp0289160                 366
## FBpp0301573                5653
## FBpp0307274                 719
## FBpp0306212                2772
## FBpp0078893                 267
## FBpp0297908                 927
## FBpp0084875                 621
## FBpp0072886               11324
## FBpp0304917                2245
## FBpp0084191                 530
## FBpp0075559                  59
## FBpp0300451                  29
## FBpp0300454                 228
## FBpp0300453                 283
## FBpp0301994                   0
## FBpp0081890                1106
## FBpp0309414                1040
## FBpp0099609                  28
## FBpp0085609                 203
## FBpp0271827                2121
## FBpp0292601                1635
## FBpp0304696                 411
## FBpp0087508                5194
## FBpp0087507                5110
## FBpp0306041                1854
## FBpp0088972                 755
## FBpp0301166                 269
## FBpp0304456                 204
## FBpp0072256                1772
## FBpp0305245                 447
## FBpp0309235               11893
## FBpp0072048                 838
## FBpp0309965                1500
## FBpp0072045                 348
## FBpp0290603                1673
## FBpp0079271                  87
## FBpp0292520                  19
## FBpp0110163                1023
## FBpp0076184                 873
## FBpp0306282                 289
## FBpp0311985                 284
## FBpp0289242                  65
## FBpp0085890                 143
## FBpp0085891                1276
## FBpp0290686                5436
## FBpp0290682                5315
## FBpp0311210                 870
## FBpp0110105                   7
## FBpp0305185                 462
## FBpp0086986                 239
## FBpp0079788                2098
## FBpp0077502                2008
## FBpp0082980                 386
## FBpp0081609                 386
## FBpp0312210                 267
## FBpp0082235                  32
## FBpp0072429                 162
## FBpp0072428                 729
## FBpp0305622                 797
## FBpp0072723                1771
## FBpp0080677                 370
## FBpp0309012                 484
## FBpp0086500                 533
## FBpp0085074                   4
## FBpp0085571                 562
## FBpp0078694                 612
## FBpp0311413                1780
## FBpp0083529                  22
## FBpp0310585                 357
## FBpp0311962                 310
## FBpp0291601                 769
## FBpp0111933                 119
## FBpp0306839                 649
## FBpp0077170                 672
## FBpp0306840                   4
## FBpp0304071                 438
## FBpp0300816                2514
## FBpp0289422                1084
## FBpp0099946                1192
## FBpp0308662                 866
## FBpp0112011                 652
## FBpp0085743                 978
## FBpp0073204                 392
## FBpp0290816                 169
## FBpp0305646                 223
## FBpp0290815                 403
## FBpp0111746                1949
## FBpp0271862                  19
## FBpp0086591                 583
## FBpp0073943                9448
## FBpp0309126                 171
## FBpp0082056                 439
## FBpp0089363                3891
## FBpp0084619                1817
## FBpp0310244                 722
## FBpp0073059                  61
## FBpp0293017                 336
## FBpp0072322                   4
## FBpp0072717                3846
## FBpp0300803                9562
## FBpp0070749                9538
## FBpp0297994                  62
## FBpp0099653                 343
## FBpp0305144                   0
## FBpp0086707                 157
## FBpp0071379                1521
## FBpp0076580                 408
## FBpp0073194                 263
## FBpp0291143               10909
## FBpp0070024                5843
## FBpp0305587                4222
## FBpp0072583                7453
## FBpp0304981                 669
## FBpp0083213                 989
## FBpp0088000                  77
## FBpp0087999                 199
## FBpp0291368                2855
## FBpp0078971                1679
## FBpp0308407                 192
## FBpp0084733                 557
## FBpp0305264                2050
## FBpp0300417                  46
## FBpp0297298                1160
## FBpp0081373                 279
## FBpp0082637                2133
## FBpp0079539                 274
## FBpp0078992                 698
## FBpp0311560                5269
## FBpp0079682                  15
## FBpp0087045                3310
## FBpp0310441                 364
## FBpp0083105                 782
## FBpp0070453                3395
## FBpp0089160                2766
## FBpp0309142                4147
## FBpp0303898                3421
## FBpp0290108                2430
## FBpp0306758                 663
## FBpp0083086                  42
## FBpp0304760                 642
## FBpp0074808                 917
## FBpp0309380                3798
## FBpp0074995                 104
## FBpp0089115                 206
## FBpp0086267                 198
## FBpp0289662                 450
## FBpp0081799                4128
## FBpp0312359                   0
## FBpp0306873                 542
## FBpp0306872                 555
## FBpp0292326                2321
## FBpp0307571                 409
## FBpp0088910                  78
## FBpp0311477                2692
## FBpp0083936                 165
## FBpp0081352                 284
## FBpp0072250                1858
## FBpp0297127                1819
## FBpp0085358                 112
## FBpp0288415                 100
## FBpp0306532                1647
## FBpp0099870                2047
## FBpp0290011                1441
## FBpp0300329                1017
## FBpp0293951                 462
## FBpp0309983                 906
## FBpp0075086                1781
## FBpp0303197                 332
## FBpp0292502                 506
## FBpp0074330                 400
## FBpp0077346                 866
## FBpp0309827                 428
## FBpp0079846                 256
## FBpp0081343                 723
## FBpp0071551                 129
## FBpp0077715                2286
## FBpp0082581                1605
## FBpp0304099                1941
## FBpp0073673                2177
## FBpp0301197                 303
## FBpp0304675                 466
## FBpp0304032                 127
## FBpp0308470                  63
## FBpp0306916                1102
## FBpp0306186                   9
## FBpp0306187                   3
## FBpp0310129                 677
## FBpp0073859                 746
## FBpp0087714                 482
## FBpp0074784                 161
## FBpp0082972                1462
## FBpp0071129                2661
## FBpp0293470                1088
## FBpp0087686                2710
## FBpp0304681                 146
## FBpp0074247                 491
## FBpp0309143                1116
## FBpp0073608                1410
## FBpp0307465                1356
## FBpp0309979                2174
## FBpp0289514                  12
## FBpp0081661                 133
## FBpp0304298                 730
## FBpp0080997                1693
## FBpp0303820                 211
## FBpp0084021                2151
## FBpp0291142                  40
## FBpp0292665                 601
## FBpp0292663                1868
## FBpp0075999                 715
## FBpp0112193                3442
## FBpp0309360                 168
## FBpp0305497                3013
## FBpp0292329                  63
## FBpp0085812                  87
## FBpp0081156                4037
## FBpp0086395                2788
## FBpp0305484                 347
## FBpp0290814                5116
## FBpp0111781                2536
## FBpp0084726                1390
## FBpp0087866                 226
## FBpp0290138                1802
## FBpp0290139                1637
## FBpp0307216                 143
## FBpp0079619                  86
## FBpp0307931                  86
## FBpp0071548                 901
## FBpp0303265                 100
## FBpp0304880                  75
## FBpp0087429                1246
## FBpp0087431                1163
## FBpp0289106                1005
## FBpp0087436                1117
## FBpp0074792                 608
## FBpp0308816                 557
## FBpp0297771                2177
## FBpp0308626                   7
## FBpp0293582                  55
## FBpp0083533                 485
## FBpp0312441                3541
## FBpp0308362                 490
## FBpp0081263                 877
## FBpp0073740                 588
## FBpp0308981                 239
## FBpp0309092                2163
## FBpp0074614                 139
## FBpp0071677                 297
## FBpp0292148                 705
## FBpp0302807                5424
## FBpp0074314                  27
## FBpp0080450                 141
## FBpp0305020                  73
## FBpp0297111                 123
## FBpp0076185                 194
## FBpp0085075                3139
## FBpp0304800                  58
## FBpp0311990                  74
## FBpp0301042                  70
## FBpp0085082                 248
## FBpp0081996                2217
## FBpp0304377                2746
## FBpp0081380                 383
## FBpp0293890                 277
## FBpp0075918               20288
## FBpp0300990                2077
## FBpp0300989                 883
## FBpp0086347                  18
## FBpp0304174                2230
## FBpp0088962                3235
## FBpp0072100                 187
## FBpp0309809                1216
## FBpp0110195                  44
## FBpp0100035                1634
## FBpp0305194                 571
## FBpp0081850                  37
## FBpp0070795                1489
## FBpp0086647                 549
## FBpp0293738                   0
## FBpp0303580                  35
## FBpp0086868                  36
## FBpp0312526                1761
## FBpp0291652               15639
## FBpp0074464                 114
## FBpp0071686                1145
## FBpp0307759               44009
## FBpp0291742                2702
## FBpp0081823               31621
## FBpp0086483                 126
## FBpp0087352                  90
## FBpp0073556                   4
## FBpp0086506                 411
## FBpp0086505                 386
## FBpp0310686                3607
## FBpp0085586               14796
## FBpp0308422                1955
## FBpp0072038                 119
## FBpp0071713                 193
## FBpp0292778                 183
## FBpp0309285                 125
## FBpp0309324                 554
## FBpp0080672               20934
## FBpp0085448                 138
## FBpp0087704                1894
## FBpp0111946                 153
## FBpp0078881                 213
## FBpp0083584                1230
## FBpp0070025                1220
## FBpp0086487                9988
## FBpp0085484               15082
## FBpp0309089                1063
## FBpp0312096                1448
## FBpp0306962                 617
## FBpp0079737                 250
## FBpp0306292                 977
## FBpp0085804                 251
## FBpp0088841                 223
## FBpp0082876                2398
## FBpp0083537                 858
## FBpp0292800                 905
## FBpp0305226                 743
## FBpp0089133                5670
## FBpp0290667                 714
## FBpp0297610                 493
## FBpp0304215                 330
## FBpp0078566                1690
## FBpp0073437                 682
## FBpp0084340                1337
## FBpp0309079                 683
## FBpp0113120                1457
## FBpp0293289                2956
## FBpp0290158                 485
## FBpp0310810                 308
## FBpp0073297                3214
## FBpp0085829                 645
## FBpp0088130                1712
## FBpp0080867                 225
## FBpp0292901                   0
## FBpp0083979                3648
## FBpp0307180                 185
## FBpp0081526                1011
## FBpp0081990                 601
## FBpp0087180                 252
## FBpp0301778                 197
## FBpp0070410                1467
## FBpp0289206                  41
## FBpp0072494               10396
## FBpp0303821                3324
## FBpp0073856                1141
## FBpp0079915                7068
## FBpp0305776                 831
## FBpp0300799                 187
## FBpp0080484                2195
## FBpp0077367                2118
## FBpp0304732                2493
## FBpp0075316                1744
## FBpp0289785                8747
## FBpp0289783                3856
## FBpp0078086                 403
## FBpp0303364                 163
## FBpp0303361                1500
## FBpp0303363                1484
## FBpp0075631                 495
## FBpp0075425                 627
## FBpp0087320                2432
## FBpp0076872                1225
## FBpp0305817                2373
## FBpp0305335                2682
## FBpp0085453                 172
## FBpp0077998                1769
## FBpp0292047                  72
## FBpp0074691                1257
## FBpp0306051                 114
## FBpp0087958                 671
## FBpp0303416                8804
## FBpp0300431                9835
## FBpp0300432                9838
## FBpp0089348                   1
## FBpp0309386                1055
## FBpp0290266                 964
## FBpp0073664                  30
## FBpp0297358                4545
## FBpp0084937                  54
## FBpp0072192                1007
## FBpp0271836                1558
## FBpp0087404                  43
## FBpp0293134                 102
## FBpp0079453                1887
## FBpp0072172                 402
## FBpp0302579                 185
## FBpp0070083                   7
## FBpp0308713                1239
## FBpp0073678                 319
## FBpp0304879               43925
## FBpp0305851                2019
## FBpp0290364                 830
## FBpp0305633                1954
## FBpp0081907                2559
## FBpp0081649                  77
## FBpp0310727                 911
## FBpp0082180                 449
## FBpp0290908                   0
## FBpp0290906                 125
## FBpp0086127                   7
## FBpp0305244                 238
## FBpp0070465                 143
## FBpp0311868                1698
## FBpp0071296                1344
## FBpp0305527                1644
## FBpp0289954                1339
## FBpp0084258                1481
## FBpp0308637                 948
## FBpp0077975                 199
## FBpp0084738                1686
## FBpp0079326                7928
## FBpp0086359                  89
## FBpp0293627                  84
## FBpp0308762                 438
## FBpp0289849                 201
## FBpp0083767                1293
## FBpp0077740                2375
## FBpp0290166                 391
## FBpp0298309                2466
## FBpp0082742                  45
## FBpp0080903               15855
## FBpp0311147                4709
## FBpp0307273                 968
## FBpp0311874                 566
## FBpp0290830                  34
## FBpp0309990                 349
## FBpp0304417               22358
## FBpp0111524                3734
## FBpp0075727                 324
## FBpp0305584                 561
## FBpp0075560                 351
## FBpp0292507                1094
## FBpp0084792                1675
## FBpp0084791                 488
## FBpp0310240                8015
## FBpp0311780                 223
## FBpp0081226                1545
## FBpp0291807                 839
## FBpp0297588                 844
## FBpp0291806                 568
## FBpp0080969                 525
## FBpp0087085                 709
## FBpp0303022                  52
## FBpp0302026                 172
## FBpp0289455                 172
## FBpp0079815                2219
## FBpp0083736                 125
## FBpp0076871                 222
## FBpp0081961                 111
## FBpp0086475                 750
## FBpp0077872                 923
## FBpp0082186                2466
## FBpp0075617                3107
## FBpp0292595                 918
## FBpp0297553                2302
## FBpp0074793                  73
## FBpp0080751                2348
## FBpp0080753                2242
## FBpp0308295                 614
## FBpp0307682                  62
## FBpp0074442                 565
## FBpp0310438                 781
## FBpp0075214                  11
## FBpp0312219                 448
## FBpp0311549                1522
## FBpp0085540                 742
## FBpp0089141               22300
## FBpp0110281               21142
## FBpp0303067                  18
## FBpp0303069                  18
## FBpp0084874                 637
## FBpp0070612                1240
## FBpp0072127                6207
## FBpp0311256                 846
## FBpp0111705                1784
## FBpp0085693                 466
## FBpp0076791                 523
## FBpp0289278                  71
## FBpp0292037                 102
## FBpp0271718                1842
## FBpp0312201                 259
## FBpp0074005                2288
## FBpp0087519                 751
## FBpp0080275                 220
## FBpp0307803                2828
## FBpp0293626                3134
## FBpp0308661                2970
## FBpp0079623                 456
## FBpp0082720                1537
## FBpp0297605                 121
## FBpp0086429                2389
## FBpp0081312                1113
## FBpp0311512                 425
## FBpp0307870                 407
## FBpp0088031                 498
## FBpp0301215                1383
## FBpp0087323                1730
## FBpp0072564                2680
## FBpp0302692                2660
## FBpp0086421                1358
## FBpp0309704                1057
## FBpp0310002                 335
## FBpp0302006                  19
## FBpp0304310                  60
## FBpp0076921                 158
## FBpp0290353                 225
## FBpp0078673                 484
## FBpp0076599                 483
## FBpp0083846                2638
## FBpp0305307                1245
## FBpp0110549                 488
## FBpp0077448                1501
## FBpp0306201                 712
## FBpp0306862                 398
## FBpp0080148                 121
## FBpp0289644                 615
## FBpp0292236                1057
## FBpp0310630                1119
## FBpp0292244                1126
## FBpp0290516                 732
## FBpp0289110                  11
## FBpp0301793                  29
## FBpp0293619                2152
## FBpp0293617                1862
## FBpp0290577                2679
## FBpp0082135                 170
## FBpp0071969                 306
## FBpp0071028                1256
## FBpp0079304                 484
## FBpp0291536                3857
## FBpp0289568                 831
## FBpp0070977                 305
## FBpp0309146                4829
## FBpp0073154                  39
## FBpp0297431                 127
## FBpp0111550                 477
## FBpp0297369                3264
## FBpp0292362                 135
## FBpp0076391                1206
## FBpp0076393                2055
## FBpp0076083                 617
## FBpp0085866                4834
## FBpp0309317                3578
## FBpp0085865                5718
## FBpp0085869                4553
## FBpp0310296                 402
## FBpp0308526                 580
## FBpp0308430                 286
## FBpp0079802                 955
## FBpp0301708                 133
## FBpp0291801                   9
## FBpp0291802                  93
## FBpp0081671                3917
## FBpp0311386                4557
## FBpp0304902                2854
## FBpp0310418                5965
## FBpp0307951                  40
## FBpp0077750                 349
## FBpp0304134                 219
## FBpp0296949                3941
## FBpp0075348                 539
## FBpp0088490                4009
## FBpp0088945                4128
## FBpp0071069                 454
## FBpp0303564               22952
## FBpp0072601                 524
## FBpp0304008                 655
## FBpp0307166                 989
## FBpp0305415                 394
## FBpp0305414                 424
## FBpp0306018                 406
## FBpp0290840                1509
## FBpp0304449                1064
## FBpp0290720                3345
## FBpp0303277                2045
## FBpp0085646                5155
## FBpp0303276                3307
## FBpp0076329                 162
## FBpp0081619                 392
## FBpp0086640                 755
## FBpp0289951                 759
## FBpp0073067                 181
## FBpp0302570                1876
## FBpp0306848                1277
## FBpp0087585                1670
## FBpp0297612                4506
## FBpp0113036                8336
## FBpp0079390                 613
## FBpp0074650                3129
## FBpp0086349                2292
## FBpp0293835                  23
## FBpp0081956                 204
## FBpp0082107                 229
## FBpp0087010                1417
## FBpp0070874                 370
## FBpp0075088                1299
## FBpp0307982               98332
## FBpp0072975                1489
## FBpp0087712                 415
## FBpp0074285                 733
## FBpp0291370                 294
## FBpp0288668                1783
## FBpp0087926                 136
## FBpp0083321                 440
## FBpp0290811                1740
## FBpp0310381                 723
## FBpp0310380                  78
## FBpp0292261                  31
## FBpp0305198                 322
## FBpp0297444                2254
## FBpp0290693                1931
## FBpp0297443                2257
## FBpp0085367                2487
## FBpp0100031                1413
## FBpp0312456                 254
## FBpp0084716                  20
## FBpp0078142                 344
## FBpp0311638                 100
## FBpp0081810                 110
## FBpp0112110                 385
## FBpp0076076                 407
## FBpp0306019                 841
## FBpp0305179                 221
## FBpp0292727                 407
## FBpp0292725                 397
## FBpp0289282                5820
## FBpp0072638                 706
## FBpp0082350                2757
## FBpp0080823                  88
## FBpp0077964                 891
## FBpp0291724                2466
## FBpp0291723                4363
## FBpp0081322                 705
## FBpp0079524                1751
## FBpp0073256               59392
## FBpp0081601               17525
## FBpp0075940                 163
## FBpp0079033                 320
## FBpp0111941                 761
## FBpp0292379                   4
## FBpp0080489                 283
## FBpp0086316                 124
## FBpp0083855                3106
## FBpp0308734                 162
## FBpp0085064                 599
## FBpp0271781                1103
## FBpp0079253                 181
## FBpp0079155                  29
## FBpp0300436                1339
## FBpp0307652                4582
## FBpp0074220                 646
## FBpp0074686                 390
## FBpp0074685                 394
## FBpp0304337                  41
## FBpp0089048                  39
## FBpp0304336                  93
## FBpp0304339                  58
## FBpp0073801                  22
## FBpp0083988                1096
## FBpp0077208                 419
## FBpp0289724                  16
## FBpp0290027                1411
## FBpp0083505                 902
## FBpp0076349                 511
## FBpp0072277               12068
## FBpp0271909                 258
## FBpp0292671                 790
## FBpp0293411                 735
## FBpp0088182                2023
## FBpp0300576                 248
## FBpp0075344                1778
## FBpp0311424                3286
## FBpp0086956                2392
## FBpp0078652                  69
## FBpp0099650                  64
## FBpp0084341                 533
## FBpp0086718                7498
## FBpp0310303                1001
## FBpp0290493                1730
## FBpp0076346                1004
## FBpp0310079                 194
## FBpp0310199                4368
## FBpp0309007                2338
## FBpp0290728                 785
## FBpp0074831                4283
## FBpp0311237                1089
## FBpp0307146                 233
## FBpp0311705               11016
## FBpp0081324                4061
## FBpp0083453                 104
## FBpp0071077                 704
## FBpp0071078                 176
## FBpp0112020                 262
## FBpp0081813                 166
## FBpp0081404                 763
## FBpp0291063                  33
## FBpp0291065                 863
## FBpp0310104                 893
## FBpp0291064                 338
## FBpp0071138                2485
## FBpp0081234                 249
## FBpp0305104                2420
## FBpp0075122                 308
## FBpp0290613                1841
## FBpp0088350                1091
## FBpp0079247                1815
## FBpp0292965                1314
## FBpp0298283               20278
## FBpp0081563                 303
## FBpp0308010                1174
## FBpp0072318                 137
## FBpp0311843                 595
## FBpp0308437                 725
## FBpp0312545                 826
## FBpp0290714                1555
## FBpp0290713                2321
## FBpp0307414                  95
## FBpp0099512                  29
## FBpp0082228                  60
## FBpp0271771                 648
## FBpp0292977                 667
## FBpp0312390                  14
## FBpp0306663                 468
## FBpp0304979                 912
## FBpp0088184                 126
## FBpp0308319               13868
## FBpp0290729                  92
## FBpp0303602                7277
## FBpp0075170                6343
## FBpp0075171                 989
## FBpp0309151                4948
## FBpp0305606                 788
## FBpp0311588                 213
## FBpp0312012                 504
## FBpp0082754                 113
## FBpp0077337                 302
## FBpp0301023                1887
## FBpp0080690                1260
## FBpp0085457                 903
## FBpp0080719                 924
## FBpp0291768                1595
## FBpp0086484                2825
## FBpp0309269                5477
## FBpp0306655                   0
## FBpp0308322                  14
## FBpp0078543                   7
## FBpp0070255                 762
## FBpp0088519                1948
## FBpp0079874                1484
## FBpp0309262                 575
## FBpp0089259                 439
## FBpp0303140                  47
## FBpp0303512                 739
## FBpp0086531                 685
## FBpp0076342                 588
## FBpp0076341                 566
## FBpp0074558                4388
## FBpp0112125                  49
## FBpp0304596                1154
## FBpp0082923                 407
## FBpp0310259                 595
## FBpp0080676                1222
## FBpp0305099                 348
## FBpp0304216                 221
## FBpp0089221                 799
## FBpp0311249                 917
## FBpp0088498                 257
## FBpp0307026                 633
## FBpp0083163                 234
## FBpp0307645                 267
## FBpp0079573                 768
## FBpp0300655               18715
## FBpp0084585                5168
## FBpp0289762                 702
## FBpp0289758                 737
## FBpp0304847                3336
## FBpp0297624                1837
## FBpp0307553                1989
## FBpp0073651                 439
## FBpp0077713                  30
## FBpp0076438                 329
## FBpp0309417                 801
## FBpp0071144                2873
## FBpp0111874                 173
## FBpp0075118                 764
## FBpp0310589                 252
## FBpp0077220                   5
## FBpp0304993               10833
## FBpp0306197                1218
## FBpp0085216                 460
## FBpp0075916                2613
## FBpp0310202                 371
## FBpp0084955                  22
## FBpp0084956                   7
## FBpp0081650                 768
## FBpp0076408                6991
## FBpp0072461                1805
## FBpp0084252                6087
## FBpp0074595                 429
## FBpp0304939                 162
## FBpp0307855                 432
## FBpp0307451                 500
## FBpp0302864                 531
## FBpp0072874                1638
## FBpp0304148                  32
## FBpp0294043                  64
## FBpp0309242                 108
## FBpp0079944                 531
## FBpp0075458                1028
## FBpp0099768                  53
## FBpp0071905                1209
## FBpp0305296                2245
## FBpp0301223               10191
## FBpp0078935                  16
## FBpp0307930                 421
## FBpp0072629                2516
## FBpp0290774                1927
## FBpp0291126                 333
## FBpp0083897                1418
## FBpp0310636                 962
## FBpp0309502                 206
## FBpp0290395                  91
## FBpp0290862                1912
## FBpp0078606                2468
## FBpp0071254                 281
## FBpp0079663               22064
## FBpp0301123                3124
## FBpp0309427                2282
## FBpp0309943                 756
## FBpp0082187                2113
## FBpp0088538                 990
## FBpp0312497                  30
## FBpp0076326                  22
## FBpp0291183                  40
## FBpp0077654                 868
## FBpp0073113                   0
## FBpp0310948                3937
## FBpp0307565                  56
## FBpp0310187                 665
## FBpp0290522                 613
## FBpp0071973                1094
## FBpp0301840                1162
## FBpp0292193                 354
## FBpp0304305                 449
## FBpp0310504                   8
## FBpp0312538                 118
## FBpp0308486                3795
## FBpp0073365                1133
## FBpp0307854                 149
## FBpp0082065                1625
## FBpp0303087                 155
## FBpp0076824                 438
## FBpp0078363                1189
## FBpp0302626               17564
## FBpp0291679                 860
## FBpp0071936                 286
## FBpp0307842                1130
## FBpp0077507                 787
## FBpp0086054                  12
## FBpp0088592                  15
## FBpp0086657                 100
## FBpp0303934                2592
## FBpp0309866                 165
## FBpp0082719                 861
## FBpp0077133                 603
## FBpp0110463                1286
## FBpp0087645                 453
## FBpp0292924                2083
## FBpp0302662                 589
## FBpp0311373                3915
## FBpp0075713                 885
## FBpp0304398                 873
## FBpp0290828                1063
## FBpp0084175                   5
## FBpp0084174                  40
## FBpp0077868                 559
## FBpp0073672                 487
## FBpp0085713                 666
## FBpp0112042                2157
## FBpp0072863                2120
## FBpp0303359                 323
## FBpp0088009                 967
## FBpp0307642                 746
## FBpp0297360                1951
## FBpp0081048                 383
## FBpp0070865                 892
## FBpp0077362                1829
## FBpp0070584                1762
## FBpp0290710                1133
## FBpp0305007                8708
## FBpp0083611                9661
## FBpp0080704                1161
## FBpp0309669                 582
## FBpp0111766                  49
## FBpp0084115                2449
## FBpp0071072                2935
## FBpp0308636                 576
## FBpp0297518                 340
## FBpp0290035                1576
## FBpp0309825                1976
## FBpp0309714                1306
## FBpp0088493                 171
## FBpp0075166                  46
## FBpp0076879                 357
## FBpp0084580                 213
## FBpp0070155                 382
## FBpp0089216                1469
## FBpp0087040                 468
## FBpp0089293                 133
## FBpp0079565                6445
## FBpp0071503                 435
## FBpp0072557                 537
## FBpp0288575                 989
## FBpp0072590                 258
## FBpp0290366                 699
## FBpp0074971                 551
## FBpp0084739                 129
## FBpp0305453                   8
## FBpp0305532                 739
## FBpp0082801                 677
## FBpp0070719                1485
## FBpp0292262                3563
## FBpp0305108                 318
## FBpp0084742                 440
## FBpp0084350                 292
## FBpp0081153              399313
## FBpp0086469                 610
## FBpp0311431                2245
## FBpp0081355                 926
## FBpp0081756                 813
## FBpp0074426                1192
## FBpp0088506                4850
## FBpp0081682                1334
## FBpp0087711                 181
## FBpp0309915                  10
## FBpp0070635                1952
## FBpp0075519                 340
## FBpp0082832                 688
## FBpp0308474                 229
## FBpp0300972                 297
## FBpp0088877                1221
## FBpp0110197                1139
## FBpp0073236                 387
## FBpp0307917                1278
## FBpp0309076                   2
## FBpp0081672                6793
## FBpp0309189                 135
## FBpp0071449                 385
## FBpp0297411                1015
## FBpp0071752                  97
## FBpp0309391                1174
## FBpp0307248                1350
## FBpp0071194                1152
## FBpp0078334                 127
## FBpp0078333                 348
## FBpp0309029                 792
## FBpp0086650                 474
## FBpp0071462                 315
## FBpp0309429                 562
## FBpp0306887                 237
## FBpp0081082                 249
## FBpp0304769                 581
## FBpp0288449                 210
## FBpp0298003                 863
## FBpp0085308                 767
## FBpp0072782                 267
## FBpp0300747                  74
## FBpp0306127                 378
## FBpp0309847                 337
## FBpp0084759                5664
## FBpp0310071                1613
## FBpp0307446                 352
## FBpp0290865                2514
## FBpp0307388                3516
## FBpp0084748                 627
## FBpp0084750                2568
## FBpp0309852                2620
## FBpp0087918                1108
## FBpp0312411                 543
## FBpp0088901                1397
## FBpp0079514                 374
## FBpp0291851                  83
## FBpp0304445                2585
## FBpp0077125                 296
## FBpp0084635                  22
## FBpp0079666                  51
## FBpp0305602                 868
## FBpp0292906                 210
## FBpp0307720                 504
## FBpp0271936                2284
## FBpp0300820                3984
## FBpp0085072                 137
## FBpp0084839                 214
## FBpp0087094                 749
## FBpp0292077                1187
## FBpp0308932                 664
## FBpp0086340                 764
## FBpp0297608                 410
## FBpp0310748                1306
## FBpp0301183                4377
## FBpp0081613                   7
## FBpp0083244                 493
## FBpp0078472                2327
## FBpp0078469                 824
## FBpp0306930                 798
## FBpp0297182                 707
## FBpp0310684                  55
## FBpp0312292                3585
## FBpp0086247                   6
## FBpp0309068                 147
## FBpp0080048                8470
## FBpp0310139                  65
## FBpp0082545                 187
## FBpp0290896                  21
## FBpp0087084                1209
## FBpp0291294                 556
## FBpp0070431                 726
## FBpp0310902                1431
## FBpp0289079                 480
## FBpp0074770                  55
## FBpp0076467                 919
## FBpp0087483                 374
## FBpp0073659                1282
## FBpp0304173                 946
## FBpp0308303                   3
## FBpp0310396                 901
## FBpp0072956                 313
## FBpp0089217                 891
## FBpp0080191                3024
## FBpp0293371                  37
## FBpp0078829                1128
## FBpp0078827                4442
## FBpp0311332                 268
## FBpp0311333                 193
## FBpp0297498                  82
## FBpp0291374                  97
## FBpp0290947                 163
## FBpp0306506               26171
## FBpp0309978                  63
## FBpp0291853                1891
## FBpp0289912                 221
## FBpp0077884                   4
## FBpp0077883                  56
## FBpp0308405                2787
## FBpp0073040                 128
## FBpp0309612                3760
## FBpp0311683                 390
## FBpp0310032                 134
## FBpp0077029                 658
## FBpp0074113                 340
## FBpp0077738                 423
## FBpp0292775                 114
## FBpp0112025               10022
## FBpp0089139                8503
## FBpp0297171                 923
## FBpp0290315                1104
## FBpp0076331                 325
## FBpp0300973                 298
## FBpp0310009                 559
## FBpp0082223                7949
## FBpp0307187                 468
## FBpp0307186                3721
## FBpp0304407                 469
## FBpp0072899                 473
## FBpp0307805                 221
## FBpp0084244                 846
## FBpp0112160                 250
## FBpp0308524                 634
## FBpp0072453                 454
## FBpp0073902                7683
## FBpp0300790                2219
## FBpp0304832                 332
## FBpp0294034                 341
## FBpp0307426                2627
## FBpp0081896                2914
## FBpp0084486                 211
## FBpp0080495                1813
## FBpp0307142                1792
## FBpp0080496                 556
## FBpp0079454                5017
## FBpp0075396                 261
## FBpp0077873                 152
## FBpp0075398                1069
## FBpp0087732                 130
## FBpp0306672                1923
## FBpp0306632                2273
## FBpp0083179                3415
## FBpp0081055                 601
## FBpp0290679                 793
## FBpp0309338                 930
## FBpp0087780                3271
## FBpp0111956                 101
## FBpp0308815                 860
## FBpp0304976                1763
## FBpp0111689                1520
## FBpp0305859                  29
## FBpp0074836                  76
## FBpp0289451                 414
## FBpp0083935                 415
## FBpp0075248                4753
## FBpp0305488                4804
## FBpp0076488                 869
## FBpp0111729                 862
## FBpp0076486                 764
## FBpp0297602                 857
## FBpp0076413                 817
## FBpp0292356                 436
## FBpp0301966                 106
## FBpp0080392                 354
## FBpp0086590                 775
## FBpp0310569                 322
## FBpp0083225                  97
## FBpp0082624                 154
## FBpp0311839                1595
## FBpp0078425                 232
## FBpp0309629                   0
## FBpp0307200                 840
## FBpp0082070                 126
## FBpp0084348                 533
## FBpp0110299                  44
## FBpp0074162                 666
## FBpp0297167                 646
## FBpp0078402                 227
## FBpp0080017                 813
## FBpp0082352                1761
## FBpp0085616                 932
## FBpp0070373                 367
## FBpp0309962                 399
## FBpp0297204                1804
## FBpp0311440                 136
## FBpp0289046                 138
## FBpp0070094                  70
## FBpp0306551                1768
## FBpp0311451                2673
## FBpp0306904                 684
## FBpp0292890                 750
## FBpp0079182                 461
## FBpp0081434                 473
## FBpp0310065                 829
## FBpp0290409                 786
## FBpp0305003                 183
## FBpp0290541                  62
## FBpp0110215                 591
## FBpp0110216                 583
## FBpp0289622                 296
## FBpp0311786                 576
## FBpp0288892                1655
## FBpp0080237                3557
## FBpp0309425                1062
## FBpp0079519                2664
## FBpp0073717                 787
## FBpp0080293                   4
## FBpp0082359                 154
## FBpp0074713                  61
## FBpp0294007                  79
## FBpp0075937                1516
## FBpp0302564                 996
## FBpp0310151                 185
## FBpp0307613                1283
## FBpp0307614                2024
## FBpp0305795                2035
## FBpp0305929                  12
## FBpp0303428               13110
## FBpp0077343                1644
## FBpp0293153                 410
## FBpp0308330                 184
## FBpp0074538                 410
## FBpp0074246                 300
## FBpp0289116                 406
## FBpp0078927                1156
## FBpp0111580                 272
## FBpp0111581                 111
## FBpp0311993                 132
## FBpp0301283                 554
## FBpp0310514                1060
## FBpp0083814                  64
## FBpp0309252                 890
## FBpp0309251                  27
## FBpp0085546                 482
## FBpp0070899                3777
## FBpp0308555                  42
## FBpp0073199               31967
## FBpp0085617                 600
## FBpp0083879                 792
## FBpp0076533                 445
## FBpp0072033                 810
## FBpp0300997                  37
## FBpp0290635                  13
## FBpp0290637                   2
## FBpp0307448                  91
## FBpp0305528                  36
## FBpp0309676                 497
## FBpp0289404                 258
## FBpp0082120                  93
## FBpp0075441                 857
## FBpp0081493                 437
## FBpp0311646                  71
## FBpp0297092                 124
## FBpp0079847                 242
## FBpp0078372                 744
## FBpp0311009                1790
## FBpp0312439                 114
## FBpp0085365                 290
## FBpp0082597                 114
## FBpp0291071                4116
## FBpp0084210                  43
## FBpp0073103                 877
## FBpp0076132                 361
## FBpp0075935                1461
## FBpp0075934                1076
## FBpp0310637                7854
## FBpp0079527                7069
## FBpp0076784                1078
## FBpp0113091                 648
## FBpp0290169                 256
## FBpp0081886                2323
## FBpp0087154                 301
## FBpp0300236                1050
## FBpp0083713                1081
## FBpp0084689                1482
## FBpp0081473                2032
## FBpp0305823                1998
## FBpp0074206                 559
## FBpp0074205                 524
## FBpp0083413                5016
## FBpp0311406                  97
## FBpp0111853                 111
## FBpp0071938                 574
## FBpp0071937                 662
## FBpp0085912                 514
## FBpp0311802                 520
## FBpp0089422                1384
## FBpp0110316                1142
## FBpp0293227                2861
## FBpp0290866                  65
## FBpp0292613                1623
## FBpp0292910                2216
## FBpp0301134                   2
## FBpp0288437                   2
## FBpp0312193                 110
## FBpp0076936                 163
## FBpp0305555                 159
## FBpp0301128                 356
## FBpp0301125                 109
## FBpp0300161                2372
## FBpp0084452                2690
## FBpp0306270                2438
## FBpp0307572                  85
## FBpp0303191                 887
## FBpp0077893                 471
## FBpp0083885                  80
## FBpp0290185                 118
## FBpp0080931                 595
## FBpp0078296                 925
## FBpp0311643               12473
## FBpp0076186                  17
## FBpp0110083                   0
## FBpp0297227                 690
## FBpp0311995                3765
## FBpp0312093                4879
## FBpp0073588                 879
## FBpp0304873                 216
## FBpp0302777                5481
## FBpp0070330               12019
## FBpp0289493                1899
## FBpp0290057                 219
## FBpp0311034                1569
## FBpp0071746                 514
## FBpp0305326                1357
## FBpp0311768                 103
## FBpp0311803                 284
## FBpp0271832                 487
## FBpp0089163                1092
## FBpp0305662                 832
## FBpp0305522                 322
## FBpp0290099                  23
## FBpp0073128                 103
## FBpp0070110                  36
## FBpp0290755                6914
## FBpp0078678                8988
## FBpp0078679                6641
## FBpp0078680                6584
## FBpp0290753               10252
## FBpp0078674               10237
## FBpp0311450                3449
## FBpp0083227                1685
## FBpp0070652                1904
## FBpp0087744                 386
## FBpp0083842                 331
## FBpp0310511                 103
## FBpp0074671                1603
## FBpp0297103                2076
## FBpp0307144                2035
## FBpp0082314                1611
## FBpp0080747                 500
## FBpp0307178                1153
## FBpp0079254                 147
## FBpp0079255                 660
## FBpp0078797                  12
## FBpp0312226               15373
## FBpp0076718                 448
## FBpp0311478                1413
## FBpp0070751                1388
## FBpp0301702                1006
## FBpp0071288                1024
## FBpp0084410                1270
## FBpp0290701                2256
## FBpp0304333                  37
## FBpp0304672                 325
## FBpp0312207                 616
## FBpp0078116                 316
## FBpp0081608                1293
## FBpp0303579                   7
## FBpp0070730                 205
## FBpp0086858                  22
## FBpp0100135                  52
## FBpp0304159                  19
## FBpp0079425                2040
## FBpp0072612                1419
## FBpp0083514                  98
## FBpp0083658                 979
## FBpp0300948                 835
## FBpp0309228                3696
## FBpp0077899                 296
## FBpp0309275                 118
## FBpp0099795                 802
## FBpp0088033                 163
## FBpp0292920                 633
## FBpp0079347                 629
## FBpp0309192                 135
## FBpp0307570                 804
## FBpp0306011                 906
## FBpp0080118                1730
## FBpp0310753                1551
## FBpp0081780                7429
## FBpp0311199                  63
## FBpp0303185                 141
## FBpp0297586               12476
## FBpp0290579                  61
## FBpp0293283                 540
## FBpp0083436                 629
## FBpp0083399                 289
## FBpp0075276                 905
## FBpp0072803                 570
## FBpp0293860                 423
## FBpp0083275                   8
## FBpp0310640                  52
## FBpp0305834                 132
## FBpp0072562                1863
## FBpp0078726                 238
## FBpp0293533                 510
## FBpp0293530                 442
## FBpp0072906                1111
## FBpp0111528                1508
## FBpp0304182                1292
## FBpp0074559                1434
## FBpp0309263                1742
## FBpp0302955                 613
## FBpp0304613                2139
## FBpp0077447                 443
## FBpp0085763                 663
## FBpp0090954                 442
## FBpp0306434                 368
## FBpp0089292                  97
## FBpp0310440                1781
## FBpp0303715                  56
## FBpp0309362                 967
## FBpp0305199                 172
## FBpp0076330                 126
## FBpp0288869                1989
## FBpp0311911                1726
## FBpp0307983                8292
## FBpp0302536                1044
## FBpp0306001                 406
## FBpp0307977                   4
## FBpp0305775                3596
## FBpp0297531                   3
## FBpp0301158                  37
## FBpp0089033                5757
## FBpp0304835                5783
## FBpp0073620                 850
## FBpp0089405                  69
## FBpp0088792                2555
## FBpp0292532                 348
## FBpp0307163                1352
## FBpp0080012                 386
## FBpp0312437                1425
## FBpp0078642               40429
## FBpp0085657                 111
## FBpp0293054                 531
## FBpp0311547                3855
## FBpp0302962                 303
## FBpp0076520                2841
## FBpp0073714                 474
## FBpp0070708                2105
## FBpp0070707                 106
## FBpp0310379                 602
## FBpp0086501                 705
## FBpp0074343                1515
## FBpp0074342                1752
## FBpp0311662                  16
## FBpp0291423                 227
## FBpp0087398                 519
## FBpp0305083                  97
## FBpp0310796                  71
## FBpp0080945                 292
## FBpp0072475                 717
## FBpp0070798              106831
## FBpp0085809                3958
## FBpp0084176                3772
## FBpp0305448                1578
## FBpp0306905                5600
## FBpp0112097                  53
## FBpp0311408                  68
## FBpp0297108                  74
## FBpp0311758                 615
## FBpp0305422                  48
## FBpp0075680                2432
## FBpp0087544                1157
## FBpp0290074                2128
## FBpp0087535                 786
## FBpp0081130                 179
## FBpp0289747                1802
## FBpp0309398                 624
## FBpp0309396                 198
## FBpp0288679                 686
## FBpp0302956                 576
## FBpp0271928                 159
## FBpp0073823                5013
## FBpp0079924                1898
## FBpp0087855                 317
## FBpp0077925                1987
## FBpp0080817                 862
## FBpp0311225                 589
## FBpp0071217                1488
## FBpp0289916                   4
## FBpp0070006                1924
## FBpp0070355                3907
## FBpp0304913                   7
## FBpp0083613                  94
## FBpp0306504                1021
## FBpp0077119                1254
## FBpp0311657                8743
## FBpp0309217                 360
## FBpp0075794                5592
## FBpp0084386                  13
## FBpp0074401                2851
## FBpp0293060                1651
## FBpp0297302                3537
## FBpp0289196                3183
## FBpp0307245                2105
## FBpp0289199                3314
## FBpp0290285                1578
## FBpp0070041                4091
## FBpp0082985                3962
## FBpp0293233                 763
## FBpp0070811                 736
## FBpp0077867                 423
## FBpp0301152               10381
## FBpp0307374                1026
## FBpp0088044                 579
## FBpp0310230                3194
## FBpp0297726                 130
## FBpp0312446                1585
## FBpp0074915                1637
## FBpp0074917                1552
## FBpp0086664                 301
## FBpp0308325                 644
## FBpp0099896                1359
## FBpp0291730               18370
## FBpp0291729                 433
## FBpp0111759                 228
## FBpp0311113                 326
## FBpp0072567                2093
## FBpp0072568                1994
## FBpp0072685                   1
## FBpp0083686                 852
## FBpp0290295                 756
## FBpp0309318                5188
## FBpp0088429                1021
## FBpp0291866                 669
## FBpp0290196                 709
## FBpp0290195                  64
## FBpp0073247                  40
## FBpp0304109                  25
## FBpp0073246                  72
## FBpp0079042                  15
## FBpp0303230                 280
## FBpp0082592                1119
## FBpp0308326                1815
## FBpp0298033                1032
## FBpp0293201                4434
## FBpp0303033                 342
## FBpp0071049                 499
## FBpp0112712                1142
## FBpp0112713                 723
## FBpp0304252                2480
## FBpp0076892                2626
## FBpp0073893                1158
## FBpp0304930                 146
## FBpp0304921                 765
## FBpp0304920                1226
## FBpp0304314                 290
## FBpp0310097                 373
## FBpp0075360                1093
## FBpp0080045                 485
## FBpp0082055                 217
## FBpp0080915                 530
## FBpp0086116                1672
## FBpp0310439                 152
## FBpp0303068                 848
## FBpp0297935                3506
## FBpp0307449                 157
## FBpp0082935                 196
## FBpp0112331                1827
## FBpp0302987                 100
## FBpp0310714                 537
## FBpp0083082                1902
## FBpp0076385              315796
## FBpp0271876                 214
## FBpp0312324                  78
## FBpp0292158                 391
## FBpp0302934                1222
## FBpp0307685                  89
## FBpp0300173                 386
## FBpp0088654                 355
## FBpp0077520                 228
## FBpp0075038                 860
## FBpp0077503                2695
## FBpp0305130                1742
## FBpp0298330                 276
## FBpp0298328                 144
## FBpp0082036                 144
## FBpp0082034                  30
## FBpp0303610                 147
## FBpp0082035                  33
## FBpp0297342                2654
## FBpp0297339                2762
## FBpp0087069                  44
## FBpp0293018                 489
## FBpp0311368                 297
## FBpp0083958                 261
## FBpp0300568                1939
## FBpp0071192               15294
## FBpp0087058               11786
## FBpp0298037                   0
## FBpp0309260                 635
## FBpp0072148                 179
## FBpp0072952                2223
## FBpp0079770                 636
## FBpp0305544                 172
## FBpp0087353                1081
## FBpp0305265                6706
## FBpp0293000                2196
## FBpp0307992                4036
## FBpp0292404                1745
## FBpp0293465                1510
## FBpp0307387                6995
## FBpp0087120                 291
## FBpp0297105                4970
## FBpp0297107                 164
## FBpp0304785                 598
## FBpp0080319                 645
## FBpp0304830                6020
## FBpp0085870                1655
## FBpp0302018                1271
## FBpp0312568                1807
## FBpp0073848                1227
## FBpp0305206                 937
## FBpp0310741                4534
## FBpp0292340                 303
## FBpp0078920                 677
## FBpp0307630                 372
## FBpp0310809                 332
## FBpp0292773                 351
## FBpp0307564                 395
## FBpp0305270                 647
## FBpp0072006                 119
## FBpp0079156                 126
## FBpp0080862                 167
## FBpp0070216                1576
## FBpp0309594                 668
## FBpp0076654                1688
## FBpp0304343               16364
## FBpp0308578                 130
## FBpp0293051               64683
## FBpp0078722                 305
## FBpp0078723                 259
## FBpp0290848                 849
## FBpp0290847                 131
## FBpp0081509                3771
## FBpp0111835                 817
## FBpp0078268               17776
## FBpp0073996                2641
## FBpp0311495                1208
## FBpp0073674                1088
## FBpp0087406                  41
## FBpp0305988                 273
## FBpp0071406                 445
## FBpp0312156                 427
## FBpp0291523                 335
## FBpp0300977                 101
## FBpp0084089                2813
## FBpp0085117                 273
## FBpp0111932                 270
## FBpp0077316                 768
## FBpp0288548                2218
## FBpp0288731                1343
## FBpp0078648                 245
## FBpp0297777                 206
## FBpp0303837                 245
## FBpp0076549                 691
## FBpp0300839                   3
## FBpp0099631                 231
## FBpp0306416                 386
## FBpp0099632                 383
## FBpp0301556                 380
## FBpp0075675                  96
## FBpp0081447                 117
## FBpp0082768                   5
## FBpp0300995                  24
## FBpp0078182                 275
## FBpp0288833                 334
## FBpp0080280               21931
## FBpp0077262                3646
## FBpp0290850                  14
## FBpp0300975                  59
## FBpp0308852                2802
## FBpp0306093                1589
## FBpp0086732                 114
## FBpp0300976                 978
## FBpp0110523                1170
## FBpp0110483                 134
## FBpp0073449               16020
## FBpp0071298                 335
## FBpp0075693                1426
## FBpp0305742                 197
## FBpp0075239                1083
## FBpp0085765                  51
## FBpp0308364                 150
## FBpp0071279               10744
## FBpp0301591                   7
## FBpp0088357                 353
## FBpp0309307                 489
## FBpp0074026                 651
## FBpp0088013                 505
## FBpp0290003                 519
## FBpp0297426                 857
## FBpp0071700                 204
## FBpp0070224                 530
## FBpp0074788                2347
## FBpp0311541                1530
## FBpp0082664                2903
## FBpp0302792                1701
## FBpp0305673                1203
## FBpp0071115                5923
## FBpp0078166                5496
## FBpp0309959                 875
## FBpp0070831                1755
## FBpp0071418                2128
## FBpp0300807                2154
## FBpp0075115                1418
## FBpp0079819                1241
## FBpp0082550                1202
## FBpp0297136                  14
## FBpp0305313                5100
## FBpp0086941                2586
## FBpp0311473                1102
## FBpp0079685                 361
## FBpp0071212               10520
## FBpp0297401                 582
## FBpp0309997                1286
## FBpp0085715                1595
## FBpp0293217                2093
## FBpp0074936                 607
## FBpp0084907                6531
## FBpp0088054                 211
## FBpp0083932                 693
## FBpp0070876                 230
## FBpp0291029                   5
## FBpp0077302                 297
## FBpp0309145                 258
## FBpp0078358                 375
## FBpp0088561                2678
## FBpp0293844                 908
## FBpp0293837                 445
## FBpp0308484                 261
## FBpp0309944                 642
## FBpp0077022                 221
## FBpp0077099                 334
## FBpp0088818                 667
## FBpp0084213                 213
## FBpp0311052                1957
## FBpp0291660                 737
## FBpp0309024                1179
## FBpp0289975                 942
## FBpp0082989               21135
## FBpp0080125                 915
## FBpp0080124                 896
## FBpp0084043                 188
## FBpp0075096                1161
## FBpp0074104                1344
## FBpp0082541                1982
## FBpp0079549                3272
## FBpp0072239                1775
## FBpp0290893               67350
## FBpp0304044                7539
## FBpp0289202                 427
## FBpp0307782                 202
## FBpp0305102                2999
## FBpp0085589                 299
## FBpp0310242                1843
## FBpp0071849                 598
## FBpp0071269                1263
## FBpp0078843                 564
## FBpp0307972                 354
## FBpp0082855                 130
## FBpp0076485                 656
## FBpp0300612              227657
## FBpp0074568                1888
## FBpp0309288                4471
## FBpp0082148                 498
## FBpp0303793                3517
## FBpp0086939                1210
## FBpp0082813                1388
## FBpp0084852              135214
## FBpp0084729                 689
## FBpp0075186                 567
## FBpp0306920                1143
## FBpp0303568                 807
## FBpp0303566                1041
## FBpp0303567                 877
## FBpp0303572                 135
## FBpp0079844                  46
## FBpp0304186                 275
## FBpp0079699                 725
## FBpp0311506                  93
## FBpp0303939                 375
## FBpp0297063                 279
## FBpp0099825                 286
## FBpp0309030               10570
## FBpp0289825                 683
## FBpp0310679                 127
## FBpp0290777                 232
## FBpp0306931                4697
## FBpp0306932                4980
## FBpp0111703                9430
## FBpp0079621                1627
## FBpp0079620                1541
## FBpp0099403                 379
## FBpp0305725                 240
## FBpp0294008                1436
## FBpp0073262                1563
## FBpp0312024                1080
## FBpp0305599                1099
## FBpp0308632                  90
## FBpp0310095                 921
## FBpp0081255                 538
## FBpp0311801                 281
## FBpp0088552                 334
## FBpp0088549                 109
## FBpp0088550                 320
## FBpp0088551                  59
## FBpp0311398                3498
## FBpp0297162                 645
## FBpp0289985                 595
## FBpp0289573                 140
## FBpp0074230                1946
## FBpp0074365                2822
## FBpp0311684                  75
## FBpp0074145                 171
## FBpp0074144                 166
## FBpp0078189                 195
## FBpp0309136                5172
## FBpp0312067                 364
## FBpp0311350               17966
## FBpp0071676                 553
## FBpp0077724                 346
## FBpp0082198                3709
## FBpp0304357                 998
## FBpp0076244                3406
## FBpp0290592                  73
## FBpp0302852               15051
## FBpp0084974                 717
## FBpp0307415                1678
## FBpp0082030                  31
## FBpp0300571                 347
## FBpp0310367                   7
## FBpp0289826                 273
## FBpp0082642                 585
## FBpp0305826                 115
## FBpp0086416                1844
## FBpp0099994                  85
## FBpp0303192                 118
## FBpp0078161                 450
## FBpp0304320                  41
## FBpp0083033                2226
## FBpp0073935                  46
## FBpp0304105                 127
## FBpp0309250                  96
## FBpp0310006                1668
## FBpp0071624                1371
## FBpp0301213                1864
## FBpp0084688                 689
## FBpp0291670                 841
## FBpp0289873                 432
## FBpp0306680                   2
## FBpp0080559                 281
## FBpp0080560                 324
## FBpp0082061                 139
## FBpp0300539                1901
## FBpp0099784                1971
## FBpp0070202                 631
## FBpp0073530                 304
## FBpp0075120                 314
## FBpp0305621                 286
## FBpp0113023                  35
## FBpp0088471                1070
## FBpp0304171                3973
## FBpp0089164                6599
## FBpp0309972                 687
## FBpp0297436                  71
## FBpp0075677                1743
## FBpp0305943                 589
## FBpp0077012                2070
## FBpp0305095               28383
## FBpp0072781                2381
## FBpp0081600                 451
## FBpp0081310                 660
## FBpp0311267                1187
## FBpp0073976                 846
## FBpp0112329                4990
## FBpp0074702                5303
## FBpp0310405                   8
## FBpp0308016                4955
## FBpp0088478                  74
## FBpp0304476                1366
## FBpp0292427                6959
## FBpp0301586                 251
## FBpp0071277                1510
## FBpp0078450                1194
## FBpp0308997                 214
## FBpp0083007                 179
## FBpp0078404                 840
## FBpp0089179                1336
## FBpp0305722                 222
## FBpp0297580                1024
## FBpp0086016                 570
## FBpp0309558                 261
## FBpp0089192                1048
## FBpp0085524                1615
## FBpp0308657                 301
## FBpp0308658                 253
## FBpp0085694                6614
## FBpp0309999                   0
## FBpp0087756                 249
## FBpp0072177               55091
## FBpp0085720              254380
## FBpp0309138                 664
## FBpp0309137                1701
## FBpp0087673                   8
## FBpp0312008                1962
## FBpp0077073                 216
## FBpp0077074                 229
## FBpp0305363                 441
## FBpp0111556                  14
## FBpp0293539                 476
## FBpp0077625                 533
## FBpp0074847                 191
## FBpp0075498                5625
## FBpp0311978               12916
## FBpp0099868                 304
## FBpp0077896                 372
## FBpp0290552                 876
## FBpp0304691                1128
## FBpp0292227                1062
## FBpp0309795                 428
## FBpp0073135                 197
## FBpp0110561                 396
## FBpp0082798                 464
## FBpp0076375                1106
## FBpp0303522                 215
## FBpp0079757                1869
## FBpp0084926                 103
## FBpp0084678                 162
## FBpp0291398                4529
## FBpp0302943                  30
## FBpp0291373                  50
## FBpp0099801                  77
## FBpp0302948                 104
## FBpp0301784                  72
## FBpp0075952                1678
## FBpp0075785                1583
## FBpp0075071                 140
## FBpp0087646                 621
## FBpp0081840                 635
## FBpp0072977                1001
## FBpp0309742                  62
## FBpp0310070               23297
## FBpp0293997                1117
## FBpp0086427                 106
## FBpp0111676                 131
## FBpp0071536                  63
## FBpp0290444                 439
## FBpp0074453                 133
## FBpp0081280                 910
## FBpp0307212                1298
## FBpp0076461                 360
## FBpp0073874                 304
## FBpp0310509                 510
## FBpp0307038                 872
## FBpp0301792                   0
## FBpp0086503                 988
## FBpp0086669                2566
## FBpp0300754                1703
## FBpp0288750                 317
## FBpp0309393                 131
## FBpp0309394                 321
## FBpp0292044                 186
## FBpp0311694                 843
## FBpp0300207                  47
## FBpp0306617                8701
## FBpp0081412                1530
## FBpp0311120                  31
## FBpp0089099                  31
## FBpp0306427                1484
## FBpp0297861                 618
## FBpp0312556                  87
## FBpp0072175                 549
## FBpp0072176                 619
## FBpp0309841                2773
## FBpp0081010                1067
## FBpp0305300                1514
## FBpp0081031                2463
## FBpp0297159                1666
## FBpp0311540                 270
## FBpp0074076                1021
## FBpp0077206                2563
## FBpp0309042                 498
## FBpp0302561                 367
## FBpp0077193                 225
## FBpp0084662                2793
## FBpp0071256                 680
## FBpp0306949                 624
## FBpp0087514                4035
## FBpp0081004                  64
## FBpp0305755                 390
## FBpp0075275                1160
## FBpp0082569                 892
## FBpp0079305                6158
## FBpp0305784                 889
## FBpp0075267                 443
## FBpp0112507                1156
## FBpp0073359                 176
## FBpp0306285                 104
## FBpp0309690                 989
## FBpp0076990                 645
## FBpp0310005                  15
## FBpp0077624                 181
## FBpp0288476                4389
## FBpp0086114                1142
## FBpp0297907               10631
## FBpp0305594                1073
## FBpp0311378                 217
## FBpp0070677                2370
## FBpp0301227                1405
## FBpp0304278                 109
## FBpp0304734               10017
## FBpp0304802               12888
## FBpp0311905                 487
## FBpp0111655                1104
## FBpp0308545                1067
## FBpp0310067                 227
## FBpp0076328               11098
## FBpp0293240                 374
## FBpp0070997                 340
## FBpp0099836                 309
## FBpp0076823                  40
## FBpp0308810                  16
## FBpp0086727                 307
## FBpp0071142                 530
## FBpp0293260                1163
## FBpp0308603                1411
## FBpp0080918                 137
## FBpp0087286                 122
## FBpp0072675                 247
## FBpp0310232               12626
## FBpp0078482                   1
## FBpp0082676                   6
## FBpp0078379                 896
## FBpp0110392                 574
## FBpp0311091                1051
## FBpp0084999                5049
## FBpp0084998                9716
## FBpp0304477               14388
## FBpp0071744                1190
## FBpp0302924                2150
## FBpp0292491                2297
## FBpp0083800                 263
## FBpp0077033                 495
## FBpp0303527                 333
## FBpp0083704                  15
## FBpp0306035                  99
## FBpp0074937                6443
## FBpp0311608                 791
## FBpp0080679                1214
## FBpp0088169                 404
## FBpp0304017                  86
## FBpp0100104                1084
## FBpp0301108                 133
## FBpp0071994                 211
## FBpp0084575                2038
## FBpp0310702                1465
## FBpp0071136                 441
## FBpp0290388                 759
## FBpp0088377                  82
## FBpp0303668                 311
## FBpp0084248                 241
## FBpp0088518                 428
## FBpp0073508                  40
## FBpp0290924                1867
## FBpp0111603                 118
## FBpp0111604                 452
## FBpp0289067                2342
## FBpp0291053                3346
## FBpp0291052                3569
## FBpp0306149                 250
## FBpp0099909                 321
## FBpp0293077                 467
## FBpp0312375                1921
## FBpp0076470                1641
## FBpp0306699                1038
## FBpp0306701                 403
## FBpp0310302                2381
## FBpp0080638                 791
## FBpp0301780                  95
## FBpp0080165                   1
## FBpp0074009                1002
## FBpp0305995               12338
## FBpp0070363                 630
## FBpp0070364                3011
## FBpp0070362                 637
## FBpp0084495                 914
## FBpp0084494                 714
## FBpp0087762                5747
## FBpp0292240                   2
## FBpp0113017                1404
## FBpp0075761                1166
## FBpp0074733                 428
## FBpp0291675                 389
## FBpp0311748                   3
## FBpp0084776                  98
## FBpp0077252                  57
## FBpp0075533                 194
## FBpp0088417                 294
## FBpp0297991                  71
## FBpp0291528                 267
## FBpp0074845                3935
## FBpp0309048                1809
## FBpp0077042                2058
## FBpp0308589                2258
## FBpp0078689                5139
## FBpp0289611                1386
## FBpp0291030                 450
## FBpp0072908                7734
## FBpp0074122                  48
## FBpp0099826                  73
## FBpp0070402                 500
## FBpp0310622                 705
## FBpp0111786                 325
## FBpp0303394                 542
## FBpp0080257                2555
## FBpp0078326                 922
## FBpp0302575                 653
## FBpp0110132                 419
## FBpp0070862                 197
## FBpp0077107                1486
## FBpp0290722                 472
## FBpp0291457                1020
## FBpp0077147                1416
## FBpp0292780                 360
## FBpp0292781                 243
## FBpp0289186                   0
## FBpp0296944               17598
## FBpp0296943                9405
## FBpp0292316                4422
## FBpp0306777                 203
## FBpp0079946                 545
## FBpp0080387                 369
## FBpp0085851                 231
## FBpp0086057                 479
## FBpp0306074                 806
## FBpp0306682                 954
## FBpp0307737                1202
## FBpp0081476                 738
## FBpp0291789                1507
## FBpp0303374                 575
## FBpp0079912                  55
## FBpp0271815                4737
## FBpp0087047                  32
## FBpp0311524                 418
## FBpp0308533                  25
## FBpp0294031                  25
## FBpp0087979                4465
## FBpp0082579                 315
## FBpp0081733                 299
## FBpp0086104                1912
## FBpp0304462                4445
## FBpp0297663                5163
## FBpp0082593                 129
## FBpp0078006                 151
## FBpp0073648                 310
## FBpp0306724                 533
## FBpp0304888                3713
## FBpp0304885                4858
## FBpp0084731                 135
## FBpp0084730                2913
## FBpp0087135                5510
## FBpp0111944                  61
## FBpp0303775                2260
## FBpp0311670                 356
## FBpp0309872                 688
## FBpp0290012                2075
## FBpp0300406                1397
## FBpp0311421                1128
## FBpp0073005                 273
## FBpp0073310                 746
## FBpp0291479                 742
## FBpp0074460                  74
## FBpp0291264                2759
## FBpp0087366                 537
## FBpp0076338                 517
## FBpp0294000                 908
## FBpp0080114                2387
## FBpp0301163                 228
## FBpp0084441                 685
## FBpp0071672                1380
## FBpp0300801                1804
## FBpp0309600                 759
## FBpp0073474                3312
## FBpp0074919                 490
## FBpp0074742                 936
## FBpp0074708                1636
## FBpp0305392                1903
## FBpp0081453               17576
## FBpp0077427                 153
## FBpp0084254                 651
## FBpp0309851                  89
## FBpp0305604                  14
## FBpp0311260                4344
## FBpp0305786                1276
## FBpp0305785               11719
## FBpp0081507                 826
## FBpp0083099                  21
## FBpp0303771                 555
## FBpp0303930                3704
## FBpp0303931                4445
## FBpp0291072                4386
## FBpp0111686                  39
## FBpp0304507                2640
## FBpp0077174                1386
## FBpp0303628               19350
## FBpp0076121                 567
## FBpp0086049                 102
## FBpp0087272                 534
## FBpp0306755                 871
## FBpp0307450                  46
## FBpp0304330                 522
## FBpp0084807                   2
## FBpp0080287                2103
## FBpp0073391                  63
## FBpp0085157                 372
## FBpp0290176                 719
## FBpp0099488                 311
## FBpp0081725                 775
## FBpp0077788                 426
## FBpp0309407                2414
## FBpp0289098                  40
## FBpp0086771                 610
## FBpp0077121                  56
## FBpp0075456                2619
## FBpp0075859                 510
## FBpp0304907                 275
## FBpp0084264                1111
## FBpp0312580                 477
## FBpp0075790                 968
## FBpp0075791                 946
## FBpp0304695                 228
## FBpp0087277                 636
## FBpp0307452                 684
## FBpp0087692                 157
## FBpp0304736                1833
## FBpp0087221                8648
## FBpp0085514                 350
## FBpp0089175                  79
## FBpp0073488                 153
## FBpp0312046                 351
## FBpp0305549                2049
## FBpp0082205               18647
## FBpp0304812                 597
## FBpp0077617                 379
## FBpp0082213                 115
## FBpp0304650                  64
## FBpp0086665                 839
## FBpp0077909                 439
## FBpp0304054                  47
## FBpp0304053                 172
## FBpp0304049                 152
## FBpp0288814                  40
## FBpp0072673                 841
## FBpp0084099                3516
## FBpp0083452                2467
## FBpp0288456                2084
## FBpp0311069                 926
## FBpp0082540                 369
## FBpp0304632                   1
## FBpp0304623                  31
## FBpp0072825                 178
## FBpp0112050                 951
## FBpp0305280                   1
## FBpp0086084                 269
## FBpp0084165                 435
## FBpp0293620                1563
## FBpp0292197                 219
## FBpp0290456                 577
## FBpp0310597                1304
## FBpp0308219                 499
## FBpp0073267                 359
## FBpp0311871                6212
## FBpp0305239                 306
## FBpp0309465                1780
## FBpp0311517                 192
## FBpp0086373                4582
## FBpp0304067                1753
## FBpp0305757                3531
## FBpp0303193                1896
## FBpp0087986                 698
## FBpp0082270                4852
## FBpp0077204                 994
## FBpp0290837                1320
## FBpp0087719                 988
## FBpp0087720                1049
## FBpp0087718                  42
## FBpp0293114                  31
## FBpp0087721                 866
## FBpp0300426                1094
## FBpp0307227                 153
## FBpp0078598                1927
## FBpp0308990                6940
## FBpp0074616                 143
## FBpp0075523                 746
## FBpp0309601                1462
## FBpp0086020                2422
## FBpp0309015                  92
## FBpp0070160                 385
## FBpp0077962                3583
## FBpp0073100                1143
## FBpp0099998                5651
## FBpp0074400                1484
## FBpp0076201                 260
## FBpp0087368                 340
## FBpp0100067                  60
## FBpp0086958                 666
## FBpp0073438                2265
## FBpp0312146                1217
## FBpp0070411                  58
## FBpp0307743                2233
## FBpp0072703                1247
## FBpp0311281                9939
## FBpp0082383                 308
## FBpp0077054                3119
## FBpp0310632                 401
## FBpp0075632                2657
## FBpp0290049                  80
## FBpp0083070                 412
## FBpp0083673                 617
## FBpp0312104               12149
## FBpp0307598                 177
## FBpp0080003                1165
## FBpp0087481                 445
## FBpp0301728                3262
## FBpp0082847                   3
## FBpp0307411                 275
## FBpp0072077                1226
## FBpp0297347               13296
## FBpp0089315                 794
## FBpp0305278                 348
## FBpp0306125                 582
## FBpp0085707                1278
## FBpp0310802                 633
## FBpp0079631                 556
## FBpp0298308                3258
## FBpp0087774                 663
## FBpp0080016                1049
## FBpp0075169                 234
## FBpp0304267                 539
## FBpp0297286                 451
## FBpp0304213                 229
## FBpp0075700                  96
## FBpp0309730                1168
## FBpp0080701                 458
## FBpp0306553                 525
## FBpp0088734                  19
## FBpp0304855                 108
## FBpp0304854                 249
## FBpp0310015                 447
## FBpp0072022                 547
## FBpp0081737                 213
## FBpp0072540                  28
## FBpp0303106                 198
## FBpp0308767                 197
## FBpp0309884                  62
## FBpp0305671                2033
## FBpp0290965                   1
## FBpp0076094                  36
## FBpp0312515                 643
## FBpp0080022                 539
## FBpp0071222                 317
## FBpp0292592                 536
## FBpp0080049                2393
## FBpp0082865                 607
## FBpp0070172                2237
## FBpp0304575                 110
## FBpp0297275                 152
## FBpp0077623                 214
## FBpp0088430                 485
## FBpp0086848                 617
## FBpp0071583                 439
## FBpp0291372               15208
## FBpp0070998                 635
## FBpp0294038                3013
## FBpp0073906                 537
## FBpp0305942                2755
## FBpp0073890                 285
## FBpp0293606                 163
## FBpp0083026                 173
## FBpp0071470                1681
## FBpp0303745                 223
## FBpp0081331                 287
## FBpp0081330                5265
## FBpp0289146                 943
## FBpp0303963                 234
## FBpp0079220                 909
## FBpp0311688                1757
## FBpp0305479                1659
## FBpp0071376                 285
## FBpp0070644                  98
## FBpp0085953                 706
## FBpp0113073                 228
## FBpp0070907                 337
## FBpp0083898                4750
## FBpp0083899                2764
## FBpp0072839               31845
## FBpp0311204                4226
## FBpp0290407                 618
## FBpp0077678                 254
## FBpp0306037                1175
## FBpp0304606                1177
## FBpp0305783                 624
## FBpp0310086                4488
## FBpp0300565                4094
## FBpp0306444                1879
## FBpp0074177                  42
## FBpp0311744                1744
## FBpp0113092                 610
## FBpp0305470                  92
## FBpp0074709                 679
## FBpp0309305                 549
## FBpp0083174                 515
## FBpp0304992                  56
## FBpp0078139                1046
## FBpp0303007                3160
## FBpp0305433                 149
## FBpp0271920                 155
## FBpp0085394                1180
## FBpp0289731                1289
## FBpp0290630                2534
## FBpp0073538                1109
## FBpp0078414                 556
## FBpp0309926                 547
## FBpp0089002                 214
## FBpp0100043                5991
## FBpp0082727                 917
## FBpp0070415                 404
## FBpp0080943                1580
## FBpp0085658                3664
## FBpp0075646                1037
## FBpp0077132                 478
## FBpp0082468                 439
## FBpp0311863                 480
## FBpp0302568                1326
## FBpp0071811                1067
## FBpp0082539                8373
## FBpp0111711                3104
## FBpp0297439                2090
## FBpp0297438                2853
## FBpp0307469                 448
## FBpp0088886                2841
## FBpp0079093                 189
## FBpp0073835                1063
## FBpp0306600                1289
## FBpp0112215                1383
## FBpp0310241                 760
## FBpp0290716                  21
## FBpp0310184                 805
## FBpp0111975                 883
## FBpp0087487                 350
## FBpp0070637                1423
## FBpp0301756                 258
## FBpp0289650                 146
## FBpp0309333                 679
## FBpp0311785                 518
## FBpp0304678                 145
## FBpp0304676                 185
## FBpp0110232                2325
## FBpp0078457                1037
## FBpp0078454                2597
## FBpp0111520                 203
## FBpp0111521                 257
## FBpp0309018                 168
## FBpp0077714                1304
## FBpp0290558                1987
## FBpp0084110                 837
## FBpp0085121                 741
## FBpp0071235                 887
## FBpp0311987               40716
## FBpp0112471                3971
## FBpp0302837                  36
## FBpp0083397                2167
## FBpp0110120                1053
## FBpp0090943                1651
## FBpp0311184                 339
## FBpp0079294                 689
## FBpp0312522                 225
## FBpp0082044                2270
## FBpp0305767                 440
## FBpp0289118                 434
## FBpp0305931                2384
## FBpp0311859                1095
## FBpp0304807                 251
## FBpp0293948                2723
## FBpp0309949                 333
## FBpp0112375                 148
## FBpp0310726                 194
## FBpp0303765                 532
## FBpp0297250                 314
## FBpp0078399               18136
## FBpp0308728                 562
## FBpp0308653                3418
## FBpp0303635                 503
## FBpp0079437                7823
## FBpp0085372                 145
## FBpp0306395                2499
## FBpp0083022                 264
## FBpp0293159                 642
## FBpp0293161                 640
## FBpp0311289                 271
## FBpp0293162                 278
## FBpp0071848                1197
## FBpp0070298                2452
## FBpp0083948                4070
## FBpp0306961                 719
## FBpp0077841                 478
## FBpp0301955                1384
## FBpp0307022                 445
## FBpp0296971                 464
## FBpp0075027                3173
## FBpp0308602                 155
## FBpp0078405                 764
## FBpp0077492                 494
## FBpp0076451                1098
## FBpp0078565                 312
## FBpp0080819                   2
## FBpp0303977                  19
## FBpp0310135                 621
## FBpp0086328                8279
## FBpp0073767                 150
## FBpp0070113                 381
## FBpp0087515                1976
## FBpp0082982                 176
## FBpp0310166                5679
## FBpp0290024                  61
## FBpp0300430                  67
## FBpp0300414                  37
## FBpp0070482                 235
## FBpp0080448                 515
## FBpp0075644                  47
## FBpp0289778                  95
## FBpp0293231                  38
## FBpp0080456                   4
## FBpp0077890                 146
## FBpp0304246                 162
## FBpp0074730                 109
## FBpp0075538                 727
## FBpp0305182               87649
## FBpp0084016                5953
## FBpp0071571                4118
## FBpp0071570                3228
## FBpp0304765                 482
## FBpp0110241                1544
## FBpp0091112                  43
## FBpp0084812                1128
## FBpp0312237                3246
## FBpp0311354                 654
## FBpp0309268                 116
## FBpp0071629                1088
## FBpp0079568                2507
## FBpp0307849                 879
## FBpp0297078                3423
## FBpp0081074                1779
## FBpp0083512                   8
## FBpp0310305                 228
## FBpp0310304                 831
## FBpp0075238                2558
## FBpp0308840                 319
## FBpp0078012                 332
## FBpp0304959                2681
## FBpp0075212                2892
## FBpp0080595                 231
## FBpp0297152                 969
## FBpp0071168                 149
## FBpp0310508                 318
## FBpp0082590                  66
## FBpp0088146                 217
## FBpp0291325               19219
## FBpp0087222                   3
## FBpp0087733                 889
## FBpp0293020                   3
## FBpp0099713                 843
## FBpp0111811                 277
## FBpp0305723                 145
## FBpp0086223                 619
## FBpp0072931                 704
## FBpp0070054                 601
## FBpp0309712                 688
## FBpp0307880                 167
## FBpp0307664                 122
## FBpp0307797                1932
## FBpp0311695                 269
## FBpp0081828                1204
## FBpp0077085                1929
## FBpp0082685                6755
## FBpp0082686                9187
## FBpp0298359                 238
## FBpp0302585                9039
## FBpp0311357                3433
## FBpp0088628                 492
## FBpp0088627                1141
## FBpp0081980                 754
## FBpp0074972                2414
## FBpp0305446                 292
## FBpp0083501                1301
## FBpp0305410                 403
## FBpp0111757                 484
## FBpp0305411                 484
## FBpp0086704                 491
## FBpp0078120                 619
## FBpp0305015                 354
## FBpp0305016                1521
## FBpp0289642                1110
## FBpp0296980                1742
## FBpp0072618                 181
## FBpp0072620                  23
## FBpp0306964                 211
## FBpp0074348                3333
## FBpp0085137                 461
## FBpp0072323                1329
## FBpp0078315                 279
## FBpp0301208                 231
## FBpp0090963                  73
## FBpp0308828                  62
## FBpp0077396                 505
## FBpp0304567                 553
## FBpp0089095                 141
## FBpp0310306                1197
## FBpp0310459                2048
## FBpp0300317                5036
## FBpp0084623               25639
## FBpp0075318                1379
## FBpp0075319                1220
## FBpp0291544                  12
## FBpp0304250                 197
## FBpp0084240                2815
## FBpp0307707                 415
## FBpp0311006                 251
## FBpp0083850                 496
## FBpp0293271                 513
## FBpp0071945                 737
## FBpp0086203                 416
## FBpp0311640               18243
## FBpp0310889                1942
## FBpp0308775                 278
## FBpp0308774                 399
## FBpp0082287                  12
## FBpp0289769                 129
## FBpp0078887                 441
## FBpp0303137                 494
## FBpp0292214                6767
## FBpp0076112                 749
## FBpp0303404                 725
## FBpp0083632                 765
## FBpp0301099                8165
## FBpp0301949                7332
## FBpp0306969                1525
## FBpp0303328                3394
## FBpp0298338                  43
## FBpp0086101                 490
## FBpp0304370                 165
## FBpp0303933                  15
## FBpp0289092                 100
## FBpp0079045                 390
## FBpp0086969                 132
## FBpp0088857                1619
## FBpp0088858                1689
## FBpp0291582                1167
## FBpp0112323                 541
## FBpp0305229                 365
## FBpp0307247                 339
## FBpp0308789                1480
## FBpp0086980                 597
## FBpp0085855                   0
## FBpp0087198                   1
## FBpp0301154                1406
## FBpp0309937                1095
## FBpp0087722                1747
## FBpp0311236                   8
## FBpp0076252                  13
## FBpp0297429                 823
## FBpp0078636                1097
## FBpp0290430                  64
## FBpp0290429                 179
## FBpp0077028                 870
## FBpp0111469                1056
## FBpp0074481                1165
## FBpp0085545                 711
## FBpp0086181                 466
## FBpp0086182                 437
## FBpp0071810                   4
## FBpp0072893                1005
## FBpp0309120                  23
## FBpp0072117                1112
## FBpp0079493                  43
## FBpp0087864                1606
## FBpp0307775                2013
## FBpp0089026                1693
## FBpp0089025                2258
## FBpp0311942                2658
## FBpp0074414                   6
## FBpp0072593                 883
## FBpp0310538                 492
## FBpp0078832                  18
## FBpp0308780                1157
## FBpp0084349                1416
## FBpp0079586                3338
## FBpp0309518                 341
## FBpp0075612               20166
## FBpp0079885                 535
## FBpp0079886                1069
## FBpp0301059                   0
## FBpp0311367                   3
## FBpp0088565                7139
## FBpp0076237                2173
## FBpp0081216                1640
## FBpp0303947                3409
## FBpp0293574                  89
## FBpp0083595                 283
## FBpp0311704                4779
## FBpp0081863                 680
## FBpp0088926               74605
## FBpp0088452               53995
## FBpp0084561                 175
## FBpp0289119                   1
## FBpp0305476                 352
## FBpp0271885                 320
## FBpp0293167                 751
## FBpp0310058               27817
## FBpp0070317                 419
## FBpp0290995                 665
## FBpp0271898                1773
## FBpp0087065                 553
## FBpp0304111                 540
## FBpp0304752                1064
## FBpp0304980                 908
## FBpp0072423                1032
## FBpp0304874                   5
## FBpp0306941                 588
## FBpp0089380                 265
## FBpp0088152                1315
## FBpp0088153                1410
## FBpp0310844                  62
## FBpp0078478                 674
## FBpp0079693                   0
## FBpp0310307                7164
## FBpp0297611                 189
## FBpp0076153               24756
## FBpp0309532               25934
## FBpp0290743                  24
## FBpp0077354                  80
## FBpp0310842                3377
## FBpp0084788                 240
## FBpp0113077                 457
## FBpp0309776                 104
## FBpp0087006                 250
## FBpp0086387                 348
## FBpp0079732                 965
## FBpp0087251                 598
## FBpp0085961                1099
## FBpp0082547                 488
## FBpp0086663                 315
## FBpp0297292                 402
## FBpp0311797                1047
## FBpp0073557                 936
## FBpp0305432                1402
## FBpp0305431                1518
## FBpp0072930                 300
## FBpp0073400                5094
## FBpp0289693                 348
## FBpp0080887                  12
## FBpp0086205                 320
## FBpp0083195                1600
## FBpp0080886                  23
## FBpp0308361                 581
## FBpp0086435                1106
## FBpp0099793                2439
## FBpp0082794                   0
## FBpp0307859                3149
## FBpp0303864                1124
## FBpp0077122                 758
## FBpp0112349                2827
## FBpp0307592                1127
## FBpp0112984                2365
## FBpp0311763                3471
## FBpp0290588                9965
## FBpp0303266               27086
## FBpp0291327                3456
## FBpp0084241                2073
## FBpp0080755                  77
## FBpp0308859                 138
## FBpp0073803                   0
## FBpp0300986                 102
## FBpp0290880                  49
## FBpp0088686                 761
## FBpp0088685                  25
## FBpp0289803                  19
## FBpp0304195                4506
## FBpp0110174                1430
## FBpp0076837               24512
## FBpp0076647                2859
## FBpp0300836                1852
## FBpp0100080                4052
## FBpp0079218                1035
## FBpp0078361                 181
## FBpp0077144                1432
## FBpp0082438                 342
## FBpp0085716                1405
## FBpp0074609                1366
## FBpp0071921                  44
## FBpp0082225                 485
## FBpp0088909                   2
## FBpp0310568                 424
## FBpp0303669                 701
## FBpp0086186                 461
## FBpp0303034                 861
## FBpp0087376                 485
## FBpp0087375                 494
## FBpp0298350                   7
## FBpp0075317                1371
## FBpp0072660                1713
## FBpp0297922                   8
## FBpp0077335                  77
## FBpp0302002                 155
## FBpp0271693                 234
## FBpp0312512                 889
## FBpp0303188                 563
## FBpp0074251                2063
## FBpp0089088               14758
## FBpp0077103                 343
## FBpp0310656                 961
## FBpp0086603                5362
## FBpp0080564                1355
## FBpp0303613                 789
## FBpp0090953                1288
## FBpp0306723                1171
## FBpp0090952                2666
## FBpp0083248               15371
## FBpp0072854                  85
## FBpp0303172                4172
## FBpp0288784               21662
## FBpp0306142                  58
## FBpp0311343                 692
## FBpp0309942                 925
## FBpp0297229                 425
## FBpp0078463                7207
## FBpp0293292                 164
## FBpp0304791                 319
## FBpp0306684                 171
## FBpp0293338                 434
## FBpp0305213                 614
## FBpp0304655                 100
## FBpp0309618                 646
## FBpp0305074                 685
## FBpp0309171                 818
## FBpp0290490                1134
## FBpp0070953                1628
## FBpp0077314                2492
## FBpp0079090                 169
## FBpp0290511                 154
## FBpp0085224                 492
## FBpp0077934                1204
## FBpp0306799                2102
## FBpp0078343                 564
## FBpp0082596                3848
## FBpp0087770                  52
## FBpp0291059                1265
## FBpp0289521                1162
## FBpp0300815                 267
## FBpp0111666                 677
## FBpp0304349                1118
## FBpp0311618                 191
## FBpp0289959                   9
## FBpp0291628                  85
## FBpp0293235                1853
## FBpp0304386                 257
## FBpp0290275                 105
## FBpp0089344                  59
## FBpp0309392                1513
## FBpp0110110                  84
## FBpp0289288                 447
## FBpp0311906                1768
## FBpp0112293                 422
## FBpp0112292                2013
## FBpp0310349                 527
## FBpp0076099                 168
## FBpp0307929                  86
## FBpp0071478                1657
## FBpp0292511                5295
## FBpp0303451                 811
## FBpp0073084                 112
## FBpp0073627                 405
## FBpp0077836                5717
## FBpp0304201                3869
## FBpp0074107                 553
## FBpp0310364                 188
## FBpp0080011                  55
## FBpp0076343                 737
## FBpp0305262                  87
## FBpp0075645                  76
## FBpp0290798                 411
## FBpp0293213                 782
## FBpp0078625                  38
## FBpp0312080               39444
## FBpp0303072                 494
## FBpp0306751                 502
## FBpp0085281                  11
## FBpp0293600                 400
## FBpp0099843                5338
## FBpp0310025                1779
## FBpp0305548                 235
## FBpp0071593                2169
## FBpp0111700                 545
## FBpp0306893                1540
## FBpp0071469                   4
## FBpp0081027                 678
## FBpp0086629                1825
## FBpp0303879                 255
## FBpp0303082                 326
## FBpp0072072                2384
## FBpp0079589                 576
## FBpp0307453                 260
## FBpp0077339                  39
## FBpp0293270                 404
## FBpp0304988                6615
## FBpp0301574                  27
## FBpp0081989                 427
## FBpp0304694                 781
## FBpp0305426                 609
## FBpp0307147                 371
## FBpp0082591                 659
## FBpp0080826                 822
## FBpp0310843                 386
## FBpp0303944                 573
## FBpp0083249                  31
## FBpp0303595                   7
## FBpp0294023                 539
## FBpp0294024                 581
## FBpp0077911                 908
## FBpp0080700                 498
## FBpp0306780                 776
## FBpp0312224               40617
## FBpp0099770                2844
## FBpp0071316                 265
## FBpp0304066                9023
## FBpp0078756                4783
## FBpp0312149                1787
## FBpp0083656               10587
## FBpp0271854                 236
## FBpp0084161               11723
## FBpp0301218                   7
## FBpp0083373                1462
## FBpp0083975                 327
## FBpp0078663                2999
## FBpp0082329                 193
## FBpp0086786                 728
## FBpp0070517                  28
## FBpp0291631                 101
## FBpp0085353                 607
## FBpp0085351                 608
## FBpp0079614                 325
## FBpp0078265                4622
## FBpp0075168                2850
## FBpp0081602                 247
## FBpp0074055                1551
## FBpp0075707                7275
## FBpp0303809                3237
## FBpp0086322                 839
## FBpp0074017               31067
## FBpp0307742                2349
## FBpp0307741                 643
## FBpp0310028                6845
## FBpp0075837                  42
## FBpp0289271                 261
## FBpp0071530                1777
## FBpp0305284                3435
## FBpp0087196                1013
## FBpp0087534                1922
## FBpp0081209                 419
## FBpp0073009                 644
## FBpp0308705                 735
## FBpp0083503                 625
## FBpp0080449                4188
## FBpp0078376                6012
## FBpp0072135                3533
## FBpp0082599                  37
## FBpp0073293                 173
## FBpp0291316                 370
## FBpp0080889                 794
## FBpp0305150                  63
## FBpp0309296                  11
## FBpp0309298                   4
## FBpp0297427                 549
## FBpp0076695                 366
## FBpp0304366               10711
## FBpp0302767                3822
## FBpp0073562                3320
## FBpp0303477                 163
## FBpp0303476                 126
## FBpp0088396               11404
## FBpp0110350                  36
## FBpp0110346                  61
## FBpp0297695                  23
## FBpp0110394                  22
## FBpp0088091                  64
## FBpp0288481                  96
## FBpp0305289                1471
## FBpp0071505                  39
## FBpp0071507                  43
## FBpp0290663                8953
## FBpp0305367                  10
## FBpp0306714                1167
## FBpp0306430                 727
## FBpp0071509                3490
## FBpp0302593                 186
## FBpp0077167                1935
## FBpp0113027                 129
## FBpp0304593               19777
## FBpp0306438                1749
## FBpp0308267                4137
## FBpp0304638                5496
## FBpp0304642                5090
## FBpp0309276                 683
## FBpp0306903                7979
## FBpp0294020                1876
## FBpp0087126                 325
## FBpp0075323                   1
## FBpp0301732                1430
## FBpp0292147                1007
## FBpp0307963                  97
## FBpp0303370                1212
## FBpp0111712                 281
## FBpp0070791                 197
## FBpp0086110                1602
## FBpp0085449                 131
## FBpp0311537                5584
## FBpp0074717                  99
## FBpp0310174                 572
## FBpp0298370                  48
## FBpp0310631                 978
## FBpp0087347                 377
## FBpp0072151                2976
## FBpp0080521                1825
## FBpp0312423                   8
## FBpp0079964                  44
## FBpp0073643                  26
## FBpp0306704                  16
## FBpp0309257                1278
## FBpp0099426                 676
## FBpp0303989                 279
## FBpp0303990                 367
## FBpp0309988                 424
## FBpp0304594                1067
## FBpp0079702                  61
## FBpp0079183                 250
## FBpp0086904                1322
## FBpp0311996                  30
## FBpp0080774                   6
## FBpp0072016                  31
## FBpp0079780                1963
## FBpp0083678                  76
## FBpp0073016                 778
## FBpp0297442                3266
## FBpp0305318                1906
## FBpp0303571                 373
## FBpp0113010                  76
## FBpp0307736                 454
## FBpp0307747                 151
## FBpp0311844                1131
## FBpp0083168                 223
## FBpp0312381                 777
## FBpp0306611                4157
## FBpp0310682                1063
## FBpp0076363                1760
## FBpp0311505                  84
## FBpp0311562                 779
## FBpp0304270                 829
## FBpp0305302                6152
## FBpp0305303                3392
## FBpp0080408                 472
## FBpp0080256                 280
## FBpp0297243                 996
## FBpp0082250                  26
## FBpp0306911                2664
## FBpp0312031                2351
## FBpp0303895                   6
## FBpp0100187               38614
## FBpp0087031                   1
## FBpp0309567                1556
## FBpp0076183                 193
## FBpp0298026                 491
## FBpp0082767                 545
## FBpp0071259                 532
## FBpp0306412                1757
## FBpp0086643                   6
## FBpp0306740                4594
## FBpp0073750                 371
## FBpp0308667                 232
## FBpp0303465                1352
## FBpp0077210                1161
## FBpp0306846                   8
## FBpp0310039                1287
## FBpp0081592                2128
## FBpp0075923                  40
## FBpp0082068                  32
## FBpp0082600                  25
## FBpp0079399                 250
## FBpp0309930                2025
## FBpp0078385                   0
## FBpp0084162                 456
## FBpp0303849                 604
## FBpp0310826                1050
## FBpp0088080                1843
## FBpp0310275                   2
## FBpp0074825                6369
## FBpp0071732                  11
## FBpp0086096                 561
## FBpp0083371               21317
## FBpp0070594                   1
## FBpp0078319                 832
## FBpp0312482                  17
## FBpp0087236                  22
## FBpp0080532                 737
## FBpp0306146                  18
## FBpp0071275                 841
## FBpp0076732                  21
## FBpp0071892                4199
## FBpp0082867                2317
## FBpp0081659                5694
## FBpp0304538               11341
## FBpp0307011                 820
## FBpp0082996                 667
## FBpp0089177                  61
## FBpp0071940                2370
## FBpp0303003                1237
## FBpp0075833                 385
## FBpp0300826                 505
## FBpp0303827                5018
## FBpp0306442                6220
## FBpp0308324               15333
## FBpp0085917               69936
## FBpp0305067                2726
## FBpp0306023                 478
## FBpp0305158                 179
## FBpp0307150                  20
## FBpp0110455                   4
## FBpp0305520                1780
## FBpp0076459                 132
## FBpp0085065                1809
## FBpp0082121                2162
## FBpp0073110                  25
## FBpp0088021                3166
## FBpp0086098               11336
## FBpp0084782                 197
## FBpp0084842                1382
## FBpp0082984                 763
## FBpp0087647                 417
## FBpp0099972                3240
## FBpp0080622                 724
## FBpp0087092                 360
## FBpp0075466                 126
## FBpp0290948                 340
## FBpp0087941                 179
## FBpp0070359                   1
## FBpp0078894                 178
## FBpp0303176                  53
## FBpp0305308                1464
## FBpp0071145                 136
## FBpp0084989                 764
## FBpp0070143               19390
## FBpp0073458                1044
## FBpp0085430                 952
## FBpp0084172                 361
## FBpp0112128                 180
## FBpp0293236                7011
## FBpp0305603                7154
## FBpp0291704                3380
## FBpp0072083                2273
## FBpp0307562                1332
## FBpp0298351                2638
## FBpp0070875                   0
## FBpp0072518                 439
## FBpp0084050                 301
## FBpp0080659                6858
## FBpp0307999                  72
## FBpp0304005                2568
## FBpp0078929                 739
## FBpp0079470                 320
## FBpp0074092                  13
## FBpp0310415                   3
## FBpp0305562               19810
## FBpp0305250                  47
## FBpp0081390                  16
## FBpp0076723                 509
## FBpp0089108                   6
## FBpp0085122                 141
## FBpp0099934                1148
## FBpp0071897                 835
## FBpp0072129                  92
## FBpp0304566                2223
## FBpp0309606                3565
## FBpp0079233                5985
## FBpp0086627                9182
## FBpp0301157                1899
## FBpp0307729                 575
## FBpp0100180               48139
## FBpp0071847                 811
## FBpp0111920                4538
## FBpp0070860                 521
## FBpp0087118                 912
## FBpp0079495               10669
## FBpp0087367                2484
## FBpp0304236                1002
## FBpp0077247                   6
## FBpp0084959               19521
## FBpp0087013                 689
## FBpp0307127                1683
## FBpp0289480                 662
## FBpp0309705                1752
## FBpp0083124                 395
## FBpp0304919                   4
## FBpp0080715                   6
## FBpp0310904                 120
## FBpp0290000                1007
## FBpp0112156                  38
## FBpp0289181                 346
## FBpp0084950                6219
## FBpp0070817                  37
## FBpp0086701               34557
## FBpp0072144                1116
## FBpp0310558                 890
## FBpp0074662                5475
## FBpp0311889                 245
## FBpp0071223                 304
## FBpp0305858               19374
## FBpp0079979                2486
## FBpp0311474                5819
## FBpp0087346                 418
## FBpp0310165                 772
## FBpp0306592                 283
## FBpp0309710                  11
## FBpp0306426                1159
## FBpp0087463                 205
## FBpp0303030                 628
## FBpp0087870                4181
## FBpp0305141                3376
## FBpp0291643                1065
## FBpp0078061                  23
## FBpp0310943                 400
## FBpp0307559                   2
## FBpp0085875                 355
## FBpp0074191                 394
## FBpp0310331                 359
## FBpp0083687                  39
## FBpp0076686                 334
## FBpp0307934                 243
## FBpp0099679                 922
## FBpp0311983                2241
## FBpp0297140                1497
## FBpp0070703                 155
## FBpp0074863                 119
## FBpp0308792                 332
## FBpp0070295                 246
## FBpp0083415                1623
## FBpp0077308                6383
## FBpp0075970                 191
## FBpp0084774                3587
## FBpp0079375                1488
## FBpp0311396                 793
## FBpp0311481               12600
## FBpp0310769                  40
## FBpp0304189                2116
## FBpp0072941                 653
## FBpp0085718                   1
## FBpp0083843                2755
## FBpp0079324                2699
## FBpp0298273                  75
## FBpp0073316                 817
## FBpp0081245                 413
## FBpp0312000                3387
## FBpp0302563                 489
## FBpp0305750                 539
## FBpp0077735                 845
## FBpp0076704                  29
## FBpp0085690                1907
## FBpp0113033                  27
## FBpp0072463                 348
## FBpp0082735                5648
## FBpp0079219                 884
## FBpp0300512               16633
## FBpp0087734                 682
## FBpp0070716                 951
## FBpp0311458               18289
## FBpp0306837               24963
## FBpp0081350                 206
## FBpp0077145                 246
## FBpp0308582                 669
## FBpp0083962                 242
## FBpp0086468                2800
## FBpp0089414                 169
## FBpp0072691                1667
## FBpp0075349                2629
## FBpp0071497                2144
## FBpp0078246              224738
## FBpp0290448                1193
## FBpp0080120                  65
## FBpp0081371                5523
## FBpp0309989                 903
## FBpp0081372                   2
## FBpp0292398                2782
## FBpp0306090                  15
## FBpp0309009                   1
## FBpp0292508                 454
## FBpp0078447                1606
## FBpp0083799                 328
## FBpp0305395                  76
## FBpp0070476                 211
## FBpp0082137                 763
## FBpp0071748                1578
## FBpp0099824               12658
## FBpp0305495                 274
## FBpp0073430                2630
## FBpp0311691                3631
## FBpp0312315                 299
## FBpp0071846               19056
## FBpp0290696                1091
## FBpp0084505                   5
## FBpp0310390                7812
## FBpp0311384                 838
## FBpp0080628                  92
## FBpp0070949                 577
## FBpp0074909                1471
## FBpp0085902                1288
## FBpp0077538                 364
## FBpp0291744                  10
## FBpp0311526                 217
## FBpp0083928                2919
## FBpp0305267                3267
## FBpp0071451                2743
## FBpp0291478                 945
## FBpp0308496                4324
## FBpp0080648                 995
## FBpp0302815                 219
## FBpp0087969                2603
## FBpp0079041                1338
## FBpp0309738                5520
## FBpp0312179                 176
## FBpp0300667                2008
## FBpp0301600                4062
## FBpp0309234                2065
## FBpp0292605                   9
## FBpp0309201               21921
## FBpp0309477                3908
## FBpp0310207                1276
## FBpp0306002                1089
## FBpp0077806                 204
## FBpp0083451                1187
## FBpp0311966                   6
## FBpp0086465                2014
## FBpp0306603               24835
## FBpp0303937                 840
## FBpp0076098                 459
## FBpp0072349                  12
## FBpp0304055                   1
## FBpp0305777                 463
## FBpp0087957                  60
## FBpp0311728                   9
## FBpp0085166               36267
## FBpp0084911                2662
## FBpp0081867                  48
## FBpp0305959                 489
## FBpp0311265                  22
## FBpp0085260                1340
## FBpp0075034                2129
## FBpp0113056                 410
## FBpp0079492                   3
## FBpp0305334                2240
## FBpp0081860                  71
## FBpp0310433                2205
## FBpp0084528                 403
## FBpp0075382                2866
## FBpp0312205                1296
## FBpp0306036                1238
## FBpp0087859                1302
## FBpp0308731                1678
## FBpp0076545                 155
## FBpp0304646                1709
## FBpp0304645                3799
## FBpp0305717               12441
## FBpp0086269               18892
## FBpp0304381                  22
## FBpp0312542                  44
## FBpp0307760                 939
## FBpp0289972                 751
## FBpp0086235                   1
## FBpp0305836                 992
## FBpp0311779                9226
## FBpp0081879                1792
## FBpp0072146                 381
## FBpp0074693                   0
## FBpp0087241               11941
## FBpp0302861                 468
## FBpp0077933                   0
## FBpp0306868                1361
## FBpp0301709                   0
## FBpp0311461                4481
## FBpp0112463                 256
## FBpp0085703                4440
## FBpp0291141                   2
## FBpp0112271                   0
## FBpp0078701                   0
## FBpp0080335                 732
## FBpp0087227                1097
## FBpp0071825                 197
## FBpp0297872                   0
## FBpp0309036                 368
## FBpp0086841                 580
## FBpp0111805                 534
## FBpp0083502                 972
## FBpp0075854                1474
## FBpp0110410                 950
## FBpp0072029                  28
## FBpp0077839                 425
## FBpp0071981                   6
## FBpp0086795                 731
## FBpp0073974                1117
## FBpp0311639                 249
## FBpp0070306                 416
## FBpp0292879                  35
## FBpp0078416               17068
## FBpp0312189                 861
## FBpp0074513                1262
## FBpp0306392                 561
## FBpp0083972                 525
## FBpp0292258                 279
## FBpp0079892                2218
## FBpp0113050                   8
## FBpp0073074                   8
## FBpp0300658                  16
## FBpp0086334                  40
## FBpp0085560                  28
## FBpp0099382                   6
## FBpp0290229                 146
## FBpp0072687               25746
## FBpp0073805                 285
## FBpp0076655                   5
## FBpp0308558                1005
## FBpp0083411                 500
## FBpp0071703                2586
## FBpp0071794               18926
## FBpp0082953                 300
## FBpp0076134                4881
## FBpp0300391                1080
## FBpp0082877                1578
## FBpp0073058                 254
## FBpp0079550                1162
## FBpp0306730                  60
## FBpp0309239                 125
## FBpp0086067                2372
## FBpp0075676                 443
## FBpp0113028                   2
## FBpp0290083                4596
## FBpp0310192                1873
## FBpp0081744                 364
## FBpp0076643                3483
## FBpp0091107                  89
## FBpp0301542                   0
## FBpp0086666                  24
## FBpp0306203                 184
## FBpp0311405                 186
## FBpp0292380                 216
## FBpp0308546                 209
## FBpp0311452               22218
## FBpp0088502                 146
## FBpp0305517                 145
## FBpp0071476                 959
## FBpp0076142                 909
## FBpp0304016                 335
## FBpp0110438                 771
## FBpp0079352                 256
## FBpp0073459                 336
## FBpp0310669                   5
## FBpp0308423                  37
## FBpp0309555                 194
## FBpp0304126                  21
## FBpp0291626                 576
## FBpp0310050                   0
## FBpp0071587                2615
## FBpp0085636                1309
## FBpp0110260                 620
## FBpp0100141                1999
## FBpp0303962                1485
## FBpp0083028                 183
## FBpp0072458                  99
## FBpp0082758                 999
## FBpp0312460                 451
## FBpp0297522               20161
## FBpp0070058               17936
## FBpp0086844                2028
## FBpp0304101                2244
## FBpp0081205                  24
## FBpp0081302                 588
## FBpp0311550                2345
## FBpp0077189                  20
## FBpp0300893                2850
## FBpp0082462                 731
## FBpp0088522                1664
## FBpp0075043                1174
## FBpp0070924                 274
## FBpp0112365                  98
## FBpp0077004                 833
## FBpp0070244                 354
## FBpp0312506                 434
## FBpp0312112                  12
## FBpp0311959                1906
## FBpp0070326                  12
## FBpp0087232                 705
## FBpp0305840                 946
## FBpp0081068                1320
## FBpp0290355                   2
## FBpp0082932                   6
## FBpp0311991                 753
## FBpp0311484                  32
## FBpp0290699                   6
## FBpp0303153                   0
## FBpp0073983                1657
## FBpp0070102                 678
## FBpp0075184                   8
## FBpp0082549                 384
## FBpp0075717                   4
## FBpp0084456                   2
## FBpp0084457                   0
## FBpp0082125                   0
## FBpp0289083                 560
## FBpp0308369                 525
## FBpp0077763                 882
## FBpp0289361                  76
## FBpp0082328                  12
## FBpp0312508                 365
## FBpp0309829                1266
## FBpp0077149                 469
## FBpp0070855                   8
## FBpp0083084                  93
## FBpp0099895                 730
## FBpp0085725                1166
## FBpp0301972                  20
## FBpp0082745                   1
## FBpp0110263                   0
## FBpp0084017                 518
## FBpp0081276                   0
## FBpp0072096                 878
## FBpp0305981                   1
## FBpp0087511                1459
## FBpp0304014                  22
## FBpp0076458                 285
## FBpp0311371                1976
## FBpp0298340                   3
## FBpp0307926                   3
## FBpp0082782                   6
## FBpp0305852                3915
## FBpp0079577                 521
## FBpp0080024                 392
## FBpp0072021                 461
## FBpp0099646                 342
## FBpp0310878                 318
## FBpp0297621                1314
## FBpp0304354                 434
## FBpp0290593                  12
## FBpp0302908                   0
## FBpp0304385                 476
## FBpp0297132               17259
## FBpp0082231                 156
## FBpp0309066                 882
## FBpp0289888                 516
## FBpp0081442                 124
## FBpp0081444                 970
## FBpp0311070                   0
## FBpp0075261                 107
## FBpp0306219                1521
## FBpp0309195                 791
## FBpp0306864                   0
## FBpp0302818                 518
## FBpp0288705                 565
## FBpp0304290                 824
## FBpp0081533                 356
## FBpp0312451                 639
## FBpp0088368                3306
## FBpp0303612                   8
## FBpp0088872                  11
## FBpp0312215                  11
## FBpp0311394                8296
## FBpp0084626                 699
## FBpp0071535                 192
## FBpp0290642                 408
## FBpp0080789                   5
## FBpp0293880                   4
## FBpp0086751                  47
## FBpp0075284                 890
## FBpp0310472                 566
## FBpp0074694                   0
## FBpp0086582                2213
## FBpp0083006                   0
## FBpp0308760                   1
## FBpp0083549                   8
## FBpp0303960                   0
## FBpp0305746                  39
## FBpp0073104                  60
## FBpp0078655                9988
## FBpp0303791                 323
## FBpp0309221                 354
## FBpp0099899               11224
## FBpp0086024                 454
## FBpp0292215                3369
## FBpp0080390                 536
## FBpp0304878                 845
## FBpp0074807                   9
## FBpp0304388                3814
## FBpp0071818                 466
## FBpp0309313                   2
## FBpp0309311                   0
## FBpp0112504                1974
## FBpp0301986                   5
## FBpp0309678                1376
## FBpp0311114                 195
## FBpp0312005                 395
## FBpp0311872                  64
## FBpp0289675                 511
## FBpp0300789                   0
## FBpp0078449                2850
## FBpp0070930                  28
## FBpp0077885                 431
## FBpp0077676                 210
## FBpp0079575                 276
## FBpp0074843                1997
## FBpp0305515                2059
## FBpp0311613               38970
## FBpp0084117                 221
## FBpp0081123                 185
## FBpp0311123                   7
## FBpp0110094                   1
## FBpp0082645                1462
## FBpp0312034                 622
## FBpp0082770                 446
## FBpp0072184                 469
## FBpp0077892                   7
## FBpp0085097                   5
## FBpp0087052                 208
## FBpp0288689                   0
## FBpp0288671                 363
## FBpp0078422                 635
## FBpp0079060                   4
## FBpp0293494                  16
## FBpp0309279                   0
## FBpp0306674                   0
## FBpp0304588                  37
## FBpp0312566                  31
## FBpp0308353                   0
## FBpp0078260                   1
## FBpp0087858                   2
## FBpp0306658                  43
## FBpp0099812                1155
## FBpp0072187                 367
## FBpp0297308                   0
## FBpp0297309                   0
## FBpp0297314                   1
## FBpp0083817                   1
## FBpp0307641                2209
## FBpp0088678                 790
## FBpp0088679                 663
## FBpp0099722                 740
## FBpp0071813                 685
## FBpp0088602                  14
## FBpp0312318                  24
## FBpp0307010                  40
## FBpp0083131                 203
## FBpp0304368                 911
## FBpp0303494                3253
## FBpp0292100                   3
## FBpp0305596                  36
## FBpp0290238                   8
## FBpp0309933                1071
## FBpp0084714                 185
## FBpp0088656                1504
## FBpp0088085                  65
## FBpp0080722                 917
## FBpp0075139                 504
## FBpp0303088                  30
## FBpp0303373                   1
## FBpp0306805                   3
## FBpp0076651                   5
## FBpp0099504                  66
## FBpp0292305                   0
## FBpp0305227                 212
## FBpp0305170                  20
## FBpp0078810                1093
## FBpp0077511                1713
## FBpp0075764               21052
## FBpp0305677                1988
## FBpp0306622                1342
## FBpp0076621                   1
## FBpp0075684                 574
## FBpp0083861                4351
## FBpp0087985                 127
## FBpp0290861                   0
## FBpp0074949                1364
## FBpp0070826                   0
## FBpp0075561                1298
## FBpp0311917                2785
## FBpp0084027                  15
## FBpp0086314                 642
## FBpp0303596                  48
## FBpp0071087                  54
## FBpp0290569                 928
## FBpp0076792                 352
## FBpp0306232               79569
## FBpp0082190                  32
## FBpp0311276                5697
## FBpp0087259                   7
## FBpp0303108                1593
## FBpp0306279                5352
## FBpp0292882                2270
## FBpp0304342                 238
## FBpp0293781                 210
## FBpp0073384                   3
## FBpp0303666                4431
## FBpp0072477                 455
## FBpp0070979                  16
## FBpp0308286                   5
## FBpp0297316                   1
## FBpp0308011                 514
## FBpp0311226                2147
## FBpp0085489                2380
## FBpp0306720                   2
## FBpp0304563                3802
## FBpp0075749                  73
## FBpp0311514                 902
## FBpp0079833                 401
## FBpp0070543                 590
## FBpp0303999                 247
## FBpp0112315                 302
## FBpp0112316                 238
## FBpp0086393                 368
## FBpp0309967                 437
## FBpp0303780               13364
## FBpp0087004                 722
## FBpp0303294                 725
## FBpp0082543                 598
## FBpp0080694                 139
## FBpp0292321                 476
## FBpp0310359                   3
## FBpp0310587                   2
## FBpp0073802                   2
## FBpp0076155                2031
## FBpp0082297                 822
## FBpp0291576                 772
## FBpp0088620                 429
## FBpp0090982                 317
## FBpp0292059                  41
## FBpp0304400                2634
## FBpp0075715                1618
## FBpp0304082                 727
## FBpp0070457                1111
## FBpp0301687                  64
## FBpp0309589                  49
## FBpp0075245                   3
## FBpp0308313                  30
## FBpp0305844                 844
## FBpp0086266               24349
## FBpp0080698                  72
## FBpp0086361                   2
## FBpp0075431                1394
## FBpp0290046                3258
## FBpp0080721                 292
## FBpp0310397                 925
## FBpp0112197                 216
## FBpp0077230                 105
## FBpp0304503                 552
## FBpp0074664                1888
## FBpp0311460               29323
## FBpp0112333                 516
## FBpp0078240               21694
## FBpp0084247                 967
## FBpp0310721                 391
## FBpp0072460                 392
## FBpp0078383                1417
## FBpp0305707                  29
## FBpp0081882                3434
## FBpp0074318                 637
## FBpp0086002                6643
## FBpp0081480                 177
## FBpp0081481                 204
## FBpp0304773                   5
## FBpp0088692                 661
## FBpp0081814                 163
## FBpp0271799                  23
## FBpp0289380                   0
## FBpp0303146                2283
## FBpp0310411                5458
## FBpp0302674                 191
## FBpp0079258                 499
## FBpp0290319                   2
## FBpp0290318                   0
## FBpp0072224                 418
## FBpp0077674                 161
## FBpp0081044                   1
## FBpp0077326                 338
## FBpp0083696                  66
## FBpp0081159                1576
## FBpp0305155                   3
## FBpp0309347                   8
## FBpp0303214                2132
## FBpp0311414                  20
## FBpp0306948                  17
## FBpp0087086                7642
## FBpp0073989                2300
## FBpp0080687                5439
## FBpp0110412               18274
## FBpp0074687                 326
## FBpp0073324                  63
## FBpp0304395                   0
## FBpp0081148                 552
## FBpp0293332                  17
## FBpp0293331                   0
## FBpp0073900                 345
## FBpp0304862                  74
## FBpp0081317                1039
## FBpp0072881                 123
## FBpp0110109                  48
## FBpp0304165                1153
## FBpp0112438                1990
## FBpp0305969                   0
## FBpp0071255                  32
## FBpp0311603                2294
## FBpp0110179                 536
## FBpp0290333                1133
## FBpp0070333                 708
## FBpp0311994                 713
## FBpp0297282                3431
## FBpp0081544                3882
## FBpp0085373                1826
## FBpp0301800                 107
## FBpp0074161                2190
## FBpp0311178                1286
## FBpp0088190                1660
## FBpp0305835                2327
## FBpp0082998                 392
## FBpp0083126                 723
## FBpp0070302                 251
## FBpp0083630                1092
## FBpp0310876                 642
## FBpp0074707                 667
## FBpp0076001                 559
## FBpp0085195                6018
## FBpp0292371                 668
## FBpp0112608                 572
## FBpp0081754                1590
## FBpp0073173                1772
## FBpp0307198                 633
## FBpp0087335                 127
## FBpp0083757                  15
## FBpp0304397                   0
## FBpp0087499                 597
## FBpp0078995                   2
## FBpp0070814                 585
## FBpp0306039               24117
## FBpp0311507                  12
## FBpp0110166                 282
## FBpp0075136                 370
## FBpp0081910                 191
## FBpp0075148                 571
## FBpp0080255                 537
## FBpp0071631                2252
## FBpp0085226                  36
## FBpp0070654                  25
## FBpp0307128                 174
## FBpp0077416                 244
## FBpp0308331                   6
## FBpp0079203                 211
## FBpp0088528                  10
## FBpp0289706                1047
## FBpp0302795                1034
## FBpp0086535                1168
## FBpp0306251                 157
## FBpp0082803                 604
## FBpp0303867                 112
## FBpp0099725                1553
## FBpp0099726               10670
## FBpp0075080                 499
## FBpp0070994                1229
## FBpp0074760                  59
## FBpp0288916                 741
## FBpp0288915                 562
## FBpp0087146                 482
## FBpp0305186                 407
## FBpp0077109                   8
## FBpp0308232                2318
## FBpp0084930                 882
## FBpp0288730                 495
## FBpp0080710                 438
## FBpp0082729                 810
## FBpp0082570                 627
## FBpp0079429                 174
## FBpp0302535                  82
## FBpp0302530                1085
## FBpp0072802               31428
## FBpp0304595                 284
## FBpp0090944                 378
## FBpp0290912                1504
## FBpp0311627                 192
## FBpp0082175                  28
## FBpp0082849                 281
## FBpp0309932                 340
## FBpp0079111                1181
## FBpp0082154                  28
## FBpp0088969                   1
## FBpp0309815                2030
## FBpp0086244                 898
## FBpp0086845                 393
## FBpp0076861                 404
## FBpp0084471               16328
## FBpp0084120                1831
## FBpp0087676                 521
## FBpp0078602                3811
## FBpp0082737                 696
## FBpp0306422                1210
## FBpp0311282                2493
## FBpp0302782                  57
## FBpp0305990                2326
## FBpp0077706                   1
## FBpp0083507                 541
## FBpp0072097               20970
## FBpp0312110                5238
## FBpp0079472                2443
## FBpp0082850                 269
## FBpp0071350                 159
## FBpp0111884                  25
## FBpp0081800                 926
## FBpp0293284                   1
## FBpp0312478                3055
## FBpp0305462                2809
## FBpp0293064                 734
## FBpp0083832                   1
## FBpp0311887                2210
## FBpp0074123                 179
## FBpp0081872                  18
## FBpp0293004                 159
## FBpp0083072                  73
## FBpp0076608                 835
## FBpp0081834                2368
## FBpp0310043                  14
## FBpp0076868                 691
## FBpp0291019                2467
## FBpp0071681                 340
## FBpp0311873                 152
## FBpp0081548                 554
## FBpp0072570                1646
## FBpp0072569                1736
## FBpp0113013                 366
## FBpp0309017                 112
## FBpp0289214                  34
## FBpp0312030                 435
## FBpp0309363                 141
## FBpp0075445                  70
## FBpp0307732                1291
## FBpp0311933                  50
## FBpp0304934                 666
## FBpp0071609                 517
## FBpp0077043                 424
## FBpp0081459                 619
## FBpp0309566                1856
## FBpp0308311                 151
## FBpp0087182                 523
## FBpp0309344                  14
## FBpp0080872                1364
## FBpp0308273                 182
## FBpp0306890                 257
## FBpp0302969                  24
## FBpp0070037                4979
## FBpp0308926                1300
## FBpp0309765                 588
## FBpp0079574                 186
## FBpp0304756                  15
## FBpp0311535                 292
## FBpp0307576                  10
## FBpp0091111               11535
## FBpp0071600                1070
## FBpp0311533                1735
## FBpp0111906                1144
## FBpp0311555                 469
## FBpp0071516                1258
## FBpp0309685                1271
## FBpp0110478                 480
## FBpp0088027                 633
## FBpp0079641                 627
## FBpp0086820                 421
## FBpp0075400                  35
## FBpp0311466                   3
## FBpp0312108                  62
## FBpp0084329                1006
## FBpp0085562                 493
## FBpp0075508                1387
## FBpp0307394                   0
## FBpp0303919                 520
## FBpp0087583                7559
## FBpp0311922                 221
## FBpp0300656                 375
## FBpp0083740                1145
## FBpp0072334                 160
## FBpp0297101                   1
## FBpp0307770               19242
## FBpp0088505               15259
## FBpp0302735                 223
## FBpp0307181                1775
## FBpp0303609                   4
## FBpp0305564                  20
## FBpp0080121                 289
## FBpp0070935                 202
## FBpp0307367                 367
## FBpp0082127                  13
## FBpp0085553                   1
## FBpp0303516                2422
## FBpp0303319                  11
## FBpp0100089                  17
## FBpp0077357                 197
## FBpp0077129                1803
## FBpp0305743                 137
## FBpp0073088                 796
## FBpp0085119                 612
## FBpp0073098                1001
## FBpp0087607                1367
## FBpp0081504                 275
## FBpp0073791                2703
## FBpp0099814                 173
## FBpp0080432                  47
## FBpp0073421                  34
## FBpp0305845                 143
## FBpp0304259                   8
## FBpp0081552                2682
## FBpp0311716                  72
## FBpp0302766                1779
## FBpp0087437                1770
## FBpp0083546                 103
## FBpp0306398                  94
## FBpp0077055                   2
## FBpp0293275                1097
## FBpp0311274                 278
## FBpp0310687                4847
## FBpp0080320                2030
## FBpp0078250                   4
## FBpp0293601                  24
## FBpp0309831                  47
## FBpp0075534                 122
## FBpp0309283                 651
## FBpp0307590                 503
## FBpp0074660                 783
## FBpp0074661                 340
## FBpp0074729                  60
## FBpp0288420                   1
## FBpp0070864                 358
## FBpp0082571                6106
## FBpp0084778                 308
## FBpp0087973                 444
## FBpp0087206                  82
## FBpp0075697                 786
## FBpp0305575                   5
## FBpp0083160                 665
## FBpp0083923                 639
## FBpp0307666                1853
## FBpp0079267                 429
## FBpp0112047                 210
## FBpp0309820                 391
## FBpp0297937                   5
## FBpp0309326                   1
## FBpp0297362                2680
## FBpp0303045                   0
## FBpp0073149                 627
## FBpp0306543                  26
## FBpp0084012                 339
## FBpp0088171                   0
## FBpp0073196                 329
## FBpp0077011                 444
## FBpp0311629                 264
## FBpp0081451                 215
## FBpp0081582                1179
## FBpp0077688                   2
## FBpp0086271                   7
## FBpp0073148                 944
## FBpp0086702                  92
## FBpp0309238                 818
## FBpp0292109                  54
## FBpp0075864                  13
## FBpp0075938                 808
## FBpp0081187                  19
## FBpp0303075                 124
## FBpp0311631                  12
## FBpp0309034                1651
## FBpp0292351                 181
## FBpp0084155                 134
## FBpp0113041                5652
## FBpp0310080                 363
## FBpp0085481                  21
## FBpp0304799                  55
## FBpp0082242                  63
## FBpp0072848                2388
## FBpp0291346                2011
## FBpp0087939                2566
## FBpp0088018                1801
## FBpp0083610                   6
## FBpp0074525                 876
## FBpp0305702                9926
## FBpp0297081               11230
## FBpp0073982                 499
## FBpp0303585                 480
## FBpp0087055                 161
## FBpp0290487                 104
## FBpp0072112                1310
## FBpp0288766                1958
## FBpp0309306                 728
## FBpp0293874                 324
## FBpp0073828                1546
## FBpp0271901                 809
## FBpp0309625                 157
## FBpp0289952                1185
## FBpp0083665                 526
## FBpp0080855                 186
## FBpp0079584                 182
## FBpp0078388                 367
## FBpp0306644                3007
## FBpp0297504                  80
## FBpp0072641                 141
## FBpp0083238                1565
## FBpp0073120                 316
## FBpp0311482                 663
## FBpp0082326                  34
## FBpp0309739                7045
## FBpp0078431                2104
## FBpp0084545                 910
## FBpp0311161                 536
## FBpp0079801                 301
## FBpp0312025                   3
## FBpp0305994                 302
## FBpp0301711                 228
## FBpp0307618                 579
## FBpp0304260              247837
## FBpp0082510                 583
## FBpp0311982                5297
## FBpp0293863                 320
## FBpp0303400                 576
## FBpp0079576                1394
## FBpp0305024                1157
## FBpp0288974                1210
## FBpp0309996                1016
## FBpp0312428                1058
## FBpp0084196                  43
## FBpp0079616                 866
## FBpp0311678                  11
## FBpp0086115                 622
## FBpp0072060                1414
## FBpp0309617                  71
## FBpp0083005                 280
## FBpp0308258                 427
## FBpp0301991                 571
## FBpp0081593                 472
## FBpp0086597                 173
## FBpp0309051                  75
## FBpp0289447                  35
## FBpp0307794                  19
## FBpp0296928                 794
## FBpp0082392                 554
## FBpp0291505                  17
## FBpp0086207                2625
## FBpp0076789                 102
## FBpp0078124                 335
## FBpp0310069                1275
## FBpp0290546                1268
## FBpp0310068                  10
## FBpp0081288                 139
## FBpp0308778                 274
## FBpp0073440                 289
## FBpp0078054                 652
## FBpp0072145                 683
## FBpp0289505                  24
## FBpp0082973                1205
## FBpp0289444                2243
## FBpp0293211                 827
## FBpp0099977                 610
## FBpp0291491                2811
## FBpp0089109                2190
## FBpp0311799               13959
## FBpp0312573                 781
## FBpp0088139                  87
## FBpp0304824                 970
## FBpp0289452                  15
## FBpp0086323                 396
## FBpp0310042                   6
## FBpp0111303                 719
## FBpp0087479               22819
## FBpp0072720                 539
## FBpp0073017                2485
## FBpp0291732                  11
## FBpp0288739                 109
## FBpp0087636                 307
## FBpp0089247                   4
## FBpp0307408                  92
## FBpp0297354                  85
## FBpp0083969                 239
## FBpp0311543                 690
## FBpp0086984                 509
## FBpp0081401                9510
## FBpp0083134                 523
## FBpp0083769                 217
## FBpp0302571                 798
## FBpp0076705                  16
## FBpp0076656                  17
## FBpp0088287                1522
## FBpp0310817                  46
## FBpp0303487                 288
## FBpp0074381                1700
## FBpp0309384                5882
## FBpp0304108                1356
## FBpp0288791                 155
## FBpp0288779                3137
## FBpp0291674                  69
## FBpp0309447                 444
## FBpp0087062                1167
## FBpp0099892                6662
## FBpp0306952                   7
## FBpp0080778                 125
## FBpp0304364                  14
## FBpp0111294                1112
## FBpp0075930                5346
## FBpp0080264                  12
## FBpp0081704                1608
## FBpp0309074                1335
## FBpp0086261                2369
## FBpp0088881                6888
## FBpp0300167                1919
## FBpp0300169                1974
## FBpp0290365                   1
## FBpp0070723                 419
## FBpp0083818                 203
## FBpp0082111                1662
## FBpp0079171                2365
## FBpp0308244                  25
## FBpp0311888                 630
## FBpp0304265               15561
## FBpp0305258                1179
## FBpp0298346                 438
## FBpp0077081                1713
## FBpp0077537                  61
## FBpp0086741                  70
## FBpp0290377                 188
## FBpp0290380                 159
## FBpp0089007                  31
## FBpp0075250                 367
## FBpp0085264                 628
## FBpp0309483                 909
## FBpp0074213                1471
## FBpp0085466                 100
## FBpp0099923                3325
## FBpp0085204                 910
## FBpp0070748                 373
## FBpp0081096                 648
## FBpp0074261                1432
## FBpp0305700                 570
## FBpp0311816                4698
## FBpp0307190                   0
## FBpp0085255                1434
## FBpp0112117                4703
## FBpp0086965                 265
## FBpp0291553                  31
## FBpp0309175                 700
## FBpp0083076                 734
## FBpp0072116                 632
## FBpp0075581                3940
## FBpp0307389                 951
## FBpp0088990                1880
## FBpp0304361                 276
## FBpp0073235                 304
## FBpp0084466                 759
## FBpp0077277                1482
## FBpp0312199                3385
## FBpp0086767                 345
## FBpp0303181                   4
## FBpp0298306                  35
## FBpp0088517                6377
## FBpp0289815               23298
## FBpp0311531                1554
## FBpp0076833                5206
## FBpp0311387                5958
## FBpp0293147                  11
## FBpp0293145                   1
## FBpp0293149                  92
## FBpp0304061                 302
## FBpp0075395                 652
## FBpp0075104                3445
## FBpp0309389                  35
## FBpp0309390                  58
## FBpp0086674                 641
## FBpp0072035               11691
## FBpp0082655                1448
## FBpp0309448                4879
## FBpp0087524                 492
## FBpp0304214                1734
## FBpp0070262                1564
## FBpp0297643                 817
## FBpp0291922                 703
## FBpp0291923                 491
## FBpp0291924                 419
## FBpp0305374                 288
## FBpp0305376                 159
## FBpp0086994                3661
## FBpp0087716                2077
## FBpp0308386                1901
## FBpp0271912                1000
## FBpp0074146                  77
## FBpp0110565                  95
## FBpp0293109                  35
## FBpp0304253                1295
## FBpp0306915                3640
## FBpp0112205                  27
## FBpp0073797                   5
## FBpp0075119                1793
## FBpp0079091                1043
## FBpp0100186               48432
## FBpp0072004                 454
## FBpp0078710                 346
## FBpp0311825               13761
## FBpp0073572                 694
## FBpp0304749                 799
## FBpp0298366                 667
## FBpp0076359               21887
## FBpp0075013                 275
## FBpp0075012                3346
## FBpp0293864                 651
## FBpp0310529                 535
## FBpp0311454               10138
## FBpp0110337                1149
## FBpp0290817                1015
## FBpp0085317                 321
## FBpp0308451                2193
## FBpp0303631                 413
## FBpp0079073                1279
## FBpp0071427                8090
## FBpp0288543                 221
## FBpp0304443                 556
## FBpp0075260                 609
## FBpp0071178                 561
## FBpp0087399                  31
## FBpp0290496                  74
## FBpp0305501                1679
## FBpp0082883                 113
## FBpp0308306                 566
## FBpp0072834                  31
## FBpp0072833                  95
## FBpp0312410                 218
## FBpp0310321                1155
## FBpp0070129                 291
## FBpp0099820                 272
## FBpp0085775                1740
## FBpp0309737                 467
## FBpp0307647                2666
## FBpp0088269                1193
## FBpp0312218                 316
## FBpp0311936                 925
## FBpp0303776                 488
## FBpp0077652                2005
## FBpp0070760                3781
## FBpp0304264                 155
## FBpp0308772                2664
## FBpp0291114                  68
## FBpp0291113                  70
## FBpp0291632               51594
## FBpp0075718                 629
## FBpp0081545                1806
## FBpp0071461                7360
## FBpp0071459                1037
## FBpp0292237                 110
## FBpp0073029                 512
## FBpp0073200                  15
## FBpp0303833                 471
```

#gene_mean_expression for sm_male_hdhorn

```r
compute_matrix_row_mean = function(matrix, isLog2 = FALSE){
     apply(matrix,1, compute_mean, isLog2)
}

sm_gene_mean_expression = compute_matrix_row_mean(sm_male_hdhorn)
head(gene_mean_expression) #same gene, accross samples.
```

```
## FBpp0087248 FBpp0293785 FBpp0080383 FBpp0077879 FBpp0311746 FBpp0289081 
##    23.45455  3446.90909    79.54545   139.21818   145.09091  1485.90909
```

#gene_mean_expression for lg_male_hdhorn

```r
compute_matrix_row_mean = function(matrix, isLog2 = FALSE){
  apply(matrix,1, compute_mean, isLog2)
}

lg_gene_mean_expression = compute_matrix_row_mean(lg_male_hdhorn)
head(gene_mean_expression) #same gene, accross samples.
```

```
## FBpp0087248 FBpp0293785 FBpp0080383 FBpp0077879 FBpp0311746 FBpp0289081 
##    23.45455  3446.90909    79.54545   139.21818   145.09091  1485.90909
```

#head and tail of small and large male headhorn to see if it worked

```r
head(lg_gene_mean_expression)
```

```
## FBpp0087248 FBpp0293785 FBpp0080383 FBpp0077879 FBpp0311746 FBpp0289081 
##       21.25     6722.50       66.00        6.25      145.00     1159.50
```

```r
tail(lg_gene_mean_expression)
```

```
## FBpp0071461 FBpp0071459 FBpp0292237 FBpp0073029 FBpp0073200 FBpp0303833 
##     6440.50      939.25       96.00      673.50        7.25      427.25
```

```r
head(sm_gene_mean_expression)
```

```
## FBpp0087248 FBpp0293785 FBpp0080383 FBpp0077879 FBpp0311746 FBpp0289081 
##       37.50     3313.25       71.75       22.50      119.00     1179.00
```

```r
tail(sm_gene_mean_expression)
```

```
## FBpp0071461 FBpp0071459 FBpp0292237 FBpp0073029 FBpp0073200 FBpp0303833 
##     3969.50      710.25       75.50      601.25        4.00      427.50
```

#difference between lg_male_horn and sm_male_horn

```r
diff <- cbind(lg_gene_mean_expression - sm_gene_mean_expression)
tail(diff)
```

```
##                [,1]
## FBpp0071461 2471.00
## FBpp0071459  229.00
## FBpp0292237   20.50
## FBpp0073029   72.25
## FBpp0073200    3.25
## FBpp0303833   -0.25
```

8. Using the basic plot function (although you can use ggplot2 if you prefer), plot the mean expression of each gene on the X axis, and the difference in expression values on the Y axis. Now repeat, but with log2 transformed data. This is the basic idea of a MAplot.

#basic plot of mean expression of each gene on X axis & difference in expression values on Y axis

```r
plot(gene_mean_expression, diff)
```

![](Final_week3_assignment_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
#I donot understand why it shows the positive values
```

#ggplot for gene expression

```r
require(ggplot2)
```

```
## Loading required package: ggplot2
```

```r
ggplot(data = rna_counts, aes(x = gene_mean_expression, y = diff)) +
  geom_point() +
  labs(title = "Basic plot for gene mean expression")
```

![](Final_week3_assignment_files/figure-html/unnamed-chunk-15-1.png)<!-- -->
  
#log2 transformations of mean gene expressions and difference of expreesion values

```r
log2_gene_mean_expression <- log2(gene_mean_expression) #remove -Inf
log2_diff_mean_expression <- log2(diff) #remove -Inf
```

```
## Warning: NaNs produced
```

#basic plot of log2 transformed values

```r
plot(x = log2_gene_mean_expression, y = log2_diff_mean_expression)
```

![](Final_week3_assignment_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

#ggplot for gene expression

```r
ggplot(data = rna_counts, aes(x = log2(gene_mean_expression), y = log2(diff))) +
  geom_point() +
  labs(title = "Basic plot for gene mean expression")
```

```
## Warning in FUN(X[[i]], ...): NaNs produced

## Warning in FUN(X[[i]], ...): NaNs produced
```

```
## Warning: Removed 1025 rows containing missing values (geom_point).
```

![](Final_week3_assignment_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

#Bonus question. What other way might you be able to do these operations (could be a tidyverse way, or a more Rish vectorized way)?

*Answer*: I am not sure if you asked to perfome both the ways. I tried tidyverse way, but I was able to perfom Rish vectorized way more convinently and preferably than the tidyverse way.
