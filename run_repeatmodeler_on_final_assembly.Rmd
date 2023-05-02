---
title: "Run RepeatModeler on Final BFR Assembly"
author: "Jason Toy"
date: "1/28/2023"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Start tmux session

```{bash eval = FALSE}
tmux new -s repeatmodeler
```

## Allocate large node

```{bash eval = FALSE}
salloc --partition=256x44 --nodes=1 --exclusive
ssh $SLURM_NODELIST
```

## Activate RepeatModeler (v1.3)

```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/repeatmodeler
```

## Create working directory

```{bash eval = FALSE}
mkdir /hb/groups/bernardi_lab/jason/BFR_repeatmodeler

cd /hb/groups/bernardi_lab/jason/BFR_repeatmodeler/
```

## Create database from assembly fasta file

```{bash eval = FALSE}
REF=/hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/BFR_ref_nuc_mt_final_v2.fasta

BuildDatabase -name bfr_db $REF
```

## Run RepeatModeler using created database

```{bash eval = FALSE}
RepeatModeler -database bfr_db -pa 42 -LTRStruct > out.log
```

Runtime: 112 hrs!

## Run RepeatMasker using RepeatModeler results

```{bash eval = FALSE}
conda activate /hb/groups/bernardi_lab/programs/repeatmasker

REF=/hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/BFR_ref_nuc_mt_final_v2.fasta
LIB=/hb/groups/bernardi_lab/jason/BFR_repeatmodeler/bfr_db-families.fa

RepeatMasker -s -a -xsmall -pa 42 -lib $LIB -dir RepeatMaskerOutput $REF > repeatmasker_out.log


#-s = Slow search; 0-5% more sensitive, 2-3 times slower than default
#-a = Writes alignments in .align output file
#-xsmall = Returns repetitive regions in lowercase (rest capitals) rather than masked
#-pa = The number of processors to use in parallel
#-lib = custom repeat library from RepeatModeler in fasta format

#-species = Specify the species or clade of the input sequence which is used to determine a taxon-based repeat library for use from Repbase. The species name must be a valid NCBI Taxonomy Database species name and be contained in the RepeatMasker repeat database. This will significantly underestimate the repeat content of the genome so it is recommended to use the custom library created with RepeatModeler instead (using the -lib option).

```

## View table summary

```{bash eval = FALSE}
cd RepeatMaskerOutput

cat BFR_ref_nuc_mt_final_v2.fasta.tbl
```

    ==================================================
    file name: BFR_ref_nuc_mt_final_v2.fasta
    sequences:          1004
    total length:  595951805 bp  (595951805 bp excl N/X-runs)
    GC level:         41.87 %
    bases masked:  125790847 bp ( 21.11 %)
    ==================================================
                   number of      length   percentage
                   elements*    occupied  of sequence
    --------------------------------------------------
    Retroelements        86469     21077852 bp    3.54 %
       SINEs:             9351      1041800 bp    0.17 %
       Penelope           1018       162995 bp    0.03 %
       LINEs:            44752     12474655 bp    2.09 %
        CRE/SLACS            0            0 bp    0.00 %
         L2/CR1/Rex      34775     10327620 bp    1.73 %
         R1/LOA/Jockey    2894       538148 bp    0.09 %
         R2/R4/NeSL        230        66329 bp    0.01 %
         RTE/Bov-B        3248       967986 bp    0.16 %
         L1/CIN4          2112       235645 bp    0.04 %
       LTR elements:     32366      7561397 bp    1.27 %
         BEL/Pao           125        49992 bp    0.01 %
         Ty1/Copia           0            0 bp    0.00 %
         Gypsy/DIRS1      1509       981424 bp    0.16 %
           Retroviral      317       226583 bp    0.04 %

    DNA transposons     122321     28289612 bp    4.75 %
       hobo-Activator    42220      7604541 bp    1.28 %
       Tc1-IS630-Pogo    47897     14656190 bp    2.46 %
       En-Spm                0            0 bp    0.00 %
       MuDR-IS905            0            0 bp    0.00 %
       PiggyBac            113        46327 bp    0.01 %
       Tourist/Harbinger 11650      2690290 bp    0.45 %
       Other (Mirage,     2037       387008 bp    0.06 %
        P-element, Transib)

    Rolling-circles       3541       946114 bp    0.16 %

    Unclassified:       317191     46959181 bp    7.88 %

    Total interspersed repeats:    96326645 bp   16.16 %


    Small RNA:            1930       813972 bp    0.14 %

    Satellites:             37         3076 bp    0.00 %
    Simple repeats:     420466     25546168 bp    4.29 %
    Low complexity:      44448      2291663 bp    0.38 %
    ==================================================

    * most repeats fragmented by insertions or deletions
      have been counted as one element


    RepeatMasker version 4.1.2-p1 , sensitive mode

    run with rmblastn version 2.11.0+
    The query was compared to classified sequences in ".../bfr_db-families.fa"
    FamDB:

## Create contig lengths file for genome for use in next step

```{bash eval = FALSE}
conda activate /hb/groups/bernardi_lab/programs/seqkit/ 

seqkit fx2tab --name --length $REF > BFR_ref_nuc_mt_final_v2_lengths.tsv
```

Do this next step on the head node. Otherwise it will have issues finding perl modules:

## Add pm library location from conda build to PERLLIB environment variable so perl knows where to find them

```{bash eval = FALSE}
export PERLLIB="$PERLLIB:/hb/groups/bernardi_lab/programs/repeatmasker/share/RepeatMasker"
```

## Summarize the RepeatMasker results with buildSummary.pl utility

```{bash eval = FALSE}
conda activate /hb/groups/bernardi_lab/programs/repeatmasker

buildSummary.pl -genome BFR_ref_nuc_mt_final_v2_lengths.tsv -useAbsoluteGenomeSize BFR_ref_nuc_mt_final_v2.fasta.out > repeatmasker_summary
 
head -74 repeatmasker_summary > repeatmasker_summary_repeat_classes.txt
```

    Repeat Classes
    ==============
    Total Sequences: 1004
    Total Length: 595951805 bp
    Class                  Count        bpMasked    %masked
    =====                  =====        ========     =======
    DNA                    8994         1373980      0.23%
        CMC-Chapaev-3      430          89983        0.02%
        CMC-EnSpm          796          132577       0.02%
        Crypton-A          1678         260468       0.04%
        Crypton-V          966          157472       0.03%
        IS3EU              2754         454237       0.08%
        Kolobok-T2         174          33988        0.01%
        MULE-MuDR          280          77249        0.01%
        Maverick           724          147974       0.02%
        Merlin             141          22930        0.00%
        P                  2038         387008       0.06%
        PIF-Harbinger      11688        2690290      0.45%
        PiggyBac           113          46327        0.01%
        TcMar              104          3507         0.00%
        TcMar-Fot1         51           12358        0.00%
        TcMar-ISRm11       98           32796        0.01%
        TcMar-Mariner      11042        2489547      0.42%
        TcMar-Tc1          35600        11695020     1.96%
        TcMar-Tc2          76           27082        0.00%
        TcMar-Tigger       1311         395880       0.07%
        Zator              647          100909       0.02%
        Zisupton           836          53489        0.01%
        hAT                1361         131190       0.02%
        hAT-Ac             28328        5328612      0.89%
        hAT-Blackjack      277          53210        0.01%
        hAT-Charlie        9948         1638080      0.27%
        hAT-Tip100         1124         208177       0.03%
        hAT-hAT5           977          191814       0.03%
        hAT-hAT6           340          53458        0.01%
    LINE                   93           53178        0.01%
        Dong-R4            62           22588        0.00%
        I                  2894         538148       0.09%
        I-Jockey           65           40087        0.01%
        L1                 2112         235645       0.04%
        L1-Tx1             255          65980        0.01%
        L2                 23197        6801870      1.14%
        Penelope           1018         162995       0.03%
        Proto2             62           16687        0.00%
        R2-Hero            168          43741        0.01%
        RTE-BovB           3021         896697       0.15%
        RTE-RTE            23           23971        0.00%
        RTE-X              204          47318        0.01%
        Rex-Babar          11578        3525750      0.59%
    LTR                    56           20549        0.00%
        ERV1               317          226583       0.04%
        Gypsy              1509         981424       0.16%
        Ngaro              11477        2377153      0.40%
        Pao                125          49992        0.01%
        Unknown            18882        3905696      0.66%
    RC                     --           --           --
        Helitron           3597         946114       0.16%
    Retroposon             79           12155        0.00%
    SINE                   --           --           --
        5S-Deu-L2          195          27250        0.00%
        tRNA               1134         136791       0.02%
        tRNA-Core-RTE      698          71738        0.01%
        tRNA-V-RTE         7326         806021       0.14%
    Unknown                318293       46947026     7.88%
                          ---------------------------------
        total interspersed 531336       97272759     16.32%

    Low_complexity         44448        2291663      0.38%
    Satellite              37           3076         0.00%
    Simple_repeat          420466       25546168     4.29%
    rRNA                   659          658915       0.11%
    tRNA                   142          18266        0.00%
    ---------------------------------------------------------
    Total                  997088       125790847    21.11%

```{bash eval = FALSE}
tail -1005 repeatmasker_summary > repeatmasker_summary_by_contig.txt
```

    By Sequence
    ===========
    Seq            Count      bpMasked
    =====          =====      ========
          Contig_324   739        104599
          Contig_236   1286       168751
          Contig_263   883        79462
          Contig_793   7          2583
           Contig_90   3224       491415
          Contig_600   21         2324
          Contig_105   2887       480319
          Contig_532   101        30068
          Contig_276   962        173851
          Contig_494   192        47598
          Contig_191   1420       235293
          Contig_609   86         18595
          Contig_635   54         15364
           Contig_46   6308       906948
          Contig_393   393        54373
          Contig_278   900        152677
          Contig_835   9          1540
          Contig_585   101        16529
           Contig_34   5746       746307
          Contig_641   63         20558
          Contig_484   238        33900
          Contig_382   480        77993
          Contig_309   582        66453
          Contig_724   16         16230
          Contig_482   253        36634
          Contig_157   1009       119577
           Contig_89   3805       344421
          Contig_965   12         358
          Contig_568   188        16124
          Contig_118   2237       339786
          Contig_775   4          2267
          Contig_862   2          6367
          Contig_533   147        28394
          Contig_314   723        130907
          Contig_172   1545       248169
            Contig_5   16295      1653641
          Contig_254   710        64119
          Contig_800   6          538
          Contig_344   533        75727
          Contig_709   32         9617
          Contig_162   1636       167586
          Contig_487   267        45339
          Contig_949   5          1349
          Contig_879   3          5970
          Contig_457   261        52830
          Contig_411   366        58013
          Contig_838   1          7305
          Contig_787   1          9824
          Contig_243   812        126003
          Contig_476   263        64986
           Contig_61   3450       392348
          Contig_944   1          4134
          Contig_110   2013       253950
          Contig_881   1          5965
           Contig_77   2262       235012
          Contig_632   57         12862
          Contig_594   116        14773
          Contig_421   408        77268
          Contig_994   2          1690
          Contig_260   1016       80865
           Contig_84   3554       377162
          Contig_811   4          227
          Contig_202   1340       310808
          Contig_553   180        19987
          Contig_497   208        22500
          Contig_227   1061       116403
          Contig_703   28         5473
          Contig_672   23         19170
          Contig_515   257        22435
          Contig_177   1738       207828
          Contig_737   13         1129
          Contig_526   99         29818
          Contig_778   51         5198
          Contig_740   1          13685
          Contig_545   154        17748
          Contig_624   80         7805
          Contig_242   1190       99979
          Contig_859   8          4866
          Contig_441   280        34619
          Contig_943   19         635
           Contig_22   7719       1011995
          Contig_231   1485       149889
          Contig_431   419        28877
          Contig_853   5          4923
          Contig_849   9          4806
          Contig_567   85         17235
          Contig_828   2          7759
          Contig_302   655        72409
          Contig_690   22         5772
          Contig_187   1084       96786
          Contig_293   859        141152
          Contig_685   29         22835
          Contig_126   2639       319468
          Contig_925   8          2983
          Contig_902   5          5085
          Contig_780   2          10534
          Contig_616   70         17500
          Contig_563   92         17745
          Contig_288   426        52275
          Contig_153   1484       184543
          Contig_462   236        45297
          Contig_974   5          206
          Contig_970   39         3210
          Contig_757   19         5944
          Contig_506   204        32924
          Contig_387   498        126726
          Contig_354   592        49268
          Contig_539   189        33434
          Contig_138   2005       456653
          Contig_755   1          12449
          Contig_863   19         4426
            Contig_8   11995      1743689
          Contig_343   453        78024
          Contig_213   1233       148332
          Contig_818   7          2941
          Contig_732   3          14699
           Contig_70   3705       347246
          Contig_410   257        26565
          Contig_758   12         12144
          Contig_334   527        56578
          Contig_645   47         16122
          Contig_483   210        15731
          Contig_654   66         9983
          Contig_287   819        179283
          Contig_380   458        67752
          Contig_781   27         5559
          Contig_964   1          3396
          Contig_633   129        9791
          Contig_975   13         322
          Contig_798   6          8580
          Contig_150   2450       271066
          Contig_280   850        206923
          Contig_392   359        94264
          Contig_142   1367       186685
          Contig_388   664        56303
          Contig_590   116        26666
          Contig_209   1108       202533
          Contig_882   14         1821
          Contig_437   321        50242
          Contig_586   129        20685
          Contig_127   2374       276343
          Contig_830   6          6210
          Contig_244   916        99979
          Contig_192   1437       270901
          Contig_556   108        15905
          Contig_971   15         434
          Contig_948   4          299
          Contig_936   1          4472
          Contig_660   71         10136
          Contig_480   222        43977
          Contig_702   36         6898
          Contig_486   281        46814
          Contig_836   3          7381
          Contig_189   1210       97312
          Contig_988   2          2122
          Contig_181   1331       165830
          Contig_543   162        32707
          Contig_361   344        38105
          Contig_713   40         9062
           Contig_33   7449       692837
          Contig_640   62         17860
          Contig_812   1          8520
          Contig_444   297        59212
           Contig_60   4439       761041
          Contig_449   317        47104
          Contig_160   1644       271198
          Contig_870   1          6111
          Contig_719   24         4671
          Contig_572   133        30111
          Contig_691   29         5949
          Contig_522   161        28495
         Contig_1000   1          1296
           Contig_21   10023      889820
          Contig_942   3          4205
          Contig_841   58         7002
          Contig_339   815        65456
           Contig_97   1893       203471
          Contig_825   11         4644
          Contig_141   2778       250077
          Contig_329   494        106068
          Contig_851   3          6452
          Contig_104   3970       295794
          Contig_802   17         8954
           Contig_54   3626       347657
          Contig_892   22         3352
          Contig_889   13         4147
          Contig_458   219        68590
          Contig_848   4          3361
          Contig_525   193        30234
          Contig_459   348        58992
          Contig_326   496        54635
          Contig_473   232        69284
          Contig_178   1281       167111
          Contig_601   122        16208
          Contig_980   3          2543
          Contig_241   1507       168147
          Contig_188   1327       192975
          Contig_359   322        72571
           Contig_67   3073       294408
          Contig_323   1155       86946
          Contig_128   1872       224509
          Contig_436   285        45431
          Contig_540   184        29768
           Contig_78   3720       637346
          Contig_115   2566       230283
          Contig_799   14         9028
          Contig_915   2          3121
          Contig_924   1          4881
          Contig_469   224        38432
          Contig_510   166        17324
          Contig_503   103        12359
          Contig_564   87         32518
           Contig_47   3906       370992
          Contig_374   594        48930
          Contig_852   4          1057
          Contig_463   264        48438
          Contig_819   11         7903
          Contig_461   273        46404
          Contig_277   1222       92583
          Contig_139   3184       261008
          Contig_794   1          9627
          Contig_987   1          2206
          Contig_810   1          8576
          Contig_518   108        35949
          Contig_435   426        41047
          Contig_303   647        89423
          Contig_901   6          3664
          Contig_638   61         8836
          Contig_659   33         29016
          Contig_817   9          7831
          Contig_679   49         8662
          Contig_941   5          4248
          Contig_756   16         12342
          Contig_887   10         3584
          Contig_864   5          614
          Contig_768   27         3319
          Contig_984   10         281
          Contig_689   23         22364
          Contig_565   116        27694
          Contig_628   79         8825
          Contig_164   2325       232457
          Contig_653   60         4134
          Contig_786   11         10026
          Contig_614   72         22003
          Contig_947   1          4079
           Contig_44   5614       797462
          Contig_701   35         8249
          Contig_643   63         7977
          Contig_766   21         11318
          Contig_900   1          5278
          Contig_854   2          6614
          Contig_677   25         24285
          Contig_517   190        47947
          Contig_375   384        76894
          Contig_345   551        43165
          Contig_855   3          1483
          Contig_630   66         13507
           Contig_56   6256       587022
          Contig_499   123        33590
          Contig_143   2713       232718
          Contig_621   53         15968
          Contig_464   229        46426
          Contig_245   842        147244
          Contig_582   104        17503
          Contig_413   356        78569
          Contig_927   3          4814
          Contig_320   813        78723
          Contig_932   14         3313
          Contig_591   99         15536
          Contig_331   800        64963
          Contig_294   760        167660
          Contig_972   11         328
           Contig_32   8027       1000286
          Contig_680   72         8571
          Contig_400   473        110139
          Contig_225   777        87983
          Contig_729   29         7791
           Contig_18   8027       934993
          Contig_403   395        70011
          Contig_578   100        9801
          Contig_583   119        26387
          Contig_665   37         4502
          Contig_443   191        62415
          Contig_829   13         1864
          Contig_951   2          3889
          Contig_865   3          6269
          Contig_224   1137       194050
          Contig_238   700        68056
          Contig_826   11         7897
          Contig_504   168        50247
          Contig_417   609        59096
          Contig_708   28         7632
          Contig_352   579        116995
          Contig_132   1923       223605
          Contig_489   345        30299
          Contig_742   24         8158
          Contig_573   114        19270
          Contig_842   11         3650
           Contig_39   4403       505045
          Contig_140   1537       149069
          Contig_662   54         10136
          Contig_552   109        8791
          Contig_743   24         3344
          Contig_626   41         13852
          Contig_788   2          9759
          Contig_325   544        61493
          Contig_555   85         26091
          Contig_496   192        51979
           Contig_15   7898       989107
               mtDNA   3          329
          Contig_639   62         15794
          Contig_893   2          5616
           Contig_83   2057       272683
          Contig_257   892        76425
          Contig_906   11         2645
          Contig_721   31         3725
          Contig_559   115        20554
          Contig_434   289        55413
          Contig_315   668        97075
          Contig_699   36         8345
          Contig_921   9          316
          Contig_529   128        38391
          Contig_608   39         44682
          Contig_281   931        104286
          Contig_297   701        138655
          Contig_237   1091       172951
          Contig_692   28         4322
          Contig_557   101        18086
          Contig_107   2472       375538
          Contig_967   7          384
          Contig_674   57         11075
          Contig_535   130        27779
          Contig_363   457        68054
          Contig_523   160        17701
          Contig_791   117        9638
          Contig_250   833        146930
          Contig_264   947        178785
          Contig_197   1239       275893
          Contig_541   113        41210
          Contig_715   2          18063
          Contig_179   2452       161817
          Contig_809   9          8578
           Contig_74   2226       272525
          Contig_460   277        59699
          Contig_271   1247       115351
          Contig_923   2          4882
          Contig_647   47         9002
          Contig_212   1082       220156
          Contig_837   5          1007
          Contig_735   34         8917
            Contig_7   15111      2048153
          Contig_648   80         16148
           Contig_23   8241       1036530
          Contig_381   436        74621
           Contig_62   3604       396616
          Contig_117   1489       185245
          Contig_581   68         13251
          Contig_613   51         21100
          Contig_279   818        161235
          Contig_577   101        10815
          Contig_501   230        42003
          Contig_475   135        14605
          Contig_465   299        26969
          Contig_678   39         8210
          Contig_412   402        34192
          Contig_346   642        74699
          Contig_866   7          5976
          Contig_379   474        73230
          Contig_991   1          2044
          Contig_723   2          16628
          Contig_642   48         7019
          Contig_332   522        55582
          Contig_198   1599       254004
           Contig_75   4595       695625
          Contig_966   1          3369
          Contig_446   289        36938
          Contig_376   584        88991
          Contig_888   4          5789
          Contig_301   477        75328
          Contig_765   17         5724
          Contig_615   118        20504
          Contig_833   1          7397
          Contig_300   520        82392
          Contig_681   35         7587
          Contig_161   1454       315192
          Contig_240   670        74158
          Contig_401   394        99251
          Contig_566   70         18232
          Contig_485   249        39346
          Contig_144   1058       115377
          Contig_599   81         12811
          Contig_330   1015       83397
          Contig_592   48         3705
          Contig_389   469        80325
          Contig_112   3442       238019
          Contig_738   5          604
           Contig_95   4124       348207
          Contig_512   222        37001
          Contig_131   2161       305437
          Contig_688   24         10092
          Contig_171   1544       301244
          Contig_938   11         389
          Contig_571   84         14262
           Contig_76   3620       374526
          Contig_569   74         13821
          Contig_950   3          3934
          Contig_750   23         10422
          Contig_985   10         530
          Contig_373   497        68343
          Contig_295   652        71534
          Contig_574   76         10706
          Contig_360   593        83862
          Contig_584   37         8176
          Contig_316   916        90078
          Contig_847   3          6790
          Contig_416   327        46757
          Contig_262   768        146783
          Contig_629   87         13537
          Contig_940   10         464
           Contig_55   3679       418089
          Contig_394   544        45855
          Contig_813   9          2305
          Contig_922   9          4902
          Contig_125   1812       216335
          Contig_730   16         10362
          Contig_429   355        48305
          Contig_186   1036       105332
          Contig_488   377        42038
          Contig_880   2          167
          Contig_707   23         3141
           Contig_45   5207       744863
          Contig_408   353        38152
           Contig_80   3693       279611
          Contig_353   500        107656
          Contig_760   15         11975
          Contig_675   28         7682
          Contig_230   1080       107576
          Contig_989   1          2093
          Contig_714   38         8860
          Contig_627   62         23472
          Contig_546   40         6880
           Contig_98   2696       272323
          Contig_137   2756       275477
          Contig_256   1671       108989
            Contig_1   19568      2399977
           Contig_14   8819       1136866
          Contig_973   1          3059
          Contig_337   388        43781
          Contig_538   158        30922
          Contig_509   170        23608
          Contig_524   218        34485
          Contig_355   699        59958
          Contig_920   10         3241
          Contig_534   122        15882
          Contig_159   2558       234031
          Contig_452   309        55891
          Contig_384   356        50519
           Contig_31   5047       546694
          Contig_698   29         7142
          Contig_816   9          3233
          Contig_106   2749       425558
          Contig_916   10         303
          Contig_741   27         4748
          Contig_700   30         4074
          Contig_203   1079       101708
          Contig_792   11         9542
          Contig_551   123        22978
          Contig_282   476        35558
          Contig_904   1          452
           Contig_41   4578       518612
          Contig_939   57         4269
           Contig_57   3790       406607
          Contig_222   923        134713
          Contig_154   1475       200677
          Contig_646   63         8929
          Contig_364   473        109574
          Contig_894   2          5551
          Contig_169   1114       152129
          Contig_211   1117       108601
          Contig_265   840        122785
          Contig_542   160        26361
          Contig_663   56         10406
           Contig_86   3036       473250
            Contig_9   11254      1553563
           Contig_85   2933       495336
          Contig_885   12         5767
          Contig_997   2          39
          Contig_667   36         3334
          Contig_290   825        179438
          Contig_834   3          7543
          Contig_784   13         10332
           Contig_13   9988       988541
          Contig_471   312        28746
          Contig_693   54         7230
          Contig_575   123        34250
          Contig_350   466        48147
          Contig_514   131        29436
          Contig_929   3          3256
          Contig_895   3          5320
          Contig_771   9          11038
          Contig_710   14         18548
          Contig_598   83         21157
          Contig_219   1229       176621
          Contig_502   172        37577
          Contig_983   1          2458
          Contig_824   14         1556
          Contig_109   2201       265589
           Contig_30   5930       709272
          Contig_306   988        83257
          Contig_347   466        61481
          Contig_367   386        55138
          Contig_736   15         5037
          Contig_232   1173       189356
          Contig_867   14         3029
          Contig_845   3          6838
           Contig_72   2914       366361
          Contig_839   13         1601
          Contig_589   56         6085
          Contig_251   820        99442
          Contig_840   10         2712
          Contig_223   1557       127848
          Contig_896   5          5413
          Contig_299   767        115306
          Contig_918   9          2222
          Contig_580   113        25962
          Contig_365   742        62394
           Contig_38   5732       538723
          Contig_953   1          3892
          Contig_433   354        72523
          Contig_466   274        42973
          Contig_875   3          5987
          Contig_619   72         20792
          Contig_317   778        110013
          Contig_910   4          475
           Contig_37   9255       825670
          Contig_340   392        65242
          Contig_969   1          3249
          Contig_790   28         1511
           Contig_52   4648       470815
          Contig_745   15         13049
           Contig_19   8093       1090792
           Contig_25   8223       1217568
          Contig_907   2          391
          Contig_422   367        46459
          Contig_163   2266       184146
          Contig_173   1369       142177
          Contig_519   197        39395
          Contig_990   1          2047
          Contig_152   1670       355653
          Contig_873   9          6021
          Contig_717   45         6522
          Contig_763   7          11726
          Contig_284   754        71339
          Contig_156   1431       200729
           Contig_65   4321       406755
           Contig_93   2927       246801
          Contig_687   17         5133
          Contig_377   290        40372
          Contig_631   36         8097
          Contig_194   1530       238959
          Contig_383   499        78687
          Contig_428   160        18234
          Contig_246   956        164849
          Contig_114   2588       322174
          Contig_606   67         15001
          Contig_166   1672       181279
          Contig_205   1948       159823
          Contig_100   3133       537416
           Contig_16   9389       1404192
          Contig_548   166        35133
          Contig_611   86         15384
          Contig_134   1355       139394
          Contig_749   14         3275
          Contig_283   738        91212
          Contig_807   17         3654
          Contig_370   565        50274
          Contig_806   4          2554
          Contig_725   24         3224
          Contig_579   88         25102
          Contig_477   223        56193
           Contig_88   3165       418346
          Contig_652   36         13099
          Contig_356   643        93401
          Contig_814   14         8324
          Contig_908   2          5198
          Contig_328   608        133243
          Contig_409   368        67017
          Contig_490   246        55152
           Contig_58   3817       405264
          Contig_402   344        72798
          Contig_120   1774       196795
          Contig_905   15         3200
          Contig_195   1019       228111
          Contig_993   1          1638
          Contig_405   370        97954
          Contig_122   1250       172092
          Contig_720   31         5695
          Contig_396   610        53182
          Contig_210   1575       159172
          Contig_101   2801       450992
          Contig_657   43         22270
          Contig_595   29         2763
          Contig_623   69         12221
          Contig_676   63         10891
          Contig_978   2          2605
          Contig_273   1114       117360
          Contig_831   10         7705
          Contig_774   29         5953
          Contig_448   328        56170
          Contig_498   160        46045
          Contig_451   302        65899
          Contig_706   1          19066
          Contig_759   28         6576
          Contig_183   884        97390
          Contig_804   18         8820
          Contig_311   754        92217
          Contig_310   567        64407
          Contig_259   692        87705
          Contig_658   53         17248
          Contig_129   2613       247814
          Contig_378   508        57569
          Contig_185   1295       248358
          Contig_291   882        106553
          Contig_772   1          10997
          Contig_351   691        70704
          Contig_108   2377       228899
          Contig_767   9          11240
          Contig_823   3          8003
          Contig_751   9          2502
           Contig_73   4929       450215
          Contig_926   14         473
          Contig_124   1243       160367
          Contig_992   2          1923
          Contig_946   29         4093
          Contig_671   28         2801
          Contig_220   597        76287
          Contig_716   32         5618
          Contig_664   66         13074
          Contig_560   63         16502
           Contig_59   4784       780618
          Contig_318   600        64373
          Contig_348   519        138066
          Contig_415   255        41256
          Contig_423   302        94263
          Contig_217   1088       102064
          Contig_868   1          6133
           Contig_82   4282       444385
            Contig_2   18698      2071547
          Contig_233   582        81266
          Contig_327   319        43948
          Contig_180   2221       201505
          Contig_145   1524       189658
          Contig_537   136        12382
           Contig_40   5881       902422
          Contig_218   1202       240634
          Contig_620   78         20604
          Contig_883   6          3031
          Contig_846   3          6745
          Contig_453   413        43606
          Contig_933   3          1370
          Contig_558   164        13839
          Contig_911   2          5185
          Contig_668   43         6638
          Contig_204   853        78144
          Contig_312   666        124317
          Contig_636   57         11795
          Contig_305   566        72436
          Contig_174   1733       214195
          Contig_815   6          1169
          Contig_722   26         3762
          Contig_527   179        38962
          Contig_857   8          6466
          Contig_151   2039       167409
          Contig_432   329        89945
          Contig_917   4          3951
          Contig_493   183        46820
          Contig_133   1909       265177
           Contig_12   10313      1106429
          Contig_968   6          226
          Contig_253   975        78632
           Contig_53   4959       517024
          Contig_956   2          3768
          Contig_785   9          10146
          Contig_744   30         4066
          Contig_808   3          7948
          Contig_612   74         18599
          Contig_979   1          2597
           Contig_91   4374       459901
          Contig_119   2098       213043
          Contig_335   472        84321
          Contig_491   171        22731
          Contig_544   85         11551
          Contig_856   7          800
          Contig_165   1245       133229
          Contig_130   1751       214917
          Contig_196   774        70453
          Contig_123   3250       357209
          Contig_897   6          620
          Contig_511   228        20493
           Contig_63   3580       516649
          Contig_239   1588       126785
          Contig_651   33         9356
          Contig_338   611        136425
          Contig_547   138        21368
           Contig_64   2731       315836
          Contig_267   526        68994
          Contig_876   10         3814
          Contig_478   208        41892
           Contig_87   3418       546049
           Contig_17   7971       942908
          Contig_576   95         32699
          Contig_322   533        60224
           Contig_50   4388       444440
          Contig_597   86         25161
          Contig_247   1499       113458
          Contig_111   3442       322204
          Contig_146   1733       178565
          Contig_366   449        111214
          Contig_779   17         6960
          Contig_445   322        70316
          Contig_170   1866       303264
          Contig_996   1          1649
          Contig_607   38         30422
          Contig_507   209        39172
          Contig_952   6          3791
          Contig_705   65         4944
          Contig_570   76         27809
          Contig_440   244        20560
          Contig_308   591        50568
          Contig_694   35         5918
          Contig_649   49         15797
          Contig_270   610        57134
          Contig_258   823        152523
          Contig_414   321        59899
           Contig_96   2627       233170
          Contig_285   572        78366
           Contig_49   4269       633723
          Contig_184   1777       217605
          Contig_513   187        36675
          Contig_726   18         15807
          Contig_296   697        114882
          Contig_395   530        67975
          Contig_886   25         3473
          Contig_682   43         9689
          Contig_957   2          3724
          Contig_587   101        20938
          Contig_342   645        122868
          Contig_492   136        37031
          Contig_871   10         5984
          Contig_981   12         443
          Contig_661   63         9932
          Contig_747   27         5054
           Contig_28   8033       1008866
          Contig_762   13         6700
          Contig_261   711        95528
          Contig_248   925        146112
          Contig_644   64         15272
          Contig_843   11         6720
          Contig_274   840        111659
            Contig_3   19832      2555994
          Contig_175   1636       207285
          Contig_858   88         6480
          Contig_207   1473       143454
          Contig_216   1198       184541
          Contig_695   37         6412
          Contig_148   1666       146963
          Contig_121   1936       359197
          Contig_385   410        29052
          Contig_427   255        42206
          Contig_773   13         10830
          Contig_752   10         12651
          Contig_596   108        19572
          Contig_479   189        30587
          Contig_914   1          5136
          Contig_362   349        43616
          Contig_505   190        36123
          Contig_931   20         674
            Contig_6   12940      1486447
           Contig_92   3481       406152
          Contig_999   6          173
          Contig_937   2          4461
           Contig_79   3875       461680
          Contig_822   79         8036
          Contig_424   373        43298
          Contig_963   46         3500
          Contig_455   261        53059
          Contig_731   26         7051
          Contig_228   1039       92016
          Contig_454   246        24127
          Contig_874   10         1234
          Contig_221   1177       256843
          Contig_269   926        155138
           Contig_11   11302      1165858
          Contig_252   788        123060
          Contig_234   1038       165280
          Contig_149   2056       371405
          Contig_349   559        92955
          Contig_960   54         3644
          Contig_208   1259       223814
          Contig_686   45         8867
          Contig_438   195        25897
          Contig_909   1          5188
          Contig_625   68         20244
           Contig_24   6953       894817
          Contig_406   252        36240
          Contig_796   17         4347
          Contig_655   60         7192
          Contig_201   760        73443
          Contig_508   186        44134
          Contig_884   3          3394
           Contig_71   4107       693890
           Contig_43   4371       427740
          Contig_789   16         3042
          Contig_561   143        26019
          Contig_734   1          14637
          Contig_268   989        149890
           Contig_81   3291       321037
          Contig_617   33         5578
           Contig_36   6081       728969
          Contig_439   310        68946
          Contig_976   1          2720
          Contig_728   19         15260
           Contig_51   5535       483568
          Contig_521   126        39312
          Contig_168   1954       186763
          Contig_467   203        59366
           Contig_48   4248       435327
          Contig_214   1410       210571
          Contig_890   1          5787
          Contig_358   365        59378
          Contig_669   69         9696
          Contig_313   618        122845
           Contig_27   6339       779238
           Contig_68   2086       231453
          Contig_801   1          9071
          Contig_199   1273       210451
          Contig_604   71         12638
          Contig_912   12         5152
          Contig_712   26         3047
          Contig_683   20         3151
          Contig_531   181        21103
          Contig_945   8          209
          Contig_116   2291       423565
          Contig_995   1          1651
          Contig_637   71         15447
          Contig_136   1129       109768
          Contig_418   474        41790
          Contig_474   252        49915
          Contig_670   44         15593
          Contig_704   35         11080
          Contig_753   171        12662
          Contig_304   920        96275
          Contig_986   34         2358
           Contig_20   7695       868718
          Contig_470   268        27295
          Contig_955   4          3784
          Contig_371   335        102827
          Contig_158   1987       146204
          Contig_399   397        59819
          Contig_958   1          3756
          Contig_666   29         26530
          Contig_602   68         21114
          Contig_850   20         1943
          Contig_319   846        70536
          Contig_398   342        70114
          Contig_472   202        35734
          Contig_368   231        24226
          Contig_286   781        67453
          Contig_877   10         2718
          Contig_782   38         4674
          Contig_103   1738       138281
          Contig_961   13         377
          Contig_903   3          5081
          Contig_872   13         5251
          Contig_229   783        102286
          Contig_550   171        15192
          Contig_593   74         18066
          Contig_588   66         17885
          Contig_783   12         10254
          Contig_962   4          3542
          Contig_820   18         3536
          Contig_391   388        73622
          Contig_226   1112       117834
          Contig_776   33         4220
          Contig_562   76         20248
          Contig_407   545        71448
          Contig_390   429        75509
          Contig_861   13         5154
          Contig_696   3          20294
          Contig_167   2145       203682
          Contig_684   46         14125
          Contig_844   9          6917
           Contig_42   3631       422916
          Contig_536   231        23848
          Contig_919   7          799
          Contig_733   35         3585
          Contig_307   620        70691
          Contig_397   502        47693
          Contig_430   324        85283
          Contig_266   917        111455
          Contig_190   973        109588
          Contig_795   10         1749
          Contig_770   9          11154
          Contig_200   1662       153801
            Contig_4   14989      1968749
          Contig_718   23         8717
          Contig_275   903        123294
          Contig_739   13         2126
          Contig_235   680        75577
          Contig_450   244        67921
          Contig_426   258        68066
           Contig_26   6209       494007
          Contig_147   1838       330273
          Contig_298   742        111248
          Contig_336   656        114613
           Contig_99   2676       349618
          Contig_935   8          4333
          Contig_697   39         6623
          Contig_673   1          24841
          Contig_748   23         3665
          Contig_369   587        58791
          Contig_869   1          6070
          Contig_255   1046       237235
           Contig_10   10793      1271180
          Contig_656   93         13250
          Contig_176   1453       148333
           Contig_35   9329       888991
          Contig_495   127        23955
          Contig_321   545        73487
          Contig_898   1          5397
          Contig_761   3          116
           Contig_66   3879       387081
          Contig_520   210        37875
          Contig_891   3          4370
          Contig_982   2          280
          Contig_420   318        49990
          Contig_456   165        16951
          Contig_386   432        64671
          Contig_860   21         727
           Contig_29   6595       860237
          Contig_341   472        54260
          Contig_899   2          5326
          Contig_805   9          6041
          Contig_193   741        62706
          Contig_610   52         13341
          Contig_549   93         20807
          Contig_998   2          1420
          Contig_764   1          11635
          Contig_500   192        25949
          Contig_419   267        26318
          Contig_481   182        18445
          Contig_113   3088       270132
          Contig_930   9          292
          Contig_746   41         5796
          Contig_155   925        108843
          Contig_959   1          3700
          Contig_634   73         13035
          Contig_215   891        119496
          Contig_333   506        93034
          Contig_530   148        40358
          Contig_289   450        63420
          Contig_769   23         4183
          Contig_650   59         5786
          Contig_102   2264       255029
          Contig_272   680        75309
          Contig_754   15         2273
          Contig_372   489        98360
          Contig_928   2          4819
          Contig_913   1          32
          Contig_727   1          15825
          Contig_618   55         23320
           Contig_94   2156       223491
          Contig_821   13         1247
          Contig_977   17         662
          Contig_954   1          3874
          Contig_934   5          4554
          Contig_827   3          7807
          Contig_442   229        52175
          Contig_468   205        46089
          Contig_447   282        59690
          Contig_425   372        59179
          Contig_528   151        31510
          Contig_206   1693       160937
          Contig_135   2308       227826
          Contig_605   69         13431
          Contig_797   22         5108
          Contig_832   30         2025
          Contig_249   1052       192986
          Contig_622   78         9177
          Contig_516   197        34462
          Contig_292   791        150740
          Contig_803   10         8925
          Contig_777   11         10485
           Contig_69   4452       612231
          Contig_711   30         3795
          Contig_357   478        90378
          Contig_404   354        62526
          Contig_554   106        17839
          Contig_878   7          5982
          Contig_603   75         15766
          Contig_182   1370       171338

## Calculate the Kimura divergence values with calcDivergenceFromAlign.pl utility

```{bash eval = FALSE}
calcDivergenceFromAlign.pl -a divergence_align_output -s summary_file BFR_ref_nuc_mt_final_v2.fasta.align

# -a designates name of new output align file
# -s designates name of output summary file
```

```{bash eval = FALSE}
head divergence_align_output
```

    1943 16.47 9.01 4.66 Contig_1 2292 2724 (12323366) C rnd-5_family-584#Unknown (21) 512 62 1 m_b1s001i0

      Contig_1            2292 CGACCCGGTATCACGGCAAAGCGTGAAATAGACACGCCGGTCCACTACAG 2341
                                i                                    v
    C rnd-5_family-        512 CAACCCGGTATCACGGCAAAGCGTGAAATAGACACGCCCGTCCACTACAG 463

      Contig_1            2342 CAAGACGTGAAATAACCACACTTTGCGCTGCCAAAGCGTGAAAGTT---- 2387
                                  i      i   i    i     v    i v i           ----
    C rnd-5_family-        462 CAAAACGTGAGATAGCCACGCTTTGGGCTGTCTAGGCGTGAAAGTTACAC 413

```{bash eval = FALSE}
head summary_file
```

    Weighted average Kimura divergence for each repeat family
    Class   Repeat  absLen  wellCharLen     Kimura%
    -----   ------  ------  -----------     -------
    DNA     rnd-1_family-25 338350  333299  11.71
    DNA     rnd-1_family-62 22055   21706   15.46
    DNA     rnd-4_family-337        198511  195143  15.48
    DNA     rnd-4_family-39 280533  265662  15.73

## Plot the repeat landscape

```{bash eval = FALSE}
createRepeatLandscape.pl -div summary_file -g 595951805 > BFR_repeat_landscape.html
```

![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/repeatmasker/custom_repeatmodeler_library/BFR_repeat_landscape_repeatmodeler.png)
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/repeatmasker/custom_repeatmodeler_library/BFR_repeat_landscape_repeatmodeler_2.png)
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/repeatmasker/custom_repeatmodeler_library/BFR_repeat_landscape_repeatmodeler_3.png)
