# Analysis summary

## TCGA performance (tumour = positive class)

| panel                  | model     | eval                                   |   n |   sensitivity |   specificity |   accuracy |   auc | auc_ci                                   |
|:-----------------------|:----------|:---------------------------------------|----:|--------------:|--------------:|-----------:|------:|:-----------------------------------------|
| blood_compatible_11    | logreg    | TCGA 5x10-fold CV (locked panel)       | 333 |         0.976 |         0.974 |      0.976 | 0.995 | (0.9877765388046388, 0.9997323818019624) |
| blood_compatible_11    | logreg    | TCGA single 70/30 split                | 100 |         1.000 |         1.000 |      1.000 | 1.000 | (1.0, 1.0)                               |
| blood_compatible_11    | mlp       | TCGA 5x10-fold CV (locked panel)       | 333 |         0.976 |         0.974 |      0.976 | 0.996 | (0.9893844781445139, 0.9998215878679749) |
| blood_compatible_11    | mlp       | TCGA single 70/30 split                | 100 |         0.989 |         1.000 |      0.990 | 1.000 | (1.0, 1.0)                               |
| blood_compatible_11    | rf        | TCGA 5x10-fold CV (locked panel)       | 333 |         0.993 |         0.921 |      0.985 | 0.997 | (0.992593666369313, 1.0)                 |
| blood_compatible_11    | rf        | TCGA single 70/30 split                | 100 |         1.000 |         1.000 |      1.000 | 1.000 | (1.0, 1.0)                               |
| blood_compatible_11    | mean_beta | TCGA (no training; threshold = Youden) | 333 |         0.976 |         1.000 |      0.979 | 0.996 | (0.9897390722569134, 0.9997323818019626) |
| screened_11            | logreg    | TCGA 5x10-fold CV (locked panel)       | 333 |         0.986 |         0.974 |      0.985 | 0.999 | (0.9982158786797501, 1.0)                |
| screened_11            | logreg    | TCGA single 70/30 split                | 100 |         1.000 |         1.000 |      1.000 | 1.000 | (1.0, 1.0)                               |
| screened_11            | mlp       | TCGA 5x10-fold CV (locked panel)       | 333 |         0.986 |         1.000 |      0.988 | 1.000 | (0.9983942908117751, 1.0)                |
| screened_11            | mlp       | TCGA single 70/30 split                | 100 |         0.989 |         1.000 |      0.990 | 1.000 | (1.0, 1.0)                               |
| screened_11            | rf        | TCGA 5x10-fold CV (locked panel)       | 333 |         1.000 |         0.947 |      0.994 | 1.000 | (0.9981266726137378, 1.0)                |
| screened_11            | rf        | TCGA single 70/30 split                | 100 |         1.000 |         1.000 |      1.000 | 1.000 | (1.0, 1.0)                               |
| screened_11            | mean_beta | TCGA (no training; threshold = Youden) | 333 |         0.973 |         1.000 |      0.976 | 0.999 | (0.9969669937555753, 1.0)                |
| network_prioritized_11 | logreg    | TCGA 5x10-fold CV (locked panel)       | 333 |         0.953 |         1.000 |      0.958 | 0.992 | (0.984121320249777, 0.9974130240856378)  |
| network_prioritized_11 | logreg    | TCGA single 70/30 split                | 100 |         0.944 |         1.000 |      0.950 | 0.993 | (0.9775280898876404, 1.0)                |
| network_prioritized_11 | mlp       | TCGA 5x10-fold CV (locked panel)       | 333 |         0.963 |         1.000 |      0.967 | 0.997 | (0.9920584299732381, 0.9995539696699376) |
| network_prioritized_11 | mlp       | TCGA single 70/30 split                | 100 |         0.966 |         1.000 |      0.970 | 0.999 | (0.9938712972420838, 1.0)                |
| network_prioritized_11 | rf        | TCGA 5x10-fold CV (locked panel)       | 333 |         0.997 |         0.974 |      0.994 | 1.000 | (0.9994647636039251, 1.0)                |
| network_prioritized_11 | rf        | TCGA single 70/30 split                | 100 |         1.000 |         1.000 |      1.000 | 1.000 | (1.0, 1.0)                               |
| network_prioritized_11 | mean_beta | TCGA (no training; threshold = Youden) | 333 |         0.963 |         0.974 |      0.964 | 0.988 | (0.9775178412132024, 0.995807314897413)  |
| SEPT9_promoter         | logreg    | TCGA 5x10-fold CV (locked panel)       | 333 |         0.969 |         0.974 |      0.970 | 0.996 | (0.9872390722569138, 1.0)                |
| SEPT9_promoter         | logreg    | TCGA single 70/30 split                | 100 |         0.944 |         1.000 |      0.950 | 0.999 | (0.9938712972420838, 1.0)                |
| SEPT9_promoter         | mlp       | TCGA 5x10-fold CV (locked panel)       | 333 |         0.976 |         0.947 |      0.973 | 0.996 | (0.9909009812667262, 0.9996431757359501) |
| SEPT9_promoter         | mlp       | TCGA single 70/30 split                | 100 |         0.944 |         1.000 |      0.950 | 0.998 | (0.9908069458631256, 1.0)                |
| SEPT9_promoter         | rf        | TCGA 5x10-fold CV (locked panel)       | 333 |         0.993 |         0.789 |      0.970 | 0.984 | (0.9557515611061551, 0.9991971454058876) |
| SEPT9_promoter         | rf        | TCGA single 70/30 split                | 100 |         0.989 |         0.818 |      0.970 | 0.997 | (0.9877425944841675, 1.0)                |
| SEPT9_promoter         | mean_beta | TCGA (no training; threshold = Youden) | 333 |         0.780 |         0.974 |      0.802 | 0.852 | (0.8104348795718109, 0.8880463871543265) |

## External validation, frozen TCGA model

| panel                  | model     | dataset   | subset            |   n_pos |   n_neg |   tp |   fn |   tn |   fp |   sensitivity |   specificity |   accuracy |   auc |   unlabelled_called_positive |   unlabelled_n |
|:-----------------------|:----------|:----------|:------------------|--------:|--------:|-----:|-----:|-----:|-----:|--------------:|--------------:|-----------:|------:|-----------------------------:|---------------:|
| blood_compatible_11    | logreg    | GSE48684  | all               |      64 |      41 |   57 |    7 |   40 |    1 |         0.891 |         0.976 |      0.924 | 0.981 |                       37.000 |         42.000 |
| blood_compatible_11    | logreg    | GSE131013 | all               |      96 |     144 |   91 |    5 |  136 |    8 |         0.948 |         0.944 |      0.946 | 0.976 |                      nan     |        nan     |
| blood_compatible_11    | logreg    | GSE131013 | tumor_vs_adjacent |      96 |      96 |   91 |    5 |   90 |    6 |         0.948 |         0.938 |      0.943 | 0.971 |                      nan     |        nan     |
| blood_compatible_11    | logreg    | GSE131013 | tumor_vs_healthy  |      96 |      48 |   91 |    5 |   46 |    2 |         0.948 |         0.958 |      0.951 | 0.985 |                      nan     |        nan     |
| blood_compatible_11    | mlp       | GSE48684  | all               |      64 |      41 |   55 |    9 |   40 |    1 |         0.859 |         0.976 |      0.905 | 0.982 |                       36.000 |         42.000 |
| blood_compatible_11    | mlp       | GSE131013 | all               |      96 |     144 |   91 |    5 |  137 |    7 |         0.948 |         0.951 |      0.950 | 0.972 |                      nan     |        nan     |
| blood_compatible_11    | mlp       | GSE131013 | tumor_vs_adjacent |      96 |      96 |   91 |    5 |   90 |    6 |         0.948 |         0.938 |      0.943 | 0.968 |                      nan     |        nan     |
| blood_compatible_11    | mlp       | GSE131013 | tumor_vs_healthy  |      96 |      48 |   91 |    5 |   47 |    1 |         0.948 |         0.979 |      0.958 | 0.980 |                      nan     |        nan     |
| blood_compatible_11    | rf        | GSE48684  | all               |      64 |      41 |   61 |    3 |   37 |    4 |         0.953 |         0.902 |      0.933 | 0.981 |                       41.000 |         42.000 |
| blood_compatible_11    | rf        | GSE131013 | all               |      96 |     144 |   93 |    3 |  113 |   31 |         0.969 |         0.785 |      0.858 | 0.961 |                      nan     |        nan     |
| blood_compatible_11    | rf        | GSE131013 | tumor_vs_adjacent |      96 |      96 |   93 |    3 |   75 |   21 |         0.969 |         0.781 |      0.875 | 0.955 |                      nan     |        nan     |
| blood_compatible_11    | rf        | GSE131013 | tumor_vs_healthy  |      96 |      48 |   93 |    3 |   38 |   10 |         0.969 |         0.792 |      0.910 | 0.972 |                      nan     |        nan     |
| blood_compatible_11    | mean_beta | GSE48684  | all               |      64 |      41 |   56 |    8 |   40 |    1 |         0.875 |         0.976 |      0.914 | 0.978 |                      nan     |        nan     |
| blood_compatible_11    | mean_beta | GSE131013 | all               |      96 |     144 |   91 |    5 |  136 |    8 |         0.948 |         0.944 |      0.946 | 0.974 |                      nan     |        nan     |
| screened_11            | logreg    | GSE48684  | all               |      64 |      41 |   60 |    4 |   40 |    1 |         0.938 |         0.976 |      0.952 | 0.977 |                       38.000 |         42.000 |
| screened_11            | logreg    | GSE131013 | all               |      96 |     144 |   91 |    5 |  137 |    7 |         0.948 |         0.951 |      0.950 | 0.978 |                      nan     |        nan     |
| screened_11            | logreg    | GSE131013 | tumor_vs_adjacent |      96 |      96 |   91 |    5 |   90 |    6 |         0.948 |         0.938 |      0.943 | 0.973 |                      nan     |        nan     |
| screened_11            | logreg    | GSE131013 | tumor_vs_healthy  |      96 |      48 |   91 |    5 |   47 |    1 |         0.948 |         0.979 |      0.958 | 0.987 |                      nan     |        nan     |
| screened_11            | mlp       | GSE48684  | all               |      64 |      41 |   60 |    4 |   39 |    2 |         0.938 |         0.951 |      0.943 | 0.968 |                       36.000 |         42.000 |
| screened_11            | mlp       | GSE131013 | all               |      96 |     144 |   89 |    7 |  139 |    5 |         0.927 |         0.965 |      0.950 | 0.975 |                      nan     |        nan     |
| screened_11            | mlp       | GSE131013 | tumor_vs_adjacent |      96 |      96 |   89 |    7 |   92 |    4 |         0.927 |         0.958 |      0.943 | 0.970 |                      nan     |        nan     |
| screened_11            | mlp       | GSE131013 | tumor_vs_healthy  |      96 |      48 |   89 |    7 |   47 |    1 |         0.927 |         0.979 |      0.944 | 0.986 |                      nan     |        nan     |
| screened_11            | rf        | GSE48684  | all               |      64 |      41 |   62 |    2 |   22 |   19 |         0.969 |         0.537 |      0.800 | 0.973 |                       41.000 |         42.000 |
| screened_11            | rf        | GSE131013 | all               |      96 |     144 |   92 |    4 |  135 |    9 |         0.958 |         0.938 |      0.946 | 0.969 |                      nan     |        nan     |
| screened_11            | rf        | GSE131013 | tumor_vs_adjacent |      96 |      96 |   92 |    4 |   89 |    7 |         0.958 |         0.927 |      0.943 | 0.965 |                      nan     |        nan     |
| screened_11            | rf        | GSE131013 | tumor_vs_healthy  |      96 |      48 |   92 |    4 |   46 |    2 |         0.958 |         0.958 |      0.958 | 0.975 |                      nan     |        nan     |
| screened_11            | mean_beta | GSE48684  | all               |      64 |      41 |   56 |    8 |   41 |    0 |         0.875 |         1.000 |      0.924 | 0.976 |                      nan     |        nan     |
| screened_11            | mean_beta | GSE131013 | all               |      96 |     144 |   90 |    6 |  137 |    7 |         0.938 |         0.951 |      0.946 | 0.982 |                      nan     |        nan     |
| network_prioritized_11 | logreg    | GSE48684  | all               |      64 |      41 |   61 |    3 |   37 |    4 |         0.953 |         0.902 |      0.933 | 0.976 |                       36.000 |         42.000 |
| network_prioritized_11 | logreg    | GSE131013 | all               |      96 |     144 |   90 |    6 |  135 |    9 |         0.938 |         0.938 |      0.938 | 0.975 |                      nan     |        nan     |
| network_prioritized_11 | logreg    | GSE131013 | tumor_vs_adjacent |      96 |      96 |   90 |    6 |   89 |    7 |         0.938 |         0.927 |      0.932 | 0.970 |                      nan     |        nan     |
| network_prioritized_11 | logreg    | GSE131013 | tumor_vs_healthy  |      96 |      48 |   90 |    6 |   46 |    2 |         0.938 |         0.958 |      0.944 | 0.984 |                      nan     |        nan     |
| network_prioritized_11 | mlp       | GSE48684  | all               |      64 |      41 |   59 |    5 |   39 |    2 |         0.922 |         0.951 |      0.933 | 0.990 |                       39.000 |         42.000 |
| network_prioritized_11 | mlp       | GSE131013 | all               |      96 |     144 |   91 |    5 |  134 |   10 |         0.948 |         0.931 |      0.938 | 0.974 |                      nan     |        nan     |
| network_prioritized_11 | mlp       | GSE131013 | tumor_vs_adjacent |      96 |      96 |   91 |    5 |   89 |    7 |         0.948 |         0.927 |      0.938 | 0.968 |                      nan     |        nan     |
| network_prioritized_11 | mlp       | GSE131013 | tumor_vs_healthy  |      96 |      48 |   91 |    5 |   45 |    3 |         0.948 |         0.938 |      0.944 | 0.987 |                      nan     |        nan     |
| network_prioritized_11 | rf        | GSE48684  | all               |      64 |      41 |   63 |    1 |    5 |   36 |         0.984 |         0.122 |      0.648 | 0.923 |                       41.000 |         42.000 |
| network_prioritized_11 | rf        | GSE131013 | all               |      96 |     144 |   94 |    2 |  128 |   16 |         0.979 |         0.889 |      0.925 | 0.961 |                      nan     |        nan     |
| network_prioritized_11 | rf        | GSE131013 | tumor_vs_adjacent |      96 |      96 |   94 |    2 |   85 |   11 |         0.979 |         0.885 |      0.932 | 0.957 |                      nan     |        nan     |
| network_prioritized_11 | rf        | GSE131013 | tumor_vs_healthy  |      96 |      48 |   94 |    2 |   43 |    5 |         0.979 |         0.896 |      0.951 | 0.967 |                      nan     |        nan     |
| network_prioritized_11 | mean_beta | GSE48684  | all               |      64 |      41 |   59 |    5 |   19 |   22 |         0.922 |         0.463 |      0.743 | 0.924 |                      nan     |        nan     |
| network_prioritized_11 | mean_beta | GSE131013 | all               |      96 |     144 |   90 |    6 |  135 |    9 |         0.938 |         0.938 |      0.938 | 0.965 |                      nan     |        nan     |
| SEPT9_promoter         | logreg    | GSE48684  | all               |      64 |      41 |   63 |    1 |   35 |    6 |         0.984 |         0.854 |      0.933 | 0.968 |                       32.000 |         42.000 |
| SEPT9_promoter         | logreg    | GSE131013 | all               |      96 |     144 |   93 |    3 |  125 |   19 |         0.969 |         0.868 |      0.908 | 0.956 |                      nan     |        nan     |
| SEPT9_promoter         | logreg    | GSE131013 | tumor_vs_adjacent |      96 |      96 |   93 |    3 |   86 |   10 |         0.969 |         0.896 |      0.932 | 0.953 |                      nan     |        nan     |
| SEPT9_promoter         | logreg    | GSE131013 | tumor_vs_healthy  |      96 |      48 |   93 |    3 |   39 |    9 |         0.969 |         0.812 |      0.917 | 0.962 |                      nan     |        nan     |
| SEPT9_promoter         | mlp       | GSE48684  | all               |      64 |      41 |   63 |    1 |   24 |   17 |         0.984 |         0.585 |      0.829 | 0.929 |                       26.000 |         42.000 |
| SEPT9_promoter         | mlp       | GSE131013 | all               |      96 |     144 |   94 |    2 |   99 |   45 |         0.979 |         0.688 |      0.804 | 0.948 |                      nan     |        nan     |
| SEPT9_promoter         | mlp       | GSE131013 | tumor_vs_adjacent |      96 |      96 |   94 |    2 |   75 |   21 |         0.979 |         0.781 |      0.880 | 0.958 |                      nan     |        nan     |
| SEPT9_promoter         | mlp       | GSE131013 | tumor_vs_healthy  |      96 |      48 |   94 |    2 |   24 |   24 |         0.979 |         0.500 |      0.819 | 0.928 |                      nan     |        nan     |
| SEPT9_promoter         | rf        | GSE48684  | all               |      64 |      41 |   64 |    0 |    0 |   41 |         1.000 |         0.000 |      0.610 | 0.912 |                       42.000 |         42.000 |
| SEPT9_promoter         | rf        | GSE131013 | all               |      96 |     144 |   95 |    1 |   85 |   59 |         0.990 |         0.590 |      0.750 | 0.973 |                      nan     |        nan     |
| SEPT9_promoter         | rf        | GSE131013 | tumor_vs_adjacent |      96 |      96 |   95 |    1 |   61 |   35 |         0.990 |         0.635 |      0.812 | 0.968 |                      nan     |        nan     |
| SEPT9_promoter         | rf        | GSE131013 | tumor_vs_healthy  |      96 |      48 |   95 |    1 |   24 |   24 |         0.990 |         0.500 |      0.826 | 0.982 |                      nan     |        nan     |
| SEPT9_promoter         | mean_beta | GSE48684  | all               |      64 |      41 |   63 |    1 |    0 |   41 |         0.984 |         0.000 |      0.600 | 0.843 |                      nan     |        nan     |
| SEPT9_promoter         | mean_beta | GSE131013 | all               |      96 |     144 |   76 |   20 |  135 |    9 |         0.792 |         0.938 |      0.879 | 0.884 |                      nan     |        nan     |

## Nested CV (whole selection funnel redone inside each outer fold)

- with blood step, logreg: AUC 0.994 (95% CI 0.984, 1.000), sens 97.6%, spec 97.4%, comparator-CpG overlap per fold [1, 1, 1, 1, 1]
- with blood step, mlp: AUC 0.996 (95% CI 0.989, 1.000), sens 96.3%, spec 97.4%, comparator-CpG overlap per fold [1, 1, 1, 1, 1]
- without blood step, logreg: AUC 0.999 (95% CI 0.997, 1.000), sens 98.6%, spec 97.4%, comparator-CpG overlap per fold [0, 0, 1, 1, 1]
- without blood step, mlp: AUC 0.999 (95% CI 0.997, 1.000), sens 98.0%, spec 100.0%, comparator-CpG overlap per fold [0, 0, 1, 1, 1]

## Random-panel benchmark

Null median AUC 0.996 (5–95%: 0.983–1.000) over 1000 random panels from 19841 hyper DMCs.
- SEPT9_promoter: AUC 0.995, empirical p = 0.559
- blood_compatible_11: AUC 0.996, empirical p = 0.504
- network_prioritized_11: AUC 0.991, empirical p = 0.756
- screened_11: AUC 0.999, empirical p = 0.156
## Traceability of the 11 comparator CpGs

| probe      | annotation                                  | island   | zhou_masked   | A_passes_qc   |   A_delta |       A_fdr | A_status   | A_primary_candidate   | A_extended_candidate   | B_passes_qc   |   B_delta |       B_fdr | B_status   | B_primary_candidate   | B_extended_candidate   | B_extended_specific   | B_extended_excluded_by   |   A_primary_specific | A_primary_excluded_by   |   A_extended_specific | A_extended_excluded_by   |   B_primary_specific | B_primary_excluded_by   |
|:-----------|:--------------------------------------------|:---------|:--------------|:--------------|----------:|------------:|:-----------|:----------------------|:-----------------------|:--------------|----------:|------------:|:-----------|:----------------------|:-----------------------|:----------------------|:-------------------------|---------------------:|:------------------------|----------------------:|:-------------------------|---------------------:|:------------------------|
| cg11909912 | LOC100130148:Body;MAPT:5'UTR                | Island   | False         | True          |    0.3254 | 1.92358e-39 | hyper      | False                 | False                  | True          |    0.2669 | 1.18818e-17 | hyper      | False                 | True                   | True                  |                          |                  nan | nan                     |                   nan | nan                      |                  nan | nan                     |
| cg02450613 | GDNF:TSS1500                                | Island   | False         | True          |    0.4228 | 1.30317e-38 | hyper      | True                  | True                   | True          |    0.4277 | 4.40758e-17 | hyper      | True                  | True                   | True                  |                          |                    1 |                         |                     1 |                          |                    1 |                         |
| cg09980522 | GRIA4:1stExon;GRIA4:Body                    | S_Shore  | False         | True          |    0.3976 | 3.55487e-32 | hyper      | False                 | True                   | True          |    0.3395 | 5.66576e-12 | hyper      | False                 | True                   | True                  |                          |                  nan | nan                     |                     1 |                          |                  nan | nan                     |
| cg06999323 | NEFL:TSS1500                                | S_Shore  | False         | True          |    0.3383 | 3.00825e-34 | hyper      | True                  | True                   | True          |    0.2954 | 9.86418e-13 | hyper      | True                  | True                   | True                  |                          |                    1 |                         |                     1 |                          |                    1 |                         |
| cg15032098 | UCHL1:TSS1500                               | N_Shore  | False         | True          |    0.2999 | 2.05173e-26 | hyper      | True                  | True                   | True          |    0.2985 | 1.24552e-12 | hyper      | True                  | True                   | True                  |                          |                    1 |                         |                     1 |                          |                    1 |                         |
| cg05659265 | GNAO1:1stExon;GNAO1:5'UTR;LOC283856:TSS1500 | Island   | False         | True          |    0.2302 | 6.87937e-23 | hyper      | False                 | True                   | True          |    0.2383 | 4.92949e-11 | hyper      | False                 | True                   | True                  |                          |                  nan | nan                     |                     1 |                          |                  nan | nan                     |
| cg08976810 | FGF10:TSS1500                               | OpenSea  | False         | True          |    0.4084 | 8.64439e-38 | hyper      | True                  | True                   | True          |    0.219  | 1.06774e-13 | hyper      | True                  | True                   | True                  |                          |                    1 |                         |                     1 |                          |                    1 |                         |
| cg03203223 | RELN:TSS1500                                | Island   | False         | True          |    0.4099 | 5.78755e-35 | hyper      | True                  | True                   | True          |    0.4107 | 6.20743e-15 | hyper      | True                  | True                   | True                  |                          |                    1 |                         |                     1 |                          |                    1 |                         |
| cg16368442 | GRIN2A:TSS1500;GRIN2A:TSS200                | Island   | False         | True          |    0.4988 | 1.70718e-42 | hyper      | True                  | True                   | True          |    0.4759 | 2.88518e-19 | hyper      | True                  | True                   | True                  |                          |                    1 |                         |                     1 |                          |                    1 |                         |
| cg01035160 | SNCA:5'UTR;SNCA:TSS1500;SNCA:TSS200         | Island   | False         | True          |    0.2965 | 2.14179e-13 | hyper      | True                  | True                   | True          |    0.2906 | 0.000103963 | hyper      | True                  | True                   | True                  |                          |                    1 |                         |                     1 |                          |                    1 |                         |
| cg03742003 | SLIT2:1stExon;SLIT2:5'UTR                   | Island   | False         | True          |    0.3525 | 4.30124e-39 | hyper      | False                 | True                   | True          |    0.3176 | 5.48947e-17 | hyper      | False                 | True                   | True                  |                          |                  nan | nan                     |                     1 |                          |                  nan | nan                     |

Determinism check (two identical runs): True

## Blood background and plasma cfDNA

Whole blood GSE40279 (n = 656; 328 set the rule, the rest are simulated), sorted leukocytes GSE35069 (n = 60), plasma cfDNA GSE122126 EPIC (n = 90). Blood-compatible: blood β ≤ 0.1 and tumour median − blood β ≥ 0.2.

- CRC-specific CpGs passing: 238/943
- comparator CpGs passing: 3/11 (cg02450613, cg03203223, cg16368442)
- comparator CpGs absent from EPIC: none

### Network-prioritized 11 CpGs

|            | genes   |   tcga_tumor_median |   tcga_normal_median |   wb_median |   wb_p95 |   celltype_max_median | celltype_argmax   |   blood_beta |   margin | blood_compatible   |   moss_healthy_cfdna_median |   moss_colon_cancer_cfdna_median |
|:-----------|:--------|--------------------:|---------------------:|------------:|---------:|----------------------:|:------------------|-------------:|---------:|:-------------------|----------------------------:|---------------------------------:|
| cg01035160 | SNCA    |               0.475 |                0.091 |       0.068 |    0.109 |                 0.116 | CD19+ B cells     |        0.116 |    0.359 | False              |                       0.058 |                            0.505 |
| cg02450613 | GDNF    |               0.484 |                0.026 |       0.022 |    0.043 |                 0.049 | CD19+ B cells     |        0.049 |    0.435 | True               |                       0.113 |                            0.119 |
| cg03203223 | RELN    |               0.500 |                0.015 |       0.011 |    0.035 |                 0.036 | CD19+ B cells     |        0.036 |    0.464 | True               |                       0.017 |                            0.363 |
| cg03742003 | SLIT2   |               0.532 |                0.218 |       0.215 |    0.270 |                 0.282 | CD8+ T cells      |        0.282 |    0.250 | False              |                       0.149 |                            0.478 |
| cg05659265 | GNAO1   |               0.324 |                0.117 |       0.164 |    0.206 |                 0.170 | Neutrophils       |        0.206 |    0.118 | False              |                       0.121 |                            0.185 |
| cg06999323 | NEFL    |               0.448 |                0.159 |       0.074 |    0.116 |                 0.111 | CD19+ B cells     |        0.116 |    0.332 | False              |                       0.081 |                            0.359 |
| cg08976810 | FGF10   |               0.630 |                0.404 |       0.083 |    0.120 |                 0.119 | CD19+ B cells     |        0.120 |    0.509 | False              |                       0.065 |                            0.725 |
| cg09980522 | GRIA4   |               0.619 |                0.230 |       0.062 |    0.104 |                 0.123 | CD8+ T cells      |        0.123 |    0.496 | False              |                       0.128 |                            0.445 |
| cg11909912 | MAPT    |               0.582 |                0.306 |       0.185 |    0.248 |                 0.241 | CD8+ T cells      |        0.248 |    0.334 | False              |                       0.107 |                            0.714 |
| cg15032098 | UCHL1   |               0.547 |                0.231 |       0.159 |    0.222 |                 0.203 | CD8+ T cells      |        0.222 |    0.325 | False              |                       0.144 |                            0.387 |
| cg16368442 | GRIN2A  |               0.566 |                0.044 |       0.008 |    0.023 |                 0.018 | CD8+ T cells      |        0.023 |    0.543 | True               |                       0.013 |                            0.401 |

### In-silico dilution: lowest tumour fraction with ≥80% sensitivity at 95% specificity

| panel                                |   depth |   lowest_fraction_detected |
|:-------------------------------------|--------:|---------------------------:|
| SEPT9_promoter                       |       0 |                     nan    |
| SEPT9_promoter                       |     100 |                     nan    |
| blood_compatible_11                  |       0 |                       0.1  |
| blood_compatible_11                  |     100 |                       0.1  |
| network_prioritized_11               |       0 |                       0.1  |
| network_prioritized_11               |     100 |                       0.1  |
| screened_11                          |       0 |                       0.1  |
| screened_11                          |     100 |                       0.1  |
| network_prioritized_blood_compatible |       0 |                       0.05 |
| network_prioritized_blood_compatible |     100 |                       0.1  |

### Dilution grid

| panel                                |   n_cpgs |   depth |   tumor_fraction |   auc |   sensitivity_at_spec |   threshold |   case_mean |   control_mean |
|:-------------------------------------|---------:|--------:|-----------------:|------:|----------------------:|------------:|------------:|---------------:|
| SEPT9_promoter                       |       46 |       0 |            0.005 | 0.446 |                 0.023 |       0.369 |       0.349 |          0.352 |
| SEPT9_promoter                       |       46 |       0 |            0.010 | 0.445 |                 0.022 |       0.369 |       0.349 |          0.352 |
| SEPT9_promoter                       |       46 |       0 |            0.020 | 0.444 |                 0.022 |       0.369 |       0.349 |          0.352 |
| SEPT9_promoter                       |       46 |       0 |            0.050 | 0.440 |                 0.021 |       0.369 |       0.349 |          0.352 |
| SEPT9_promoter                       |       46 |       0 |            0.100 | 0.438 |                 0.022 |       0.369 |       0.349 |          0.352 |
| SEPT9_promoter                       |       46 |       0 |            0.200 | 0.452 |                 0.031 |       0.369 |       0.349 |          0.352 |
| SEPT9_promoter                       |       46 |     100 |            0.005 | 0.456 |                 0.026 |       0.371 |       0.350 |          0.352 |
| SEPT9_promoter                       |       46 |     100 |            0.010 | 0.454 |                 0.021 |       0.371 |       0.349 |          0.352 |
| SEPT9_promoter                       |       46 |     100 |            0.020 | 0.437 |                 0.018 |       0.371 |       0.349 |          0.352 |
| SEPT9_promoter                       |       46 |     100 |            0.050 | 0.442 |                 0.020 |       0.371 |       0.349 |          0.352 |
| SEPT9_promoter                       |       46 |     100 |            0.100 | 0.442 |                 0.019 |       0.371 |       0.349 |          0.352 |
| SEPT9_promoter                       |       46 |     100 |            0.200 | 0.463 |                 0.037 |       0.371 |       0.349 |          0.352 |
| blood_compatible_11                  |       11 |       0 |            0.005 | 0.608 |                 0.053 |       0.048 |       0.031 |          0.028 |
| blood_compatible_11                  |       11 |       0 |            0.010 | 0.682 |                 0.078 |       0.048 |       0.033 |          0.028 |
| blood_compatible_11                  |       11 |       0 |            0.020 | 0.796 |                 0.180 |       0.048 |       0.039 |          0.028 |
| blood_compatible_11                  |       11 |       0 |            0.050 | 0.946 |                 0.737 |       0.048 |       0.054 |          0.028 |
| blood_compatible_11                  |       11 |       0 |            0.100 | 0.989 |                 0.949 |       0.048 |       0.081 |          0.028 |
| blood_compatible_11                  |       11 |       0 |            0.200 | 0.998 |                 0.988 |       0.048 |       0.133 |          0.028 |
| blood_compatible_11                  |       11 |     100 |            0.005 | 0.596 |                 0.050 |       0.050 |       0.031 |          0.028 |
| blood_compatible_11                  |       11 |     100 |            0.010 | 0.662 |                 0.081 |       0.050 |       0.034 |          0.028 |
| blood_compatible_11                  |       11 |     100 |            0.020 | 0.767 |                 0.141 |       0.050 |       0.039 |          0.028 |
| blood_compatible_11                  |       11 |     100 |            0.050 | 0.928 |                 0.574 |       0.050 |       0.054 |          0.028 |
| blood_compatible_11                  |       11 |     100 |            0.100 | 0.984 |                 0.925 |       0.050 |       0.080 |          0.028 |
| blood_compatible_11                  |       11 |     100 |            0.200 | 0.997 |                 0.985 |       0.050 |       0.133 |          0.028 |
| network_prioritized_11               |       11 |       0 |            0.005 | 0.606 |                 0.073 |       0.119 |       0.102 |          0.098 |
| network_prioritized_11               |       11 |       0 |            0.010 | 0.659 |                 0.095 |       0.119 |       0.104 |          0.098 |
| network_prioritized_11               |       11 |       0 |            0.020 | 0.752 |                 0.179 |       0.119 |       0.108 |          0.098 |
| network_prioritized_11               |       11 |       0 |            0.050 | 0.917 |                 0.529 |       0.119 |       0.120 |          0.098 |
| network_prioritized_11               |       11 |       0 |            0.100 | 0.984 |                 0.922 |       0.119 |       0.140 |          0.098 |
| network_prioritized_11               |       11 |       0 |            0.200 | 0.998 |                 0.985 |       0.119 |       0.181 |          0.098 |
| network_prioritized_11               |       11 |     100 |            0.005 | 0.590 |                 0.097 |       0.120 |       0.102 |          0.098 |
| network_prioritized_11               |       11 |     100 |            0.010 | 0.633 |                 0.118 |       0.120 |       0.104 |          0.098 |
| network_prioritized_11               |       11 |     100 |            0.020 | 0.701 |                 0.194 |       0.120 |       0.108 |          0.098 |
| network_prioritized_11               |       11 |     100 |            0.050 | 0.871 |                 0.483 |       0.120 |       0.120 |          0.098 |
| network_prioritized_11               |       11 |     100 |            0.100 | 0.963 |                 0.844 |       0.120 |       0.140 |          0.098 |
| network_prioritized_11               |       11 |     100 |            0.200 | 0.996 |                 0.982 |       0.120 |       0.181 |          0.098 |
| screened_11                          |       11 |       0 |            0.005 | 0.573 |                 0.096 |       0.088 |       0.064 |          0.060 |
| screened_11                          |       11 |       0 |            0.010 | 0.621 |                 0.106 |       0.088 |       0.066 |          0.060 |
| screened_11                          |       11 |       0 |            0.020 | 0.708 |                 0.147 |       0.088 |       0.072 |          0.060 |
| screened_11                          |       11 |       0 |            0.050 | 0.888 |                 0.434 |       0.088 |       0.087 |          0.060 |
| screened_11                          |       11 |       0 |            0.100 | 0.975 |                 0.898 |       0.088 |       0.114 |          0.060 |
| screened_11                          |       11 |       0 |            0.200 | 0.996 |                 0.979 |       0.088 |       0.166 |          0.060 |
| screened_11                          |       11 |     100 |            0.005 | 0.569 |                 0.088 |       0.089 |       0.064 |          0.060 |
| screened_11                          |       11 |     100 |            0.010 | 0.606 |                 0.109 |       0.089 |       0.066 |          0.060 |
| screened_11                          |       11 |     100 |            0.020 | 0.698 |                 0.153 |       0.089 |       0.072 |          0.060 |
| screened_11                          |       11 |     100 |            0.050 | 0.871 |                 0.425 |       0.089 |       0.088 |          0.060 |
| screened_11                          |       11 |     100 |            0.100 | 0.968 |                 0.865 |       0.089 |       0.114 |          0.060 |
| screened_11                          |       11 |     100 |            0.200 | 0.995 |                 0.977 |       0.089 |       0.166 |          0.060 |
| network_prioritized_blood_compatible |        3 |       0 |            0.005 | 0.679 |                 0.103 |       0.027 |       0.019 |          0.014 |
| network_prioritized_blood_compatible |        3 |       0 |            0.010 | 0.762 |                 0.155 |       0.027 |       0.021 |          0.014 |
| network_prioritized_blood_compatible |        3 |       0 |            0.020 | 0.873 |                 0.357 |       0.027 |       0.026 |          0.014 |
| network_prioritized_blood_compatible |        3 |       0 |            0.050 | 0.966 |                 0.873 |       0.027 |       0.039 |          0.014 |
| network_prioritized_blood_compatible |        3 |       0 |            0.100 | 0.984 |                 0.951 |       0.027 |       0.062 |          0.014 |
| network_prioritized_blood_compatible |        3 |       0 |            0.200 | 0.989 |                 0.970 |       0.027 |       0.108 |          0.014 |
| network_prioritized_blood_compatible |        3 |     100 |            0.005 | 0.614 |                 0.076 |       0.033 |       0.019 |          0.014 |
| network_prioritized_blood_compatible |        3 |     100 |            0.010 | 0.677 |                 0.114 |       0.033 |       0.021 |          0.014 |
| network_prioritized_blood_compatible |        3 |     100 |            0.020 | 0.772 |                 0.217 |       0.033 |       0.026 |          0.014 |
| network_prioritized_blood_compatible |        3 |     100 |            0.050 | 0.921 |                 0.598 |       0.033 |       0.040 |          0.014 |
| network_prioritized_blood_compatible |        3 |     100 |            0.100 | 0.971 |                 0.888 |       0.033 |       0.063 |          0.014 |
| network_prioritized_blood_compatible |        3 |     100 |            0.200 | 0.984 |                 0.953 |       0.033 |       0.109 |          0.014 |

### Real plasma cfDNA (z vs healthy cfDNA)

| panel                                |   n_cpgs_on_array | title                        | group               |   spike_fraction |   colon_cell_fraction |   score |   z_vs_healthy | above_all_healthy   |
|:-------------------------------------|------------------:|:-----------------------------|:--------------------|-----------------:|----------------------:|--------:|---------------:|:--------------------|
| SEPT9_promoter                       |                45 | cfDNA_young_male_pool_rep1   | healthy_cfdna       |          nan     |               nan     |   0.352 |          0.668 | False               |
| SEPT9_promoter                       |                45 | cfDNA_young_male_pool_rep2   | healthy_cfdna       |          nan     |               nan     |   0.350 |          0.195 | False               |
| SEPT9_promoter                       |                45 | cfDNA_young_female_pool_rep1 | healthy_cfdna       |          nan     |               nan     |   0.356 |          1.480 | False               |
| SEPT9_promoter                       |                45 | cfDNA_young_female_pool_rep2 | healthy_cfdna       |          nan     |               nan     |   0.352 |          0.651 | False               |
| SEPT9_promoter                       |                45 | cfDNA_old_male_pool_rep1     | healthy_cfdna       |          nan     |               nan     |   0.349 |         -0.008 | False               |
| SEPT9_promoter                       |                45 | cfDNA_old_male_pool_rep2     | healthy_cfdna       |          nan     |               nan     |   0.345 |         -0.939 | False               |
| SEPT9_promoter                       |                45 | cfDNA_old_female_pool_rep1   | healthy_cfdna       |          nan     |               nan     |   0.351 |          0.296 | False               |
| SEPT9_promoter                       |                45 | cfDNA_old_female_pool_rep2   | healthy_cfdna       |          nan     |               nan     |   0.341 |         -1.762 | False               |
| SEPT9_promoter                       |                45 | cfDNA_healthy_individual_1   | healthy_cfdna       |          nan     |               nan     |   0.353 |          0.681 | False               |
| SEPT9_promoter                       |                45 | cfDNA_colon_cancer_806       | colon_cancer_cfdna  |          nan     |               nan     |   0.374 |          5.124 | True                |
| SEPT9_promoter                       |                45 | cfDNA_colon_cancer_684       | colon_cancer_cfdna  |          nan     |               nan     |   0.362 |          2.602 | True                |
| SEPT9_promoter                       |                45 | cfDNA_colon_cancer_806_50ng  | colon_cancer_cfdna  |          nan     |               nan     |   0.366 |          3.428 | True                |
| SEPT9_promoter                       |                45 | cfDNA_lung_cancer_1          | lung_cancer_cfdna   |          nan     |               nan     |   0.339 |         -2.016 | False               |
| SEPT9_promoter                       |                45 | cfDNA_mix_5                  | cfdna_spike_in      |            0.100 |               nan     |   0.373 |          4.830 | True                |
| SEPT9_promoter                       |                45 | cfDNA_mix_6                  | cfdna_spike_in      |            0.050 |               nan     |   0.361 |          2.450 | True                |
| SEPT9_promoter                       |                45 | cfDNA_mix_7                  | cfdna_spike_in      |            0.030 |               nan     |   0.335 |         -2.830 | False               |
| SEPT9_promoter                       |                45 | cfDNA_healthy_individual_2   | healthy_cfdna       |          nan     |               nan     |   0.348 |         -0.318 | False               |
| SEPT9_promoter                       |                45 | cfDNA_lung_cancer_1478       | lung_cancer_cfdna   |          nan     |               nan     |   0.330 |         -3.901 | False               |
| SEPT9_promoter                       |                45 | cfDNA_colon_cancer_2010      | colon_cancer_cfdna  |          nan     |               nan     |   0.351 |          0.326 | False               |
| SEPT9_promoter                       |                45 | cfDNA_breast_cancer_1040     | breast_cancer_cfdna |          nan     |               nan     |   0.370 |          4.255 | True                |
| SEPT9_promoter                       |                45 | cfDNA_breast_cancer_1711     | breast_cancer_cfdna |          nan     |               nan     |   0.273 |        -15.534 | False               |
| SEPT9_promoter                       |                45 | cfDNA_mix_1                  | cfdna_spike_in      |            0.050 |               nan     |   0.368 |          3.726 | True                |
| SEPT9_promoter                       |                45 | cfDNA_healthy_individual_3   | healthy_cfdna       |          nan     |               nan     |   0.341 |         -1.658 | False               |
| SEPT9_promoter                       |                45 | cfDNA_mix_3                  | cfdna_spike_in      |            0.050 |               nan     |   0.369 |          3.960 | True                |
| SEPT9_promoter                       |                45 | cfDNA_healthy_individual_4   | healthy_cfdna       |          nan     |               nan     |   0.353 |          0.715 | False               |
| SEPT9_promoter                       |                45 | cfDNA_breast_cancer_1438     | breast_cancer_cfdna |          nan     |               nan     |   0.495 |         29.735 | True                |
| SEPT9_promoter                       |                45 | cfDNA_lung_cancer_2244       | lung_cancer_cfdna   |          nan     |               nan     |   0.341 |         -1.735 | False               |
| SEPT9_promoter                       |                45 | cfDNA_lung_cancer_1433       | lung_cancer_cfdna   |          nan     |               nan     |   0.335 |         -2.952 | False               |
| SEPT9_promoter                       |                45 | in_vitro_mix_9               | gdna_mix            |          nan     |                 0.000 |   0.346 |         -0.660 | False               |
| SEPT9_promoter                       |                45 | in_vitro_mix_10              | gdna_mix            |          nan     |                 0.000 |   0.346 |         -0.584 | False               |
| SEPT9_promoter                       |                45 | in_vitro_mix_11              | gdna_mix            |          nan     |                 0.050 |   0.334 |         -3.170 | False               |
| SEPT9_promoter                       |                45 | in_vitro_mix_12              | gdna_mix            |          nan     |                 0.100 |   0.332 |         -3.457 | False               |
| SEPT9_promoter                       |                45 | in_vitro_mix_13              | gdna_mix            |          nan     |                 0.035 |   0.333 |         -3.305 | False               |
| SEPT9_promoter                       |                45 | in_vitro_mix_14              | gdna_mix            |          nan     |                 0.000 |   0.333 |         -3.380 | False               |
| SEPT9_promoter                       |                45 | in_vitro_mix_17              | gdna_mix            |          nan     |                 0.060 |   0.327 |         -4.513 | False               |
| SEPT9_promoter                       |                45 | in_vitro_mix_16              | gdna_mix            |          nan     |                 0.000 |   0.333 |         -3.301 | False               |
| SEPT9_promoter                       |                45 | in_vitro_mix_15              | gdna_mix            |          nan     |                 0.000 |   0.333 |         -3.208 | False               |
| blood_compatible_11                  |                10 | cfDNA_young_male_pool_rep1   | healthy_cfdna       |          nan     |               nan     |   0.045 |         -0.010 | False               |
| blood_compatible_11                  |                10 | cfDNA_young_male_pool_rep2   | healthy_cfdna       |          nan     |               nan     |   0.054 |          0.659 | False               |
| blood_compatible_11                  |                10 | cfDNA_young_female_pool_rep1 | healthy_cfdna       |          nan     |               nan     |   0.027 |         -1.298 | False               |
| blood_compatible_11                  |                10 | cfDNA_young_female_pool_rep2 | healthy_cfdna       |          nan     |               nan     |   0.063 |          1.279 | False               |
| blood_compatible_11                  |                10 | cfDNA_old_male_pool_rep1     | healthy_cfdna       |          nan     |               nan     |   0.034 |         -0.737 | False               |
| blood_compatible_11                  |                10 | cfDNA_old_male_pool_rep2     | healthy_cfdna       |          nan     |               nan     |   0.047 |          0.124 | False               |
| blood_compatible_11                  |                10 | cfDNA_old_female_pool_rep1   | healthy_cfdna       |          nan     |               nan     |   0.037 |         -0.582 | False               |
| blood_compatible_11                  |                10 | cfDNA_old_female_pool_rep2   | healthy_cfdna       |          nan     |               nan     |   0.025 |         -1.436 | False               |
| blood_compatible_11                  |                10 | cfDNA_healthy_individual_1   | healthy_cfdna       |          nan     |               nan     |   0.072 |          1.911 | False               |
| blood_compatible_11                  |                10 | cfDNA_colon_cancer_806       | colon_cancer_cfdna  |          nan     |               nan     |   0.553 |         35.921 | True                |
| blood_compatible_11                  |                10 | cfDNA_colon_cancer_684       | colon_cancer_cfdna  |          nan     |               nan     |   0.709 |         46.976 | True                |
| blood_compatible_11                  |                10 | cfDNA_colon_cancer_806_50ng  | colon_cancer_cfdna  |          nan     |               nan     |   0.501 |         32.287 | True                |
| blood_compatible_11                  |                10 | cfDNA_lung_cancer_1          | lung_cancer_cfdna   |          nan     |               nan     |   0.049 |          0.307 | False               |
| blood_compatible_11                  |                10 | cfDNA_mix_5                  | cfdna_spike_in      |            0.100 |               nan     |   0.132 |          6.172 | True                |
| blood_compatible_11                  |                10 | cfDNA_mix_6                  | cfdna_spike_in      |            0.050 |               nan     |   0.084 |          2.744 | True                |
| blood_compatible_11                  |                10 | cfDNA_mix_7                  | cfdna_spike_in      |            0.030 |               nan     |   0.042 |         -0.177 | False               |
| blood_compatible_11                  |                10 | cfDNA_healthy_individual_2   | healthy_cfdna       |          nan     |               nan     |   0.037 |         -0.563 | False               |
| blood_compatible_11                  |                10 | cfDNA_lung_cancer_1478       | lung_cancer_cfdna   |          nan     |               nan     |   0.029 |         -1.126 | False               |
| blood_compatible_11                  |                10 | cfDNA_colon_cancer_2010      | colon_cancer_cfdna  |          nan     |               nan     |   0.057 |          0.885 | False               |
| blood_compatible_11                  |                10 | cfDNA_breast_cancer_1040     | breast_cancer_cfdna |          nan     |               nan     |   0.123 |          5.542 | True                |
| blood_compatible_11                  |                10 | cfDNA_breast_cancer_1711     | breast_cancer_cfdna |          nan     |               nan     |   0.029 |         -1.090 | False               |
| blood_compatible_11                  |                10 | cfDNA_mix_1                  | cfdna_spike_in      |            0.050 |               nan     |   0.107 |          4.362 | True                |
| blood_compatible_11                  |                10 | cfDNA_healthy_individual_3   | healthy_cfdna       |          nan     |               nan     |   0.046 |          0.101 | False               |
| blood_compatible_11                  |                10 | cfDNA_mix_3                  | cfdna_spike_in      |            0.050 |               nan     |   0.087 |          2.947 | True                |
| blood_compatible_11                  |                10 | cfDNA_healthy_individual_4   | healthy_cfdna       |          nan     |               nan     |   0.053 |          0.552 | False               |
| blood_compatible_11                  |                10 | cfDNA_breast_cancer_1438     | breast_cancer_cfdna |          nan     |               nan     |   0.183 |          9.744 | True                |
| blood_compatible_11                  |                10 | cfDNA_lung_cancer_2244       | lung_cancer_cfdna   |          nan     |               nan     |   0.030 |         -1.087 | False               |
| blood_compatible_11                  |                10 | cfDNA_lung_cancer_1433       | lung_cancer_cfdna   |          nan     |               nan     |   0.035 |         -0.698 | False               |
| blood_compatible_11                  |                10 | in_vitro_mix_9               | gdna_mix            |          nan     |                 0.000 |   0.056 |          0.819 | False               |
| blood_compatible_11                  |                10 | in_vitro_mix_10              | gdna_mix            |          nan     |                 0.000 |   0.052 |          0.485 | False               |
| blood_compatible_11                  |                10 | in_vitro_mix_11              | gdna_mix            |          nan     |                 0.050 |   0.070 |          1.766 | False               |
| blood_compatible_11                  |                10 | in_vitro_mix_12              | gdna_mix            |          nan     |                 0.100 |   0.060 |          1.061 | False               |
| blood_compatible_11                  |                10 | in_vitro_mix_13              | gdna_mix            |          nan     |                 0.035 |   0.055 |          0.706 | False               |
| blood_compatible_11                  |                10 | in_vitro_mix_14              | gdna_mix            |          nan     |                 0.000 |   0.063 |          1.280 | False               |
| blood_compatible_11                  |                10 | in_vitro_mix_17              | gdna_mix            |          nan     |                 0.060 |   0.058 |          0.955 | False               |
| blood_compatible_11                  |                10 | in_vitro_mix_16              | gdna_mix            |          nan     |                 0.000 |   0.056 |          0.817 | False               |
| blood_compatible_11                  |                10 | in_vitro_mix_15              | gdna_mix            |          nan     |                 0.000 |   0.054 |          0.610 | False               |
| network_prioritized_11               |                11 | cfDNA_young_male_pool_rep1   | healthy_cfdna       |          nan     |               nan     |   0.086 |         -0.795 | False               |
| network_prioritized_11               |                11 | cfDNA_young_male_pool_rep2   | healthy_cfdna       |          nan     |               nan     |   0.103 |          0.474 | False               |
| network_prioritized_11               |                11 | cfDNA_young_female_pool_rep1 | healthy_cfdna       |          nan     |               nan     |   0.084 |         -0.974 | False               |
| network_prioritized_11               |                11 | cfDNA_young_female_pool_rep2 | healthy_cfdna       |          nan     |               nan     |   0.102 |          0.351 | False               |
| network_prioritized_11               |                11 | cfDNA_old_male_pool_rep1     | healthy_cfdna       |          nan     |               nan     |   0.089 |         -0.578 | False               |
| network_prioritized_11               |                11 | cfDNA_old_male_pool_rep2     | healthy_cfdna       |          nan     |               nan     |   0.092 |         -0.383 | False               |
| network_prioritized_11               |                11 | cfDNA_old_female_pool_rep1   | healthy_cfdna       |          nan     |               nan     |   0.089 |         -0.576 | False               |
| network_prioritized_11               |                11 | cfDNA_old_female_pool_rep2   | healthy_cfdna       |          nan     |               nan     |   0.077 |         -1.476 | False               |
| network_prioritized_11               |                11 | cfDNA_healthy_individual_1   | healthy_cfdna       |          nan     |               nan     |   0.121 |          1.820 | False               |
| network_prioritized_11               |                11 | cfDNA_colon_cancer_806       | colon_cancer_cfdna  |          nan     |               nan     |   0.443 |         25.771 | True                |
| network_prioritized_11               |                11 | cfDNA_colon_cancer_684       | colon_cancer_cfdna  |          nan     |               nan     |   0.602 |         37.618 | True                |
| network_prioritized_11               |                11 | cfDNA_colon_cancer_806_50ng  | colon_cancer_cfdna  |          nan     |               nan     |   0.446 |         26.011 | True                |
| network_prioritized_11               |                11 | cfDNA_lung_cancer_1          | lung_cancer_cfdna   |          nan     |               nan     |   0.109 |          0.906 | False               |
| network_prioritized_11               |                11 | cfDNA_mix_5                  | cfdna_spike_in      |            0.100 |               nan     |   0.168 |          5.269 | True                |
| network_prioritized_11               |                11 | cfDNA_mix_6                  | cfdna_spike_in      |            0.050 |               nan     |   0.112 |          1.135 | False               |
| network_prioritized_11               |                11 | cfDNA_mix_7                  | cfdna_spike_in      |            0.030 |               nan     |   0.093 |         -0.333 | False               |
| network_prioritized_11               |                11 | cfDNA_healthy_individual_2   | healthy_cfdna       |          nan     |               nan     |   0.095 |         -0.173 | False               |
| network_prioritized_11               |                11 | cfDNA_lung_cancer_1478       | lung_cancer_cfdna   |          nan     |               nan     |   0.090 |         -0.488 | False               |
| network_prioritized_11               |                11 | cfDNA_colon_cancer_2010      | colon_cancer_cfdna  |          nan     |               nan     |   0.168 |          5.332 | True                |
| network_prioritized_11               |                11 | cfDNA_breast_cancer_1040     | breast_cancer_cfdna |          nan     |               nan     |   0.318 |         16.473 | True                |
| network_prioritized_11               |                11 | cfDNA_breast_cancer_1711     | breast_cancer_cfdna |          nan     |               nan     |   0.136 |          2.906 | True                |
| network_prioritized_11               |                11 | cfDNA_mix_1                  | cfdna_spike_in      |            0.050 |               nan     |   0.164 |          5.034 | True                |
| network_prioritized_11               |                11 | cfDNA_healthy_individual_3   | healthy_cfdna       |          nan     |               nan     |   0.115 |          1.344 | False               |
| network_prioritized_11               |                11 | cfDNA_mix_3                  | cfdna_spike_in      |            0.050 |               nan     |   0.125 |          2.081 | True                |
| network_prioritized_11               |                11 | cfDNA_healthy_individual_4   | healthy_cfdna       |          nan     |               nan     |   0.110 |          0.966 | False               |
| network_prioritized_11               |                11 | cfDNA_breast_cancer_1438     | breast_cancer_cfdna |          nan     |               nan     |   0.317 |         16.377 | True                |
| network_prioritized_11               |                11 | cfDNA_lung_cancer_2244       | lung_cancer_cfdna   |          nan     |               nan     |   0.093 |         -0.282 | False               |
| network_prioritized_11               |                11 | cfDNA_lung_cancer_1433       | lung_cancer_cfdna   |          nan     |               nan     |   0.087 |         -0.776 | False               |
| network_prioritized_11               |                11 | in_vitro_mix_9               | gdna_mix            |          nan     |                 0.000 |   0.108 |          0.792 | False               |
| network_prioritized_11               |                11 | in_vitro_mix_10              | gdna_mix            |          nan     |                 0.000 |   0.101 |          0.291 | False               |
| network_prioritized_11               |                11 | in_vitro_mix_11              | gdna_mix            |          nan     |                 0.050 |   0.114 |          1.246 | False               |
| network_prioritized_11               |                11 | in_vitro_mix_12              | gdna_mix            |          nan     |                 0.100 |   0.118 |          1.542 | False               |
| network_prioritized_11               |                11 | in_vitro_mix_13              | gdna_mix            |          nan     |                 0.035 |   0.109 |          0.877 | False               |
| network_prioritized_11               |                11 | in_vitro_mix_14              | gdna_mix            |          nan     |                 0.000 |   0.105 |          0.604 | False               |
| network_prioritized_11               |                11 | in_vitro_mix_17              | gdna_mix            |          nan     |                 0.060 |   0.123 |          1.931 | True                |
| network_prioritized_11               |                11 | in_vitro_mix_16              | gdna_mix            |          nan     |                 0.000 |   0.105 |          0.615 | False               |
| network_prioritized_11               |                11 | in_vitro_mix_15              | gdna_mix            |          nan     |                 0.000 |   0.107 |          0.773 | False               |
| screened_11                          |                10 | cfDNA_young_male_pool_rep1   | healthy_cfdna       |          nan     |               nan     |   0.069 |         -0.161 | False               |
| screened_11                          |                10 | cfDNA_young_male_pool_rep2   | healthy_cfdna       |          nan     |               nan     |   0.074 |          0.115 | False               |
| screened_11                          |                10 | cfDNA_young_female_pool_rep1 | healthy_cfdna       |          nan     |               nan     |   0.055 |         -1.033 | False               |
| screened_11                          |                10 | cfDNA_young_female_pool_rep2 | healthy_cfdna       |          nan     |               nan     |   0.087 |          0.900 | False               |
| screened_11                          |                10 | cfDNA_old_male_pool_rep1     | healthy_cfdna       |          nan     |               nan     |   0.055 |         -1.038 | False               |
| screened_11                          |                10 | cfDNA_old_male_pool_rep2     | healthy_cfdna       |          nan     |               nan     |   0.070 |         -0.104 | False               |
| screened_11                          |                10 | cfDNA_old_female_pool_rep1   | healthy_cfdna       |          nan     |               nan     |   0.052 |         -1.205 | False               |
| screened_11                          |                10 | cfDNA_old_female_pool_rep2   | healthy_cfdna       |          nan     |               nan     |   0.055 |         -1.013 | False               |
| screened_11                          |                10 | cfDNA_healthy_individual_1   | healthy_cfdna       |          nan     |               nan     |   0.108 |          2.222 | False               |
| screened_11                          |                10 | cfDNA_colon_cancer_806       | colon_cancer_cfdna  |          nan     |               nan     |   0.635 |         34.071 | True                |
| screened_11                          |                10 | cfDNA_colon_cancer_684       | colon_cancer_cfdna  |          nan     |               nan     |   0.687 |         37.186 | True                |
| screened_11                          |                10 | cfDNA_colon_cancer_806_50ng  | colon_cancer_cfdna  |          nan     |               nan     |   0.616 |         32.902 | True                |
| screened_11                          |                10 | cfDNA_lung_cancer_1          | lung_cancer_cfdna   |          nan     |               nan     |   0.079 |          0.424 | False               |
| screened_11                          |                10 | cfDNA_mix_5                  | cfdna_spike_in      |            0.100 |               nan     |   0.187 |          6.982 | True                |
| screened_11                          |                10 | cfDNA_mix_6                  | cfdna_spike_in      |            0.050 |               nan     |   0.134 |          3.761 | True                |
| screened_11                          |                10 | cfDNA_mix_7                  | cfdna_spike_in      |            0.030 |               nan     |   0.073 |          0.083 | False               |
| screened_11                          |                10 | cfDNA_healthy_individual_2   | healthy_cfdna       |          nan     |               nan     |   0.079 |          0.415 | False               |
| screened_11                          |                10 | cfDNA_lung_cancer_1478       | lung_cancer_cfdna   |          nan     |               nan     |   0.051 |         -1.242 | False               |
| screened_11                          |                10 | cfDNA_colon_cancer_2010      | colon_cancer_cfdna  |          nan     |               nan     |   0.092 |          1.201 | False               |
| screened_11                          |                10 | cfDNA_breast_cancer_1040     | breast_cancer_cfdna |          nan     |               nan     |   0.222 |          9.068 | True                |
| screened_11                          |                10 | cfDNA_breast_cancer_1711     | breast_cancer_cfdna |          nan     |               nan     |   0.048 |         -1.420 | False               |
| screened_11                          |                10 | cfDNA_mix_1                  | cfdna_spike_in      |            0.050 |               nan     |   0.164 |          5.589 | True                |
| screened_11                          |                10 | cfDNA_healthy_individual_3   | healthy_cfdna       |          nan     |               nan     |   0.079 |          0.467 | False               |
| screened_11                          |                10 | cfDNA_mix_3                  | cfdna_spike_in      |            0.050 |               nan     |   0.129 |          3.487 | True                |
| screened_11                          |                10 | cfDNA_healthy_individual_4   | healthy_cfdna       |          nan     |               nan     |   0.079 |          0.434 | False               |
| screened_11                          |                10 | cfDNA_breast_cancer_1438     | breast_cancer_cfdna |          nan     |               nan     |   0.187 |          6.962 | True                |
| screened_11                          |                10 | cfDNA_lung_cancer_2244       | lung_cancer_cfdna   |          nan     |               nan     |   0.057 |         -0.874 | False               |
| screened_11                          |                10 | cfDNA_lung_cancer_1433       | lung_cancer_cfdna   |          nan     |               nan     |   0.059 |         -0.784 | False               |
| screened_11                          |                10 | in_vitro_mix_9               | gdna_mix            |          nan     |                 0.000 |   0.099 |          1.643 | False               |
| screened_11                          |                10 | in_vitro_mix_10              | gdna_mix            |          nan     |                 0.000 |   0.089 |          1.052 | False               |
| screened_11                          |                10 | in_vitro_mix_11              | gdna_mix            |          nan     |                 0.050 |   0.114 |          2.582 | True                |
| screened_11                          |                10 | in_vitro_mix_12              | gdna_mix            |          nan     |                 0.100 |   0.099 |          1.633 | False               |
| screened_11                          |                10 | in_vitro_mix_13              | gdna_mix            |          nan     |                 0.035 |   0.095 |          1.399 | False               |
| screened_11                          |                10 | in_vitro_mix_14              | gdna_mix            |          nan     |                 0.000 |   0.106 |          2.083 | False               |
| screened_11                          |                10 | in_vitro_mix_17              | gdna_mix            |          nan     |                 0.060 |   0.101 |          1.743 | False               |
| screened_11                          |                10 | in_vitro_mix_16              | gdna_mix            |          nan     |                 0.000 |   0.101 |          1.779 | False               |
| screened_11                          |                10 | in_vitro_mix_15              | gdna_mix            |          nan     |                 0.000 |   0.104 |          1.976 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_young_male_pool_rep1   | healthy_cfdna       |          nan     |               nan     |   0.057 |          0.389 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_young_male_pool_rep2   | healthy_cfdna       |          nan     |               nan     |   0.067 |          1.211 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_young_female_pool_rep1 | healthy_cfdna       |          nan     |               nan     |   0.045 |         -0.713 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_young_female_pool_rep2 | healthy_cfdna       |          nan     |               nan     |   0.063 |          0.842 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_old_male_pool_rep1     | healthy_cfdna       |          nan     |               nan     |   0.044 |         -0.800 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_old_male_pool_rep2     | healthy_cfdna       |          nan     |               nan     |   0.048 |         -0.401 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_old_female_pool_rep1   | healthy_cfdna       |          nan     |               nan     |   0.049 |         -0.363 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_old_female_pool_rep2   | healthy_cfdna       |          nan     |               nan     |   0.029 |         -2.076 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_healthy_individual_1   | healthy_cfdna       |          nan     |               nan     |   0.070 |          1.459 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_colon_cancer_806       | colon_cancer_cfdna  |          nan     |               nan     |   0.295 |         20.987 | True                |
| network_prioritized_blood_compatible |                 3 | cfDNA_colon_cancer_684       | colon_cancer_cfdna  |          nan     |               nan     |   0.592 |         46.764 | True                |
| network_prioritized_blood_compatible |                 3 | cfDNA_colon_cancer_806_50ng  | colon_cancer_cfdna  |          nan     |               nan     |   0.258 |         17.814 | True                |
| network_prioritized_blood_compatible |                 3 | cfDNA_lung_cancer_1          | lung_cancer_cfdna   |          nan     |               nan     |   0.044 |         -0.762 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_mix_5                  | cfdna_spike_in      |            0.100 |               nan     |   0.082 |          2.507 | True                |
| network_prioritized_blood_compatible |                 3 | cfDNA_mix_6                  | cfdna_spike_in      |            0.050 |               nan     |   0.054 |          0.127 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_mix_7                  | cfdna_spike_in      |            0.030 |               nan     |   0.038 |         -1.320 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_healthy_individual_2   | healthy_cfdna       |          nan     |               nan     |   0.051 |         -0.188 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_lung_cancer_1478       | lung_cancer_cfdna   |          nan     |               nan     |   0.046 |         -0.590 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_colon_cancer_2010      | colon_cancer_cfdna  |          nan     |               nan     |   0.042 |         -0.929 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_breast_cancer_1040     | breast_cancer_cfdna |          nan     |               nan     |   0.069 |          1.381 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_breast_cancer_1711     | breast_cancer_cfdna |          nan     |               nan     |   0.022 |         -2.703 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_mix_1                  | cfdna_spike_in      |            0.050 |               nan     |   0.080 |          2.383 | True                |
| network_prioritized_blood_compatible |                 3 | cfDNA_healthy_individual_3   | healthy_cfdna       |          nan     |               nan     |   0.063 |          0.846 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_mix_3                  | cfdna_spike_in      |            0.050 |               nan     |   0.057 |          0.375 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_healthy_individual_4   | healthy_cfdna       |          nan     |               nan     |   0.050 |         -0.207 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_breast_cancer_1438     | breast_cancer_cfdna |          nan     |               nan     |   0.283 |         20.016 | True                |
| network_prioritized_blood_compatible |                 3 | cfDNA_lung_cancer_2244       | lung_cancer_cfdna   |          nan     |               nan     |   0.044 |         -0.733 | False               |
| network_prioritized_blood_compatible |                 3 | cfDNA_lung_cancer_1433       | lung_cancer_cfdna   |          nan     |               nan     |   0.042 |         -0.912 | False               |
| network_prioritized_blood_compatible |                 3 | in_vitro_mix_9               | gdna_mix            |          nan     |                 0.000 |   0.039 |         -1.191 | False               |
| network_prioritized_blood_compatible |                 3 | in_vitro_mix_10              | gdna_mix            |          nan     |                 0.000 |   0.041 |         -1.053 | False               |
| network_prioritized_blood_compatible |                 3 | in_vitro_mix_11              | gdna_mix            |          nan     |                 0.050 |   0.043 |         -0.854 | False               |
| network_prioritized_blood_compatible |                 3 | in_vitro_mix_12              | gdna_mix            |          nan     |                 0.100 |   0.040 |         -1.094 | False               |
| network_prioritized_blood_compatible |                 3 | in_vitro_mix_13              | gdna_mix            |          nan     |                 0.035 |   0.039 |         -1.185 | False               |
| network_prioritized_blood_compatible |                 3 | in_vitro_mix_14              | gdna_mix            |          nan     |                 0.000 |   0.035 |         -1.532 | False               |
| network_prioritized_blood_compatible |                 3 | in_vitro_mix_17              | gdna_mix            |          nan     |                 0.060 |   0.049 |         -0.336 | False               |
| network_prioritized_blood_compatible |                 3 | in_vitro_mix_16              | gdna_mix            |          nan     |                 0.000 |   0.036 |         -1.419 | False               |
| network_prioritized_blood_compatible |                 3 | in_vitro_mix_15              | gdna_mix            |          nan     |                 0.000 |   0.042 |         -0.922 | False               |

### cfDNA score by group

| panel                                | group                  |   count |   median |   min |   max |
|:-------------------------------------|:-----------------------|--------:|---------:|------:|------:|
| SEPT9_promoter                       | breast_cancer_cfdna    |       3 |    0.370 | 0.273 | 0.495 |
| SEPT9_promoter                       | cfdna_spike_in         |       5 |    0.368 | 0.335 | 0.373 |
| SEPT9_promoter                       | colon_cancer_cfdna     |       4 |    0.364 | 0.351 | 0.374 |
| SEPT9_promoter                       | colon_epithelium       |       3 |    0.207 | 0.206 | 0.222 |
| SEPT9_promoter                       | cup_cfdna              |       4 |    0.351 | 0.323 | 0.457 |
| SEPT9_promoter                       | gdna_mix               |       9 |    0.333 | 0.327 | 0.346 |
| SEPT9_promoter                       | healthy_cfdna          |      12 |    0.350 | 0.341 | 0.356 |
| SEPT9_promoter                       | islet_transplant_cfdna |       9 |    0.335 | 0.328 | 0.370 |
| SEPT9_promoter                       | leukocytes             |       1 |    0.326 | 0.326 | 0.326 |
| SEPT9_promoter                       | lung_cancer_cfdna      |       4 |    0.337 | 0.330 | 0.341 |
| SEPT9_promoter                       | other                  |      14 |    0.351 | 0.246 | 0.496 |
| SEPT9_promoter                       | sepsis_cfdna           |      22 |    0.342 | 0.318 | 0.450 |
| blood_compatible_11                  | breast_cancer_cfdna    |       3 |    0.123 | 0.029 | 0.183 |
| blood_compatible_11                  | cfdna_spike_in         |       5 |    0.087 | 0.042 | 0.132 |
| blood_compatible_11                  | colon_cancer_cfdna     |       4 |    0.527 | 0.057 | 0.709 |
| blood_compatible_11                  | colon_epithelium       |       3 |    0.091 | 0.046 | 0.144 |
| blood_compatible_11                  | cup_cfdna              |       4 |    0.073 | 0.038 | 0.173 |
| blood_compatible_11                  | gdna_mix               |       9 |    0.056 | 0.052 | 0.070 |
| blood_compatible_11                  | healthy_cfdna          |      12 |    0.046 | 0.025 | 0.072 |
| blood_compatible_11                  | islet_transplant_cfdna |       9 |    0.022 | 0.018 | 0.031 |
| blood_compatible_11                  | leukocytes             |       1 |    0.062 | 0.062 | 0.062 |
| blood_compatible_11                  | lung_cancer_cfdna      |       4 |    0.032 | 0.029 | 0.049 |
| blood_compatible_11                  | other                  |      14 |    0.031 | 0.021 | 0.054 |
| blood_compatible_11                  | sepsis_cfdna           |      22 |    0.046 | 0.037 | 0.331 |
| network_prioritized_11               | breast_cancer_cfdna    |       3 |    0.317 | 0.136 | 0.318 |
| network_prioritized_11               | cfdna_spike_in         |       5 |    0.125 | 0.093 | 0.168 |
| network_prioritized_11               | colon_cancer_cfdna     |       4 |    0.444 | 0.168 | 0.602 |
| network_prioritized_11               | colon_epithelium       |       3 |    0.191 | 0.172 | 0.288 |
| network_prioritized_11               | cup_cfdna              |       4 |    0.148 | 0.108 | 0.189 |
| network_prioritized_11               | gdna_mix               |       9 |    0.108 | 0.101 | 0.123 |
| network_prioritized_11               | healthy_cfdna          |      12 |    0.093 | 0.077 | 0.121 |
| network_prioritized_11               | islet_transplant_cfdna |       9 |    0.088 | 0.083 | 0.096 |
| network_prioritized_11               | leukocytes             |       1 |    0.106 | 0.106 | 0.106 |
| network_prioritized_11               | lung_cancer_cfdna      |       4 |    0.092 | 0.087 | 0.109 |
| network_prioritized_11               | other                  |      14 |    0.080 | 0.054 | 0.096 |
| network_prioritized_11               | sepsis_cfdna           |      22 |    0.108 | 0.092 | 0.348 |
| network_prioritized_blood_compatible | breast_cancer_cfdna    |       3 |    0.069 | 0.022 | 0.283 |
| network_prioritized_blood_compatible | cfdna_spike_in         |       5 |    0.057 | 0.038 | 0.082 |
| network_prioritized_blood_compatible | colon_cancer_cfdna     |       4 |    0.276 | 0.042 | 0.592 |
| network_prioritized_blood_compatible | colon_epithelium       |       3 |    0.054 | 0.027 | 0.169 |
| network_prioritized_blood_compatible | cup_cfdna              |       4 |    0.080 | 0.048 | 0.135 |
| network_prioritized_blood_compatible | gdna_mix               |       9 |    0.040 | 0.035 | 0.049 |
| network_prioritized_blood_compatible | healthy_cfdna          |      12 |    0.051 | 0.029 | 0.070 |
| network_prioritized_blood_compatible | islet_transplant_cfdna |       9 |    0.029 | 0.022 | 0.042 |
| network_prioritized_blood_compatible | leukocytes             |       1 |    0.036 | 0.036 | 0.036 |
| network_prioritized_blood_compatible | lung_cancer_cfdna      |       4 |    0.044 | 0.042 | 0.046 |
| network_prioritized_blood_compatible | other                  |      14 |    0.020 | 0.015 | 0.026 |
| network_prioritized_blood_compatible | sepsis_cfdna           |      22 |    0.047 | 0.035 | 0.303 |
| screened_11                          | breast_cancer_cfdna    |       3 |    0.187 | 0.048 | 0.222 |
| screened_11                          | cfdna_spike_in         |       5 |    0.134 | 0.073 | 0.187 |
| screened_11                          | colon_cancer_cfdna     |       4 |    0.625 | 0.092 | 0.687 |
| screened_11                          | colon_epithelium       |       3 |    0.091 | 0.044 | 0.137 |
| screened_11                          | cup_cfdna              |       4 |    0.131 | 0.071 | 0.255 |
| screened_11                          | gdna_mix               |       9 |    0.101 | 0.089 | 0.114 |
| screened_11                          | healthy_cfdna          |      12 |    0.072 | 0.052 | 0.108 |
| screened_11                          | islet_transplant_cfdna |       9 |    0.040 | 0.034 | 0.059 |
| screened_11                          | leukocytes             |       1 |    0.111 | 0.111 | 0.111 |
| screened_11                          | lung_cancer_cfdna      |       4 |    0.058 | 0.051 | 0.079 |
| screened_11                          | other                  |      14 |    0.049 | 0.034 | 0.072 |
| screened_11                          | sepsis_cfdna           |      22 |    0.073 | 0.064 | 0.355 |

### Blood-compatible panels in tissue

| panel                                | model   | dataset                       |   n |   sensitivity |   specificity |   accuracy |   auc |
|:-------------------------------------|:--------|:------------------------------|----:|--------------:|--------------:|-----------:|------:|
| network_prioritized_blood_compatible | logreg  | TCGA CV (locked panel)        | 333 |         0.942 |         1.000 |      0.949 | 0.979 |
| network_prioritized_blood_compatible | logreg  | GSE48684 (frozen TCGA model)  | 105 |         0.812 |         1.000 |      0.886 | 0.965 |
| network_prioritized_blood_compatible | logreg  | GSE131013 (frozen TCGA model) | 240 |         0.969 |         0.931 |      0.946 | 0.974 |
| network_prioritized_blood_compatible | mlp     | TCGA CV (locked panel)        | 333 |         0.936 |         1.000 |      0.943 | 0.978 |
| network_prioritized_blood_compatible | mlp     | GSE48684 (frozen TCGA model)  | 105 |         0.828 |         1.000 |      0.895 | 0.964 |
| network_prioritized_blood_compatible | mlp     | GSE131013 (frozen TCGA model) | 240 |         0.938 |         0.944 |      0.942 | 0.975 |
| blood_compatible_11                  | logreg  | TCGA CV (locked panel)        | 333 |         0.976 |         0.974 |      0.976 | 0.995 |
| blood_compatible_11                  | logreg  | GSE48684 (frozen TCGA model)  | 105 |         0.891 |         0.976 |      0.924 | 0.981 |
| blood_compatible_11                  | logreg  | GSE131013 (frozen TCGA model) | 240 |         0.948 |         0.944 |      0.946 | 0.976 |
| blood_compatible_11                  | mlp     | TCGA CV (locked panel)        | 333 |         0.976 |         0.974 |      0.976 | 0.996 |
| blood_compatible_11                  | mlp     | GSE48684 (frozen TCGA model)  | 105 |         0.859 |         0.976 |      0.905 | 0.982 |
| blood_compatible_11                  | mlp     | GSE131013 (frozen TCGA model) | 240 |         0.948 |         0.951 |      0.950 | 0.972 |

Panels: ['cg18672939', 'cg02500300', 'cg11855526', 'cg11666087', 'cg08413157', 'cg23977631', 'cg16368442', 'cg19898128', 'cg10002178', 'cg08808128', 'cg01808545']
