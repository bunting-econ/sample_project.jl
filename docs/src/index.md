# HomogeneityTestBBU.jl

HomogeneityTestBBU implements the permutation test from Bugni; Bunting and Ura (2023). By default; the test function takes the data ```X``` and a test statistic ```test_stat_fn``` and returns the p-value computed as Equation (4); Bugni; Bunting and Ura (2023). The test rejects the null of homogeneity for all significance levels greater than the p-value.

```@docs
homogeneity_test_fn
```

```@docs
euler_algorithm
```


The following code partially replicates the empirical application in Bugni; Bunting and Ura (2023).

```@setup init
import Random
using sample_project
Random.seed!(1)
S = [ 13  13  13  13  11   8   8   8   7   7   6   6   9   9   9   8   8   9   9;
 15  20  19  19  20  20  14  15  21  20  20  20  21  21  19  19  19  20  20;
 14  14  14  14  14  15  15  15  15  15  15  15  13  13  12  12  12  12  12;
 12  12  12  12  12  12  12  12  12  12  13  13  13  13  13  13  12  13  13;
 32  36  36  37  38  38  38  36  32  33  33  35  36  35  33  33  33  35  34;
 13  18  16  16  16  17  17  17  14  15  15  14  14  15  13  13  13  13  13;
  7  10  10  10  10   8   8   9  11  11  11  10  10  10   9  10  10  10  10;
 17  17  18  18  18  18  18  16  16  16  16  15  15  15  14  14  14  14  14;
 14  13  13  11  11  11  11  12  12  12  12  11  11  11  10  11  11  11  10;
 23  23  22  24  24  24  24  24  23  22  22  20  20  20  19  20  20  20  20;
  8  10  10  10  10  10  10  10  10   8   8   8   8   8   8   8   9   9  10;
 13  13  13  12  12   7   7   9  12  12  12  13  13  13  12  12  12  12  12;
 14  14  14  14  14  14  14  14  14  14  14  14  14  14  13  12  12  13  13;
 11  11  11  11  11  11  11  11  11  11   9   9   9   9   8   8   8   8   8;
 11  11  11   9   9   9   8   9   9   6   6   6   6   6   5   6   6   6   6;
 17  18  18  18  18  17  18  18  18  18  18  18  18  18  16  16  16  16  17;
 28  25  23  25  25  23  23  23  23  23  23  23  24  24  22  22  22  20  20;
 22  21  21  22  21  21  22  22  22  21  22  21  21  22  20  20  20  20  21;
 19  19  15  17  17  17  17  17  17  17  17  17  18  18  14  14  14  14  15;
 11  12  12  12  11  11  11   9   8   8   8   8   8   5   5   5   5   5   5;
 33  32  27  30  30  30  30  31  31  31  31  33  33  34  30  31  30  31  31;
 12  12  12  11  12  12  12  12  12  12  12  13  12  12  11  11  11  11  12;
 47  46  45  51  48  47  47  44  43  39  39  39  39  39  36  36  36  36  36]

function test_fn(S_data, A_data)

    S_support = sort(unique(S_data));
    S_data = hcat(S_data, A_data[:,end]);
    outcome_support = sort(unique(S_data[:, 2:end]))
  
    S_cols = size(S_data, 2);
    
    T_P = Array{Union{Missing, Float64}}(missing, size(S_data,1), size(S_support, 1), size(outcome_support, 1));
    T_P_star = Array{Union{Missing, Float64}}(missing, size(S_data,1),size(S_support, 1), size(outcome_support, 1));
    P_d = Array{Float64}(undef, size(S_support, 1), size(outcome_support, 1));
    
    for j2 in eachindex(S_support)
        for j3 in eachindex(outcome_support)
            tmp0 = (S_data[:,1:(S_cols-1)].==S_support[j2]);
            if sum(tmp0)>0
                P_d[j2,j3] = sum((S_data[:,2:S_cols][tmp0].==outcome_support[j3])) / sum(tmp0);
            end
            if P_d[j2,j3] > 0 # If P_d==0, then tpstar contribution is \divide 0, and wjd is \infty
                for j1 in axes(S_data, 1)
                    tmp1 = (S_data[j1,1:(S_cols-1)].==S_support[j2]);
                    if sum(tmp1)>0 # If tmp1 is zero, then Pjd is \divide 0
                        P_j_d = sum((S_data[j1,2:S_cols][tmp1] .== outcome_support[j3])) ./   sum(tmp1);
                        W_j_d = sum(tmp1) / P_d[j2,j3];
                        T_P[j1,j2,j3] = W_j_d * (P_j_d - P_d[j2,j3])^2;
                        if (P_j_d > 0)
                            T_P_star[j1, j2, j3] = 2 *  W_j_d * P_d[j2,j3] * P_j_d * log(P_j_d / P_d[j2,j3]);
                        end
                    end
                end
            end
        end
    end

    T_P = sum(skipmissing(T_P));
    T_P_star = sum(skipmissing(T_P_star));
    
    return T_P, T_P_star
end

```

```@example init
S_pre = S[:,1:10];
A_pre = S[:,2:11]
TS = test_fn(S_pre, A_pre)
println(TS)
```
```@example init
t = homogeneity_test_fn(X = (S_pre, A_pre),
    test_stat_fn = test_fn,
    K=10000,
    A_trivial=true)
t
```