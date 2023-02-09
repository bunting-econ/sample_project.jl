"""
    homogeneity_test_fn

This function implements the homogeneity test of Bugni, Bunting and Ura (2023).


## Arguments
- `X::Tuple`: The data in the form of tuple ``X=(S, A)``, where A may be empty.
- `test_stat_fn::Function`: A function that takes the arguments ``S`` and ``A`` (or ``S`` only, if ``A=[~]``) and returns a (vector of) test statistics.
- `K::Int`: An integer indicating the length of the MCMC chain. Defaults to 10,000.
- `verbose::Vector{String}`: Possible elements of the vector are `"P value"`, `"Test stat"`, and `"X^K"`. Defaults to empty.
- `A_trivial::Boolean`: A Boolean indicating whether ``S'=A``. Defaults to false.

## Values
- `Pvalue`: The test's p-value (Bugni, Bunting and Ura (2023), equation (4)).
- `Pvalue_chain`: The p-value computed at each of ``K`` steps (optional).
- `test_stat_chain`: The test statistic computed on ``X`` and of ``K`` steps (optional).
- `X_K': The tuple containing the the permuted data (optional).
"""
function homogeneity_test_fn(;
    X::Tuple, ## The data
    test_stat_fn::Function,  ## The test statistic
    K::Int=10000, ## The length of the MCMC chain to generate
    verbose::Vector{String} = [],
    A_trivial::Bool = false
    )

    S_data, A_data = X

    S_only = isempty(A_data)
    
    if !(isempty(A_data) | (size(A_data, 1) > 0))
        println("***Check that X is a tuple***")
    end
    
    if !all(size(S_data) .== size(A_data) ) & !isempty(A_data)
        println("***S and A are not the same dimension***")
    end
    
    if !(K>0)
        println("***Choose K > 0***")
    end

    if S_only
        test_stat = test_stat_fn(S_data)
    else
        test_stat = test_stat_fn(S_data, A_data)
    end
    numb_statistics = length(test_stat)
    
    local_n, local_t = size(S_data);

    SS_list =    hcat(S_data[:,1],S_data[:,2])     
    for t = 2:(local_t-1)
        SS_list = vcat(SS_list, hcat(S_data[:,t],S_data[:,t+1]) );
    end
    SS_list = vcat(SS_list, hcat(S_data[:,end], zeros(Int, local_n)) )
    
    S_data_permuted = copy(S_data);
    if !S_only
        A_data_permuted = copy(A_data);
    end

    perm_TestStat = fill(NaN, K, numb_statistics);
    for k = 1:K
        rand_markets = rand(1:local_n, 2);
        if rand_markets[1]!=rand_markets[2]
            S_data_mixed = [S_data_permuted[rand_markets[1],:] ;0 ;
                                        S_data_permuted[rand_markets[2],:];0];
            S_data_mixed  = euler_algorithm(S_data_mixed);
            deleteat!(S_data_mixed, [local_t+1 ;length(S_data_mixed)])
            S_data_permuted[rand_markets[1],:] = S_data_mixed[1:local_t];
            S_data_permuted[rand_markets[2],:] = S_data_mixed[local_t.+(1:local_t)];
        end

        for i in 1:local_n
            S_data_permuted[i,:] = euler_algorithm(S_data_permuted[i,:]);
        end

        if !S_only
            for ind in axes(unique(SS_list, dims=1), 1)
                ss_ind = unique(SS_list, dims=1)[ind, :]'
                if (A_trivial) || (A_trivial && ss_ind[2]==0)
                    ss_perm_ind = BitArray(undef, size(S_data_permuted))
                    S_perm_temp = hcat(S_data_permuted, zeros(Int, size(S_data_permuted, 1)))
                    for i in axes(ss_perm_ind, 1)
                        for j in axes(ss_perm_ind, 2)
                            ss_perm_ind[i,j] = all(S_perm_temp[i,j .+ [0,1]]' .== ss_ind)
                        end
                    end

                    ss_ind = all(SS_list .== ss_ind, dims= 2)
                    ss_ind = reshape(ss_ind, size(S_data))
                    A_data_permuted[ss_perm_ind] = Random.shuffle(A_data[ss_ind])
                end
            end
            if A_trivial
                A_data_permuted[:,1:(end-1)] = S_data_permuted[:,2:end]
            end
        end
        if S_only
            perm_TestStat[k,:] .= test_stat_fn(S_data_permuted)
        else
            perm_TestStat[k,:] .= test_stat_fn(S_data_permuted, A_data_permuted);
        end
    end

    pvalue = fill(NaN, numb_statistics);
    for ind in 1:numb_statistics
        pvalue[ind] = sum(perm_TestStat[:,ind] .>= test_stat[ind])/K; 
    end

    if "p value" ∈ lowercase.(verbose)
        pvalue_eachk = fill(NaN, K, numb_statistics)
        for ind in 1:numb_statistics
            for ind2 in 1:K
                pvalue_eachk[ind2, ind] = sum(perm_TestStat[1:ind2, ind] .>= test_stat[ind]) / ind2
            end
        end
    end

    t = (Pvalue = pvalue,)
    if "p value" ∈ lowercase.(verbose)
        t =  merge(t, (;Pvalue_chain = pvalue_eachk))
    end
    if "test stat" ∈ first.(lowercase.(verbose), 9)
        t =  merge(t, (;test_stat_chain = vcat([test_stat[i] for i=1:numb_statistics]', perm_TestStat)))
    end
    if "x^k" ∈ lowercase.(verbose)
        t = merge(t, (;X_K = (S_data_permuted, A_data_permuted)))
    end

    return t
end
