"""
    euler_algorithm

This function applies a Euler algorithm to the First Order Markov chain `x`

## Arguments
- `x::vector`
"""
function euler_algorithm(seq_old)
    seq_length = length(seq_old);

    mid_zero = seq_old[seq_length ÷ 2] == 0;

    # v2: include stop for algorithm as soon as known that first 0 not in middle of sequence
    if !mid_zero
        zeta = seq_old[end];
        zeta = [zeta; rand( repeat(seq_old[seq_old[[2:seq_length;1]] .== zeta], 2), 1)];
        
        while length(unique(zeta[2:end])) < length(unique(seq_old))
            zeta = [zeta; rand( repeat(seq_old[seq_old[[2:seq_length;1]] .== zeta[length(zeta)]], 2), 1)];
        end 
        zeta
    
        N_count = hcat(seq_old[1:(seq_length-1)], seq_old[2:seq_length]);
        exit_path = hcat(unique(seq_old), Vector{Union{Missing, Int8}}(missing, length(unique(seq_old))) );
        for s in axes(exit_path, 1)
            exit_path[s, 2] = zeta[findfirst(zeta[2:end] .== exit_path[s,1])];
            if exit_path[s,1] != seq_old[end]
                N_count = N_count[1:end .!= findfirst((exit_path[s,1] .== N_count[:,1]) .& (N_count[:,2] .== exit_path[s,2])),:];
            end
        end

        seq_new = seq_old[1]
        for i in 2:seq_length
            if sum(N_count[:,1] .== seq_new[i-1]) >= 1
                seq_new = [seq_new; rand(repeat(N_count[N_count[:,1].==seq_new[i-1],2],2), 1)];
                N_count = N_count[1:end .!= minimum((1:size(N_count,1))[(seq_new[i-1].==N_count[:,1]).&(seq_new[i].==N_count[:,2])]),:];
            else
                seq_new = [seq_new; exit_path[exit_path[:,1] .== seq_new[i-1], 2]];
            end
            # println(seq_new)
        end
    else ## this is for 0 in middle
        seq_new = seq_old[1]
        while ((seq_new[end] !=0) | (length(seq_new) < (seq_length ÷ 2)))
            zeta = seq_old[end];
            zeta = [zeta; rand( repeat(seq_old[seq_old[[2:seq_length;1]] .== zeta], 2), 1)];
            
            while length(unique(zeta[2:end])) < length(unique(seq_old))
                zeta = [zeta; rand( repeat(seq_old[seq_old[[2:seq_length;1]] .== zeta[length(zeta)]], 2), 1)];
            end 
        
            N_count = hcat(seq_old[1:(seq_length-1)], seq_old[2:seq_length]);
            exit_path = hcat(unique(seq_old), Vector{Union{Missing, Int8}}(missing, length(unique(seq_old))) );
            for s in axes(exit_path, 1)
                exit_path[s, 2] = zeta[findfirst(zeta[2:end] .== exit_path[s,1])];
                if exit_path[s,1] != seq_old[end]
                    N_count = N_count[1:end .!= findfirst((exit_path[s,1] .== N_count[:,1]) .& (N_count[:,2] .== exit_path[s,2])),:];
                end
            end
            
            seq_new = seq_old[1];
            i = 2;
            while ((i <= (seq_length ÷ 2)) & !any(seq_new .== 0))
                if sum(N_count[:,1] .== seq_new[i-1]) >= 1
                    seq_new = [seq_new; rand(repeat(N_count[N_count[:,1].==seq_new[i-1],2],2), 1)];
                    # N_count = N_count[1:end .!= minimum((1:size(N_count,1))[(seq_new[i-1].==N_count[:,1]).&(seq_new[i].==N_count[:,2])]),:];
                    N_count = N_count[1:end .!= findfirst((seq_new[i-1] .== N_count[:,1]) .& (seq_new[i] .== N_count[:,2])),:];
                else
                    seq_new = [seq_new; exit_path[exit_path[:,1] .== seq_new[i-1], 2]];
                end
                # println(seq_new)
                i = i + 1;    
                # println(((i <= (seq_length ÷ 2)) & !any(seq_new .== 0)))
            end
            # println(((seq_new[end] !=0) | (length(seq_new) < (seq_length ÷ 2))))
        end
        for i in (seq_length ÷ 2 +1):seq_length
            if sum(N_count[:,1] .== seq_new[i-1]) >= 1
                seq_new = [seq_new; rand(repeat(N_count[N_count[:,1].==seq_new[i-1],2],2), 1)];
                # N_count = N_count[1:end .!= minimum((1:size(N_count,1))[(seq_new[i-1].==N_count[:,1]).&(seq_new[i].==N_count[:,2])]),:];
                N_count = N_count[1:end .!= findfirst((seq_new[i-1] .== N_count[:,1]) .& (seq_new[i] .== N_count[:,2])),:];
            else
                seq_new = [seq_new; exit_path[exit_path[:,1] .== seq_new[i-1], 2]];
            end
        end
    end

    
    return seq_new
end