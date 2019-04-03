using MatrixNetworks

function run_hungarian(mat_original) 
    mat = copy(mat_original)  
    setup_time_start = time()  

    # subtract col minimum
    for c=1:size(mat)[2]
        mat[:,c] .-= minimum(mat[:,c])
    end

    # subtract row minimum
    for r=1:size(mat)[1]
        mat[r,:] .-= minimum(mat[r,:])
    end

    nrows = size(mat)[1]
    ncols = size(mat)[2]
    rows_marked = zeros(Bool,nrows)
    cols_marked = zeros(Bool,ncols)

    new_marked_row_ind = falses(nrows)
    new_marked_col_ind = falses(ncols)

    ind_marked_rows = zeros(Int64, nrows)
    ind_unmarked_rows = zeros(Int64, nrows)
    ind_marked_cols = zeros(Int64, ncols)
    ind_unmarked_cols = zeros(Int64, ncols)

    inner_while_total = 0.0
    matzeros_total_time = 0.0
    create_ind_vect_total_time = 0.0
    find_min_total_time = 0.0
    add_min_total_time = 0.0
    sub_min_total_time = 0.0
    zero_vals_total_time = 0.0
    
    matzeros = findall(iszero, mat)
    nmatzeros = length(matzeros)
    max_nmatzeros = nmatzeros

    c_outer_while = 1
    c_inner_while = 1

    println("Time for setup: ", time()-setup_time_start)

    while true
        @views matzeros_sub = matzeros[1:nmatzeros]
        new_marked_col_ind[:] .= false
        new_marked_row_ind[:] .= false

        c_outer_while += 1
       
        matching  = get_matching(matzeros_sub)

        if matching.cardinality == ncols
            println("Total inner while time: ", inner_while_total)
            println("Total matzeros_total_time: ", matzeros_total_time)
            println("Total create_ind_vect_total_time: ", create_ind_vect_total_time)
            println("Total find_min_total_time: ", find_min_total_time)
            println("Total sub_min_total_time: ", sub_min_total_time)
            println("Total add_min_total_time: ", add_min_total_time)
            println("Total zero_vals_total_time: ", zero_vals_total_time)
            println("c_outer_while: ", c_outer_while)
            println("c_inner_while: ", c_inner_while)
            return matching.match
        else
            sorted_matching = sort(matching.match)[1:matching.cardinality]

            # mark all rows which aren't matched
            mark_start = time()
            
            rows_marked[:] .= false
            cols_marked[:] .= false
            last = 0
            for r=1:matching.cardinality
                if sorted_matching[r]-last != 1
                    for ri=last+1:sorted_matching[r]-1
                        rows_marked[ri] = true
                        new_marked_row_ind[ri] = true
                    end
                    last = sorted_matching[r]
                else
                    last = sorted_matching[r]
                end
            end
            for r=last+1:nrows
                rows_marked[r] = true
                new_marked_row_ind[r] = true
            end
        end
       
        start_inner = time()
        changed = true
        while changed
            c_inner_while += 1
            changed = false
            
            # mark cols
            @inbounds for rc in matzeros_sub
                @inbounds if !cols_marked[rc[2]]
                    @inbounds if new_marked_row_ind[rc[1]]
                        cols_marked[rc[2]] = true
                        new_marked_col_ind[rc[2]] = true
                    end
                end
            end

            # mark new rows
            @inbounds for c in 1:ncols
                r = matching.match[c]
                @inbounds if new_marked_col_ind[c] && !rows_marked[r]
                    rows_marked[r] = true
                    new_marked_row_ind[r] = true
                    changed = true
                end
            end
        end
        inner_while_total += time()-start_inner
        
        start_timer = time()
        nmatzeros = 0 
        for matzero in matzeros_sub
            if (rows_marked[matzero[1]] && cols_marked[matzero[2]]) || (!rows_marked[matzero[1]] && !cols_marked[matzero[2]])
                nmatzeros += 1
                if nmatzeros > max_nmatzeros
                    max_nmatzeros += 1
                    push!(matzeros, matzero)
                else
                    matzeros[nmatzeros] = matzero
                end
            end
        end
        matzeros_total_time += time()-start_timer

        start_timer = time()
        ind_marked_rows_end = 0
        ind_unmarked_rows_end = 0
        for r=1:nrows
            if rows_marked[r]
                ind_marked_rows_end += 1
                ind_marked_rows[ind_marked_rows_end] = r
            else
                ind_unmarked_rows_end += 1
                ind_unmarked_rows[ind_unmarked_rows_end] = r
            end
        end
        
        ind_marked_cols_end = 0
        ind_unmarked_cols_end = 0
        for c=1:ncols
            if cols_marked[c]
                ind_marked_cols_end += 1
                ind_marked_cols[ind_marked_cols_end] = c
            else
                ind_unmarked_cols_end += 1
                ind_unmarked_cols[ind_unmarked_cols_end] = c
            end
        end
        
        @views ind_marked_rows_sub = ind_marked_rows[1:ind_marked_rows_end]
        @views ind_unmarked_rows_sub = ind_unmarked_rows[1:ind_unmarked_rows_end]
        @views ind_marked_cols_sub = ind_marked_cols[1:ind_marked_cols_end]
        @views ind_unmarked_cols_sub = ind_unmarked_cols[1:ind_unmarked_cols_end]
        create_ind_vect_total_time += time()-start_timer
        
        start_timer = time()
        #minimum where !cols_marked but rows_marked
        @inbounds @views min_val = minimum(mat[ind_marked_rows_sub,ind_unmarked_cols_sub])
        find_min_total_time += time()-start_timer
        
        start_timer = time()
        # subtract minimum from this mask
        @inbounds for c in ind_unmarked_cols_sub
            for r in ind_marked_rows_sub
                mat[r,c] -= min_val
                if mat[r,c] == 0
                    nmatzeros += 1
                    if nmatzeros > max_nmatzeros
                        max_nmatzeros += 1
                        push!(matzeros, CartesianIndex(r,c))
                    else
                        matzeros[nmatzeros] = CartesianIndex(r,c)
                    end
                end
            end
        end
        sub_min_total_time += time()-start_timer

        start_timer = time()        
        # add minimum where !rows_marked but cols_marked
        @inbounds @views mat[ind_unmarked_rows_sub,ind_marked_cols_sub] .+= min_val
        add_min_total_time += time()-start_timer
    end
end

function get_matching(matzeros_sub)
    # create bi partite matching graph
    ei = Vector{Int64}()
    ej = Vector{Int64}()

    for yx in matzeros_sub
        push!(ej,yx[1]) # row
        push!(ei,yx[2]) # col
    end
    
    matching = bipartite_matching(ones(Int64,length(ei)), ei, ej)
    return matching
end

function main()
    mat = [1 2 3 3 1;
           2 3 5 2 2;
           3 1 4 1 5;
           4 5 2 4 3;
           5 4 1 5 4]
    matching = hungarian(mat)
    @show matching
end

# main()