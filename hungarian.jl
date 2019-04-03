using MatrixNetworks

function run_hungarian(mat_original) 
    mat = copy(mat_original)    

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
    
    c_outer_while = 1
    c_inner_while = 1

    while true
        c_outer_while += 1
        matching = get_matching(mat)

        if matching.cardinality == ncols
            println("Total inner while time: ", inner_while_total)
            println("Total matzeros_total_time: ", matzeros_total_time)
            println("Total create_ind_vect_total_time: ", create_ind_vect_total_time)
            println("Total find_min_total_time: ", find_min_total_time)
            println("Total sub_min_total_time: ", sub_min_total_time)
            println("Total add_min_total_time: ", add_min_total_time)
            println("c_outer_while: ", c_outer_while)
            println("c_inner_while: ", c_inner_while)
            return matching.match
        else
            sorted_matching = sort(matching.match)[1:matching.cardinality]

            # mark all rows which aren't matched
            mark_start = time()
            
            rows_marked[:] .= false
            cols_marked[:] .= false
            marked_col_ind = Vector{Int64}()
            marked_row_ind = Vector{Int64}()
            new_marked_row_ind = Vector{Int64}()
            last = 0
            for r=1:matching.cardinality
                if sorted_matching[r]-last != 1
                    for ri=last+1:sorted_matching[r]-1
                        rows_marked[ri] = true
                        push!(marked_row_ind,ri)
                        push!(new_marked_row_ind,ri)
                    end
                    last = sorted_matching[r]
                else
                    last = sorted_matching[r]
                end
            end
            for r=last+1:nrows
                rows_marked[r] = true
                push!(marked_row_ind,r)
                push!(new_marked_row_ind,r)
            end
        end
       
        start_inner = time()
        changed = true
        while changed
            c_inner_while += 1
            
            changed = false
            new_marked_col_ind = Vector{Int64}()
            # mark cols
            in_start = time()
            idx_not_marked_cols = findall(x->!x,cols_marked)
            @inbounds for c in idx_not_marked_cols
                @inbounds for r in new_marked_row_ind
                    @inbounds if mat[r,c] == 0
                        cols_marked[c] = true
                        push!(marked_col_ind,c)
                        push!(new_marked_col_ind,c)
                        break
                    end
                end
            end

            new_marked_row_ind = Vector{Int64}()
            idx_not_marked_rows = findall(x->!x,rows_marked)
            # mark new rows
            @inbounds for r in idx_not_marked_rows
                @inbounds for c in new_marked_col_ind
                    @inbounds if matching.match[c] == r
                        rows_marked[r] = true
                        push!(marked_row_ind,r)
                        push!(new_marked_row_ind,r)
                        changed = true
                        break
                    end
                end
            end
        end
        inner_while_total += time()-start_inner

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
            end
        end
        sub_min_total_time += time()-start_timer

        start_timer = time()        
        # add minimum where !rows_marked but cols_marked
        @inbounds @views mat[ind_unmarked_rows_sub,ind_marked_cols_sub] .+= min_val
        add_min_total_time += time()-start_timer
    end
end

function get_matching(mat)
    # create bi partite matching graph
    zero_vals = findall(x->x==0,mat)
    ei = Vector{Int64}()
    ej = Vector{Int64}()

    for yx in zero_vals
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