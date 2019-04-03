using MatrixNetworks

function subMins!(mat)
    mcols = minimum(mat, dims=1)
    for j = 1:size(mat,2), i = 1:size(mat,1)
        mat[i,j] -= mcols[j] 
    end  
    mrows = minimum(mat, dims=2)
    for j = 1:size(mat,2), i = 1:size(mat,1)
        mat[i,j] -= mrows[i]
    end 
end 

function run_hungarian!(mat) 
    # subtract col minimum and then row minimum
    subMins!(mat)

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

    matzeros = findall(iszero, mat)
    nmatzeros = length(matzeros)
    max_nmatzeros = nmatzeros

    addition_vector = zeros(Int64, nrows)

    while true
        @views matzeros_sub = matzeros[1:nmatzeros]
        new_marked_col_ind[:] .= false
        new_marked_row_ind[:] .= false

       
        matching  = get_matching(matzeros_sub)

        if matching.cardinality == ncols
            return matching.match
        else
            sorted_matching = partialsort(matching.match, 1:matching.cardinality)

            # mark all rows which aren't matched
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
       
        changed = true
        while changed
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
        
        #minimum where !cols_marked but rows_marked
        @views mat_ind_unmarked_cols_sub = mat[:,ind_unmarked_cols_sub]
        @views rmins = vec(minimum(mat_ind_unmarked_cols_sub[ind_marked_rows_sub,:], dims=2))
        @views min_val = minimum(rmins)
        
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
       
        # add minimum where !rows_marked but cols_marked
        addition_vector .= 0
        addition_vector[ind_unmarked_rows_sub] .= min_val
        @views mat_ind_marked_cols_sub = mat[:,ind_marked_cols_sub]
        for col in eachcol(mat_ind_marked_cols_sub)
            col .+= addition_vector
        end
    end
end

function get_matching(matzeros_sub)
    # create bi partite matching graph
    nmatzeros = length(matzeros_sub)
    ei = zeros(Int64, nmatzeros)
    ej = zeros(Int64, nmatzeros)

    i = 1
    for rc in matzeros_sub
        ej[i] = rc[1] # row
        ei[i] = rc[2] # col
        i += 1
    end
    
    matching = bipartite_matching(ones(Int64,nmatzeros), ei, ej)
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