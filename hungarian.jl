using MatrixNetworks

function run_hungarian(mat_original) 
    mat = copy(mat_original)    

    # subtract col minimum
    for c=1:size(mat)[2]
        mat[:,c] -= minimum(mat[:,c])
    end

    # subtract row minimum
    for r=1:size(mat)[1]
        mat[r,:] -= minimum(mat[r,:])
    end

    nrows = size(mat)[1]
    ncols = size(mat)[2]
    rows_marked = zeros(Bool,nrows)
    cols_marked = zeros(Bool,ncols)

    mask = ones(Bool, nrows, ncols)

    col_unmask = ones(Bool, ncols)
    row_mask = ones(Bool, nrows)

    min_total = 0.0
    while true
        matching = get_matching(mat)

        if matching.cardinality == ncols
            println("Total min time: ", min_total)
            return matching.match
        else
            sorted_matching = sort(matching.match)[1:matching.cardinality]

            # mark all rows which aren't matched
            mark_start = time()
            
            rows_marked[:] = false
            cols_marked[:] = false
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
       
        changed = true
        while changed
            changed = false
            new_marked_col_ind = Vector{Int64}()
            # mark cols
            in_start = time()
            idx_not_marked_cols = find(x->!x,cols_marked)
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
            idx_not_marked_rows = find(x->!x,rows_marked)
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

        start_min = time()

        mask[:,:] = true
        mask[:,marked_col_ind] = false
        row_mask[:] = true
        row_mask[marked_row_ind] = false
        mask[row_mask,:] = false
    
        min_val = minimum(mat[mask])
        mat .-= min_val*mask

        mask[:,:] = true
        mask[marked_row_ind,:] = false
        col_unmask[:] = true
        col_unmask[marked_col_ind] = false
        mask[:,col_unmask] = false
        mat .+= min_val*mask

        min_total += time() - start_min
    end
end

function get_matching(mat)
    # create bi partite matching graph
    zero_vals = find(x->x==0,mat)
    ei = Vector{Int64}()
    ej = Vector{Int64}()

    for i in zero_vals
        # zero based
        zi = i-1
        zr = zi % size(mat)[1]
        push!(ej,zr+1) # row
        push!(ei,div((zi-zr),size(mat)[1])+1) # col
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