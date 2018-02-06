using Hungarian
include("hungarian.jl")

function benching()
    seed = srand(7)

    for n in [10, 50, 100, 200, 400, 800, 1000, 2000, 4000, 8000]
        gc()
        
        A = rand(seed, UInt16, n, n)
        A[A .> 50000] = 50000
        B = copy(A)
        println("n: ", n)
        println("Lib: ")
        @time matching_lib =  Hungarian.munkres(A)
        matching_lib = [findfirst(matching_lib[:,i].==Hungarian.STAR) for i = 1:n]
        println("Own: ")
        @time matching = run_hungarian(B)
    
        sum_lib = 0.0
        c = 1
        for r in matching_lib
            sum_lib += B[r,c]
            c += 1
        end
        sum_own = 0.0
        c = 1
        for r in matching
            sum_own += B[r,c]
            c += 1
        end
        println("Sum lib: ", sum_lib)
        println("Sum own: ", sum_own)
        @assert sum_lib == sum_own
    end
end

benching()
