# Hungarian Method

The code is explained on my blog [opensourc.es/blog/hungarian-method](http://opensourc.es/blog/hungarian-method) and the updated performance version on
[opensourc.es/blog/hungarian-performance](http://opensourc.es/blog/hungarian-performance)


If you want to run the code simply include `hungarian.jl` and run `run_hungarian(A)` where a is a `n x n` matrix. The benchmarking can be done by running `benchmark.jl`.

There you need to install [Hungarian](https://github.com/Gnimuc/Hungarian.jl) to compare this code against that library. If you find some bug or performance improvement please make a PR, issue and/or comment on my blog. Thanks in advance!

## Attention:

The method works best if you specify the type of the matrix as good as possible. For example using UInt16 instead of Float64.
Int in general is much much better than Float here. The problem is that in the code sometimes small numbers get added to some part of the matrix. This can result in an overflow and therefore produce wrong results. I'm **not** checking for overflows!
