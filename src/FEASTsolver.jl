module FEASTsolver

using LinearAlgebra, FastGaussQuadrature

export feast, cauchynodes

"""
    feast(A, M0, λmin, λmax; cpts=8, tol=1.0e-12, maxiter=20)
    feast(A, Y, zs, ws, λmin, λmax; tol=1.0e-12, maxiter=20)

Determines all the eigenvalues and vectors of an Hermitian matrix `A` inside the interval `]λmin, λmax]`. 

`tol` is the relative error tolerance. The error is estimated using the variation of the trace of the eigenvalues in `]λmin, λmax]`.
`maxiter` is the number of FEAST outer iterations.

`M0` is the size of the subspace search. It must be a bit larger than the number of eigenvalues in the interval. 
It is also possible to provide an initial matrix `Y` as the initial search basis. In that case `M0=size(Y, 2)`.

`cpts` is the number of Gauss-Legendre points using in the approximation of the Cauchy contour integral. 
The nodes `zs` and weights `ws` for the quadrature rule can also be passed explicitly.
"""
function feast(A, M0::Integer, λmin::Real, λmax::Real; cpts::Integer=8, tol::Real=1.0e-12, maxiter::Integer=20)
    # initial search space of size M0
    V0 = rand(complex(eltype(A)), size(A, 1), M0)
    # contour integral nodes and and weights
    zs, ws = cauchynodes(λmin, λmax, cpts)
    
    return feast(A, V0, zs, ws, λmin, λmax; tol=tol, maxiter=maxiter)
end

function feast(A, Y, zs, ws, λmin::Real, λmax::Real; tol::Real=1.0e-12, maxiter::Integer=20)
    
    # dimension of problem and of subspace search
    N, M0 = size(Y)
    
    # allocate output
    ritzvals = Vector{Complex{eltype(A)}}(undef, M0)
    ritzvecs = Matrix{Complex{eltype(A)}}(undef, M0, M0)
    
    # allocate computation arrays
    Q = similar(Y)
    q = Vector{eltype(Y)}(undef, N)
    AQ = Matrix{eltype(Y)}(undef, M0, M0)
    temp = similar(Y)
    BQ = Matrix{eltype(Y)}(undef, M0, M0)
    
    
    # convercenve history
    neighistory = zeros(Int, maxiter)
    errorhistory = zeros(Float64, maxiter)
    
    # error estimation
    previous = 0.0
    rerror = 0.0    
    iter = 1
    
    # evaluate and store necessary factorizations
    F = [lu(z*I - A) for z in zs]
    
    while true
        
        # apply rational filter to Y basis
        fill!(Q, zero(eltype(Q)))
        
        for i in 1:length(zs)
            for n in axes(Y, 2)
    
                ldiv!(q, F[i], view(Y, :, n))
                
                Q[:, n] .+= ws[i]*q

            end
        end
        
        # form effective matrices in reduced subspace
        mul3mat!(AQ, Q, A, Q, temp)
        mul!(BQ, Q', Q)
        
        # diagonalize subspace problem
        ritzvals, ritzvecs = eigen!(Hermitian(AQ), Hermitian(BQ), sortby = λ -> (!(λmin*(1+tol) < real(λ) <= λmax*(1+tol)), real(λ), imag(λ)))
        
        # update eigenvectors
        mul!(Y, Q, ritzvecs)
        
        # determine number of eigenvalues is inside
        neig = 0
        for λ in ritzvals
            if λmin < real(λ) <= λmax
                neig += 1
            end
        end
        neighistory[iter] = neig
        
        # estimate relative error of the trace of the eigenvalues
        tr = sum(ritzvals[1:neig])
        rerror = abs(tr - previous)/abs(tr)
        previous = sum(ritzvals[1:neig])
        errorhistory[iter] = rerror
        
        if rerror >= tol 
            if iter < maxiter # if tolerance is not satisfied and number of iterations small than max, go to next iteration
                iter += 1
            else # if tolerance is not satisfied but max number of iterations is exceded, throw warning and return result
                @warn "maximum number of iterations, $maxiter, reached"
                return ritzvals, Y, neighistory, errorhistory
            end
        else # if tolerance is satisfied, return result
            return ritzvals, Y, neighistory, errorhistory
        end
    
        end # end while

end

function mul3mat!(Out, A, B, C, temp)
    mul!(temp, B, C)
    mul!(Out, A', temp)
    return Out
end

"""
    cauchynodes(λmin, λmax, cpts)

Builds the nodes and weights for a `cpts`-points Gauss-Legendre quadrature rule along a circle in the complex plane which includes `]λmin, λmax]`.
"""
function cauchynodes(λmin, λmax, cpts)
    
    xs, ws = gausslegendre(cpts)
    
    c = (λmax + λmin)/2
    R = (λmax - λmin)/2
    
    weights = zeros(Complex{Float64}, 2*cpts)
    points = zeros(Complex{Float64}, 2*cpts)
    
    for k in 1:cpts
        weights[k] = ws[k]*R*cis(pi*(xs[k]+1)/2)/4
        points[k] = c + R*cis(pi*(xs[k]+1)/2)
        
        weights[k + cpts] = ws[k]*R*cis(-pi*(xs[k]+1)/2)/4
        points[k + cpts] = c + R*cis(-pi*(xs[k]+1)/2)
    end
    
    return points, weights
end

end
