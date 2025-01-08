function MatrixEquations.lyapc(S::LinearAlgebra.Schur, C::AbstractMatrix; adj=false)
    X = MatrixEquations.utqu(C,S.Z)
    MatrixEquations.lyapcs!(S.T, X; adj=adj)
    MatrixEquations.utqu!(X,S.Z')
    return X
end

function MatrixEquations.lyapc(A::Diagonal, C::AbstractMatrix, F::AbstractMatrix = diag(A) .+ diag(A)'; adj=false)
    return adj ? -C ./ conj(F) : -C ./ F 
end
