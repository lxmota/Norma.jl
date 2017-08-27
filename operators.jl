#
#
function odot(A, B)
    m = size(A, 1)
    n = size(A, 2)
    if m ≠ n
        error("first argument must be a square matrix")
    end
    p = size(B, 1)
    q = size(B, 2)
    if p ≠ q
        error("second argument must be a square matrix")
    end
    if m ≠ p
        error("arguments must be of the same size")
    end
    C = zeros(m, m, m, m)
    for a = 1 : m
        for b = 1 : m
            for c = 1 : m
                for d = 1 : m
                    C[a, b, c, d] = A[a, c] * B[b, d] + A[a, d] * B[b, c]
                end
            end
        end
    end
    C = 0.5 * C
    return C
end

#
#
function ox(A, B)
    m = size(A, 1)
    n = size(A, 2)
    if m ≠ n
        error("first argument must be a square matrix")
    end
    p = size(B, 1)
    q = size(B, 2)
    if p ≠ q
        error("second argument must be a square matrix")
    end
    if m ≠ p
        error("arguments must be of the same size")
    end
    C = zeros(m, m, m, m)
    for i = 1 : m
        for j = 1 : m
            C[i, j, :, :] = A[i, j] * B; 
        end
    end
    return C
end

#
#
function convecttangent(CC, S, F)
    n = size(F, 1)
    I = eye(n)
    AA = zeros(n, n, n, n)
    for i = 1 : n
        for j = 1 : n
            for k = 1 : n
                for l = 1 : n
                    s = 0.0
                    for p = 1 : n
                        for q = 1 : n
                            s = s + F[i, p] * CC[p, j, l, q] * F[k, q]
                        end
                    end
                    AA[i, j, k, l] = S[l, j] * I[i, k] + s
                end
            end
        end
    end
    return AA
end

#
#
function secondfromfourth(AA)
    dim = size(AA, 1)
    A = zeros(dim * dim)
  
    for i = 1 : dim
        for j = 1 : dim
            p = dim * (i - 1) + j
            for k = 1 : dim
                for l = 1 : dim
                    q = dim * (k - 1) + l
                    A[p, q] = AA[i, j, k, l]
                end
            end
        end
    end
    return A
end
