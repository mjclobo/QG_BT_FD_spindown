

function arakawa_jacobian_doubly_periodic(a, b, dx, dy)
    Nx, Ny = size(a)
    J = zeros(Float64, Nx, Ny)

    dx2 = 2dx
    dy2 = 2dy
    denom = 4dx * dy

    # Periodic indices in x
    ip = [i == Nx ? 1 : i + 1 for i in 1:Nx]
    im = [i == 1  ? Nx : i - 1 for i in 1:Nx]

    # Periodic indices in y
    jp = [j == Ny ? 1 : j + 1 for j in 1:Ny]
    jm = [j == 1  ? Ny : j - 1 for j in 1:Ny]

    @threads for i in 1:Nx
        i_p = ip[i]
        i_m = im[i]

        @inbounds for j in 1:Ny
            j_p = jp[j]
            j_m = jm[j]

            dady = (a[i, j_p] - a[i, j_m]) / dy2
            dbdx = (b[i_p, j] - b[i_m, j]) / dx2

            dadx = (a[i_p, j] - a[i_m, j]) / dx2
            dbdy = (b[i, j_p] - b[i, j_m]) / dy2

            J1 = dadx * dbdy - dady * dbdx

            J2 = (
                a[i_p, j] * (b[i_p, j_p] - b[i_p, j_m]) -
                a[i_m, j] * (b[i_m, j_p] - b[i_m, j_m]) -
                a[i, j_p] * (b[i_p, j_p] - b[i_m, j_p]) +
                a[i, j_m] * (b[i_p, j_m] - b[i_m, j_m])
            ) / denom

            J3 = (
                a[i_p, j_p] * b[i, j_p] - a[i_m, j_p] * b[i, j_p] -
                a[i_p, j_m] * b[i, j_m] + a[i_m, j_m] * b[i, j_m] -
                a[i_p, j_p] * b[i_p, j] + a[i_p, j_m] * b[i_p, j] +
                a[i_m, j_p] * b[i_m, j] - a[i_m, j_m] * b[i_m, j]
            ) / denom

            J[i, j] = (J1 + J2 + J3) / 3
        end
    end

    return J
end

function BT_jacobian(ah, b, kx, ky)
    # J(f,g)=∂y​[(∂x​f)g]−∂x​[(∂y​f)g]

    ∂xa = real.(ifft(im .* kx .* ah))
    ∂ya = real.(ifft(im .* ky .* ah))

    return im .* ky .* fft(∂xa .* b) .- im .* kx .* fft(∂ya .* b)
end

"""
    conservative_derivs(f, dx, dy)

Compute conservative first and second derivatives of 2D array `f` on a uniform grid.

Returns:
    fx  - ∂f/∂x (conservative, periodic)
    fy  - ∂f/∂y (conservative, periodic)
    lapf - ∇² f = ∂²f/∂x² + ∂²f/∂y²
"""
function ddx_ddy_reentrant_FD(f::Array{Float64,2}, dx::Float64, dy::Float64)
    Nx, Ny = size(f)
    
    # First derivatives (central differences, flux-conservative)
    fx = zeros(Nx, Ny)
    fy = zeros(Nx, Ny)
    
    # ∂/∂x
    for j in 1:Ny
        for i in 1:Nx
            ip = i < Nx ? i+1 : 1   # periodic
            im = i > 1  ? i-1 : Nx
            fx[i,j] = (f[ip,j] - f[im,j]) / (2*dx)
        end
    end
    
    # ∂/∂y
    for i in 1:Nx
        for j in 1:Ny
            jp = j < Ny ? j+1 : 1
            jm = j > 1  ? j-1 : Ny
            fy[i,j] = (f[i,jp] - f[i,jm]) / (2*dy)
        end
    end
    
    # # Second derivatives (Laplacian)
    # lapf = zeros(Nx, Ny)
    # for j in 1:Ny
    #     jp = j < Ny ? j+1 : 1
    #     jm = j > 1  ? j-1 : Ny
    #     for i in 1:Nx
    #         ip = i < Nx ? i+1 : 1
    #         im = i > 1  ? i-1 : Nx
    #         lapf[i,j] = (f[ip,j] - 2f[i,j] + f[im,j]) / dx^2 +
    #                     (f[i,jp] - 2f[i,j] + f[i,jm]) / dy^2
    #     end
    # end
    
    return dfdx, dfdy
end


function wavemaker_m1(ψ0, x, y, x0, y0, R, δ, t, τ)
    # x: 1D array of length Nx
    # y: 1D array of length Ny
    # returns Nx x Ny array

    Nx = length(x)
    Ny = length(y)

    # Make 2D meshgrid arrays
    X = repeat(x, 1, Ny)          # Nx x Ny, each row is x
    Y = repeat(y', Nx, 1)         # Nx x Ny, each column is y

    # Shift coordinates to patch center
    Xc = X .- x0
    Yc = Y .- y0

    # Radius and azimuth
    r = @. sqrt(Xc^2 + Yc^2)
    θ = @. atan(Yc, Xc)

    # m=1 azimuthal wavemaker at patch edge
    ψ_f = @. ψ0 * sin(2*pi*t/τ) * r * cos(θ) * exp(-((r - R)^2 / (2*δ^2)))

    return ψ_f
end


function define_reentrant_lap_op_FD(Nx, Ny, dx, dy; inv_op=true)
    N = Nx*Ny
    
    dx2 = dx^2
    dy2 = dy^2

    # function to convert (i,j) to linear index
    idx(i,j) = (j-1)*Nx + i

    # build sparse matrix
    rows = Int[]
    cols = Int[]
    vals = Float64[]

    for j in 1:Ny
        jp = j < Ny ? j+1 : 1
        jm = j > 1 ? j-1 : Ny
        for i in 1:Nx
            ip = i < Nx ? i+1 : 1
            im = i > 1 ? i-1 : Nx

            center = idx(i,j)

            push!(rows, center); push!(cols, center); push!(vals, -2/dx2 - 2/dy2 - 1/Ld^2)
            push!(rows, center); push!(cols, idx(ip,j)); push!(vals, 1/dx2)
            push!(rows, center); push!(cols, idx(im,j)); push!(vals, 1/dx2)
            push!(rows, center); push!(cols, idx(i,jp)); push!(vals, 1/dy2)
            push!(rows, center); push!(cols, idx(i,jm)); push!(vals, 1/dy2)

        end
    end

    L = sparse(rows, cols, vals, N, N)
    if inv_op==true
        L[1,:] .= 1.0
    end

    return L
end


function define_laplacian_FD(Nx, Ny, dx, dy)
    dx2 = dx^2
    dy2 = dy^2

    rows = Int[]
    cols = Int[]
    vals = Float64[]

    idx(i,j) = (j-1)*Nx + i

    for j in 1:Ny
        jp = (j == Ny) ? 1 : j+1
        jm = (j == 1)  ? Ny : j-1
        for i in 1:Nx
            ip = (i == Nx) ? 1 : i+1
            im = (i == 1)  ? Nx : i-1

            center = idx(i,j)

            push!(rows, center); push!(cols, center); push!(vals, -2/dx2 -2/dy2)
            push!(rows, center); push!(cols, idx(ip,j)); push!(vals, 1/dx2)
            push!(rows, center); push!(cols, idx(im,j)); push!(vals, 1/dx2)
            push!(rows, center); push!(cols, idx(i,jp)); push!(vals, 1/dy2)
            push!(rows, center); push!(cols, idx(i,jm)); push!(vals, 1/dy2)
        end
    end

    return sparse(rows, cols, vals, Nx*Ny, Nx*Ny)
end


# define ∇^2 - Ld^-2 operator matrix
L_re_FD_inv = define_reentrant_lap_op_FD(Nx, Ny, dx, dy)
L_re_FD     = define_reentrant_lap_op_FD(Nx, Ny, dx, dy; inv_op=false)
Lap_op      = define_laplacian_FD(Nx, Ny, dx, dy)


function rhs_BT_FD(q, t)

    Q = q .- f0
    Q .-= mean(Q)           # ensure zero mean
    Q_vec = reshape(Q, Nx*Ny)
    ψ_v = reshape(L_re_FD_inv \ Q_vec, Nx, Ny)

    ψ_f = wavemaker(ψ0, x, y, x0, y0, δx, δy, t, τ)
    # ψ_f .-= mean(ψ_f)
    ψ_f_vec = reshape(ψ_f, Nx*Ny)

    q_f = -reshape(L_re_FD * ψ_f_vec, Nx, Ny)

    dqdt = - arakawa_jacobian_doubly_periodic(ψ_v .+ ψ_f, q .+ q_f, dx, dy)

    dqdt .-= r .* q

    return dqdt

end

function rhs(qh, t)

    q = real.(ifft(qh))
    Qh = fft(q .- f0)
    ψ_vh = @. -(k2D + Ld^-2).^-1 * Qh

    ψ_f = wavemaker(ψ0, x, y, x0, y0, δx, δy, t, τ)

    # ψ_f = wavemaker_m1(ψ0, x, y, x0, y0, radius, δx, t, τ)

    dqhdt = - BT_jacobian(ψ_vh .+ fft(ψ_f), q .+ real.(ifft((k2D .+ Ld^-2) .* fft(ψ_f))), kx, ky)

    dqhdt .-= ν .* k2D .* k2D .* qh
    
    dqhdt .-= r .* qh

    return dqhdt
end


function u_from_psi(ψ)
    u = -d_dy(ψ, dy)

    ψ_hat = rfft(ψ, 1)                # (Nx/2+1, Ny)

    v_hat = 1im .* KXr .* ψ_hat       # dψ/dx in spectral space
    v = irfft(v_hat, Nx, 1) |> real   # back to x space

    return u, v
end


function u_from_psi_fft(ψh, kx, ky)
    return -real.(ifft(im .* ky .* ψh)), real.(ifft(im .* kx .* ψh))
end


function init_solver(qh0, kx, ky, k2D, Ld, ψ0, x, y, x0, y0, δx, δy)

    Nx, Ny = size(qh0)

    # -------- FFT plans --------
    fft_plan  = plan_fft(qh0;  flags=FFTW.MEASURE)
    ifft_plan = plan_ifft(qh0; flags=FFTW.MEASURE)

    # -------- Spectral-space complex arrays --------
    q       = similar(qh0)             # will store ifft(qh)
    Qh      = similar(qh0)
    ψ_vh    = similar(qh0)
    tmph    = similar(qh0)
    tmp_c   = similar(qh0)             # complex buffer for FFT inputs
    tmp2_c  = similar(qh0)             # complex buffer for Jacobian FFTs
    J       = similar(qh0)

    # -------- Physical-space real arrays --------
    ψ_f     = similar(qh0, Float64)
    tmp_r   = similar(qh0, Float64)
    tmp2_r  = similar(qh0, Float64)

    # -------- Precompute Gaussian wavemaker spatial shape --------
    G = @. exp( -((x - x0)^2 / (2δx^2)) - ((y - y0)^2 / (2δy^2)) )

    # -------- Precompute K = k2D + Ld^-2 --------
    K  = @. k2D + Ld^-2
    ikx = im .* kx
    iky = im .* ky

    return (
        fft_plan  = fft_plan,
        ifft_plan = ifft_plan,

        # complex fields
        q       = q,
        Qh      = Qh,
        ψ_vh    = ψ_vh,
        tmph    = tmph,
        tmp_c   = tmp_c,
        tmp2_c  = tmp2_c,
        J       = J,

        # real fields
        ψ_f     = ψ_f,
        tmp_r   = tmp_r,
        tmp2_r  = tmp2_r,

        # constants
        G       = G,
        K       = K,
        ikx     = ikx,
        iky     = iky,
        ψ0      = ψ0
    )
end

# ================================================================
# ==================== WAVEMAKER =================================
# ================================================================
function wavemaker!(ψ_f, ψ0, G, t, τ)
    s = sin(2π * t / τ)
    @. ψ_f = ψ0 * G * s
end

# ================================================================
# ==================== JACOBIAN ==================================
# ================================================================
function BT_jacobian!(
    J, ah, b,
    ikx, iky,
    fft_plan, ifft_plan,
    tmp_r, tmp2_r, tmph, tmp2_c
)
    # ∂x a
    @. tmph = ikx * ah
    mul!(tmp2_c, ifft_plan, tmph)
    @. tmp_r = real(tmp2_c)

    # FFT(∂x a * b)
    @. tmp2_r = tmp_r * b
    @. tmp2_c = complex(tmp2_r)
    mul!(tmph, fft_plan, tmp2_c)
    @. tmph = iky * tmph

    # ∂y a
    @. tmph = iky * ah
    mul!(tmp2_c, ifft_plan, tmph)
    @. tmp_r = real(tmp2_c)

    # FFT(∂y a * b)
    @. tmp2_r = tmp_r * b
    @. tmp2_c = complex(tmp2_r)      # keep tmp2_c unchanged
    mul!(tmph, fft_plan, tmp2_c)     # tmph = FFT(tmp2_c)
    @. tmph = ikx * tmph             # multiply by ikx

    # Combine
    @. J = tmph - tmp2_c
end

# ================================================================
# ==================== RHS FUNCTION ==============================
# ================================================================
# SPECTRAL 

function rhs!(
    dqhdt, qh, t,
    ν, r, τ, f0, kx, ky,
    params
)
    @unpack fft_plan, ifft_plan,
            q, Qh, ψ_vh, ψ_f,
            tmph, tmp_c, tmp2_c, J,
            tmp_r, tmp2_r,
            G, K, ikx, iky, ψ0 = params

    # ------------------------------------------------
    # 1. q = real(ifft(qh))
    # ------------------------------------------------
    mul!(q, ifft_plan, qh)
    @. q = real(q)

    # ------------------------------------------------
    # 2. Qh = fft(q - f0)
    # ------------------------------------------------
    @. tmp_r = q - f0
    @. tmp_c = complex(tmp_r)
    mul!(Qh, fft_plan, tmp_c)

    # ------------------------------------------------
    # 3. ψ_vh = -Qh / K
    # ------------------------------------------------
    @. ψ_vh = -Qh / K

    # ------------------------------------------------
    # 4. Wavemaker forcing ψ_f (physical space)
    # ------------------------------------------------
    wavemaker!(ψ_f, ψ0, G, t, τ)  # writes directly into ψ_f

    # ------------------------------------------------
    # 5. Spectral field for Jacobian first argument: ψ_vh + FFT(ψ_f)
    # ------------------------------------------------
    @. tmp_c = ψ_f
    mul!(tmph, fft_plan, tmp_c)    # tmph = FFT(ψ_f)
    @. tmph += ψ_vh                # spectral field

    # ------------------------------------------------
    # 6. Physical-space field for Jacobian second argument: q + ifft(K * FFT(ψ_f))
    # ------------------------------------------------
    @. tmp_c = ψ_f
    mul!(tmp2_c, fft_plan, tmp_c)  # tmp2_c = FFT(ψ_f)
    @. tmp2_c = K .* tmp2_c         # multiply by K
    mul!(tmp_c, ifft_plan, tmp2_c) # tmp_c = ifft(K * FFT(ψ_f))
    @. tmp2_r = q + real.(tmp_c)   # second argument b

    # ------------------------------------------------
    # 7. Nonlinear Jacobian
    # ------------------------------------------------
    BT_jacobian!(J, tmph, tmp2_r,
                 ikx, iky,
                 fft_plan, ifft_plan,
                 tmp_r, tmp2_r, tmph, tmp2_c)

    @. dqhdt = -J

    # ------------------------------------------------
    # 8. Linear terms
    # ------------------------------------------------
    # viscosity (optional, precompute k4D if needed)
    # @. dqhdt -= ν * k4D * qh

    # linear drag
    @. dqhdt -= r * qh

    return nothing
end


######################################################################
## Below this line are functions for an optimized and
## GPU-compatible version of the FD BT model
######################################################################

struct BTParamsOpt
    Nx::Int
    Ny::Int
    dx::Float64
    dy::Float64
    r::Float64
    f0::Float64
    L_fac          # pre-factorized Laplacian operator
    L_op::SparseMatrixCSC{Float64,Int}
    Lap_op
    ip::Vector{Int}
    im::Vector{Int}
    jp::Vector{Int}
    jm::Vector{Int}
    backend     # KA backend, e.g. CPU() or CUDABackend()
end

@kernel function arakawa_kernel!(
    J, a, b, ip, im, jp, jm, dx, dy
)
    Nx, Ny = size(a)
    dx2 = 2dx
    dy2 = 2dy
    denom = 4dx * dy

    I = @index(Global, Linear)
    total = Nx * Ny

    # Guard without return
    if I <= total
        # convert linear index to (i,j)
        j = (I - 1) ÷ Nx + 1
        i = I - (j - 1) * Nx

        i_p = ip[i];  i_m = im[i]
        j_p = jp[j];  j_m = jm[j]

        #### Central derivatives
        dady = (a[i, j_p] - a[i, j_m]) / dy2
        dbdx = (b[i_p, j] - b[i_m, j]) / dx2

        dadx = (a[i_p, j] - a[i_m, j]) / dx2
        dbdy = (b[i, j_p] - b[i, j_m]) / dy2

        J1 = dadx * dbdy - dady * dbdx

        #### Arakawa terms
        J2 = (
            a[i_p, j] * (b[i_p, j_p] - b[i_p, j_m]) -
            a[i_m, j] * (b[i_m, j_p] - b[i_m, j_m]) -
            a[i, j_p] * (b[i_p, j_p] - b[i_m, j_p]) +
            a[i, j_m] * (b[i_p, j_m] - b[i_m, j_m])
        ) / denom

        J3 = (
            a[i_p, j_p] * b[i, j_p] - a[i_m, j_p] * b[i, j_p] -
            a[i_p, j_m] * b[i, j_m] + a[i_m, j_m] * b[i, j_m] -
            a[i_p, j_p] * b[i_p, j] + a[i_p, j_m] * b[i_p, j] +
            a[i_m, j_p] * b[i_m, j] - a[i_m, j_m] * b[i_m, j]
        ) / denom

        J[i, j] = -(J1 + J2 + J3) / 3   # including minus sign because we want negative of Jaocobian
    end
end


function rhs_BT_FD!(
    dqdt, q,
    ψ_v, ψ_v_vec,
    Q, Q_vec,
    ψ_f, ψ_f_vec, q_f, ψ_sum, q_sum,
    lap_q_vec, lap2_q_vec,
    params::BTParamsOpt, t)

    @unpack Nx, Ny, dx, dy, r, f0, L_fac, L_op, ip, im, jp, jm, backend = params

    #### 1. Build potential vorticity anomaly Q = q - f0, zero mean
    @tturbo for j in 1:Ny, i in 1:Nx
        Q[i,j] = q[i,j] - f0
    end
    Q .-= sum(Q) / (Nx*Ny)

    #### 2. Invert Laplacian to get ψ_v
    @views Q_vec .= reshape(Q, :)
    ψ_v_vec .= L_fac \ Q_vec
    @views ψ_v .= reshape(ψ_v_vec, Nx, Ny)

    #### 3. Forcing streamfunction
    ψ_f .= wavemaker(ψ0, x, y, x0, y0, δx, δy, t, τ)

    #### 4. q_f = -L ψ_f
    @views ψ_f_vec .= reshape(ψ_f, :)
    q_f_vec = -(L_op * ψ_f_vec)
    @views q_f .= reshape(q_f_vec, Nx, Ny)

    #### 5. Arakawa Jacobian J(ψ_v + ψ_f, q + q_f)
    # Inputs must be pre-summed
    ψ_sum .= ψ_v .+ ψ_f
    q_sum .= q    .+ q_f    

    k = arakawa_kernel!(params.backend)

    event = k(
        dqdt,
        ψ_sum,      # (already contains ψ_v + ψ_f)
        q_sum,        # (already contains q + q_f)
        params.ip, params.im, params.jp, params.jm,
        params.dx, params.dy;
        ndrange = params.Nx * params.Ny
    )

    KernelAbstractions.synchronize(params.backend)

    #### 6. Linear drag
    @tturbo for j in 1:Ny, i in 1:Nx
        dqdt[i,j] -= r * q[i,j]
    end

    # ---- hyperviscosity (4th order) ----

    # 1. Compute Laplacian(q)
    lap_q_vec .= params.Lap_op * reshape(q, :)

    # 2. Compute Laplacian(Laplacian(q)) = ∇⁴ q
    lap2_q_vec .= params.Lap_op * lap_q_vec

    # 3. Add hyperviscous tendency
    dqdt_vec = reshape(dqdt, :)
    dqdt_vec .-= ν * lap2_q_vec
    dqdt .= reshape(dqdt_vec, Nx, Ny)

    return dqdt
end


