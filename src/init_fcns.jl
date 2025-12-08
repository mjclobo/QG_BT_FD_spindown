################################################################################
# Initialize
################################################################################

struct BTParams
    Nx::Int
    Ny::Int
    nt::Int
    Lx::Float64
    Ly::Float64
    dt::Float64
    beta::Float64
    Ld::Float64
    ν::Float64
    r::Float64
    U0::Float64
end

################################################################################
# Functions for restart and building psi of time
################################################################################


function extract_time(filename; file_strings=false)
    # Find the starting index of the search string
    idx_start = findfirst("_t", filename)[end] + 1
    idx_end   = findfirst(".jld", filename)[1] - 1

    time = filename[idx_start:idx_end]
    if file_strings==false
        return parse(Float64, time)
    else
        return filename[1:idx_start-1], filename[idx_end+1:end]
    end
end

function define_t_of_saved_files(Ny, save_path)
    ψfiles = readdir(save_path)
    t_array = zeros(length(ψfiles))

    for (i, file) in enumerate(ψfiles)
        t_array[i] = extract_time(save_path*file)
    end

    return sort(t_array)
end

function construct_psi_of_t(t_array, Ny_saved, save_path)
    ψ1_of_t = zeros(Nx, Ny_saved, length(t_array))
    ψ2_of_t = zeros(Nx, Ny_saved, length(t_array))
    file_start, file_end = extract_time(readdir(save_path)[1]; file_strings=true)

    for (i, t) in enumerate(t_array)
        savedata = load(save_path * file_start * string(t) * file_end)

        ψ1_of_t[:,:,i] = savedata["jld_data"]["ψ1"]
        ψ2_of_t[:,:,i] = savedata["jld_data"]["ψ2"]
    end
    return ψ1_of_t, ψ2_of_t
end

################################################################################
# Basic time stepping loop
################################################################################

function run_BT_model(qh, t, params; timestepper="RK4", output_every=500)
    start_time = time()
    # To turn off the PyPlot GUI
    PyPlot.pygui(false)

    # calculating y indices for prescribed meridional width of save domain
    save_ind_start = floor(Int, Ny * y_width / 2)
    save_ind_end   = floor(Int, Ny * (1 - y_width / 2))

    # defining a couple of counters
    cnt=1
    ell=1

    dqhdt = similar(qh)

    for n = 1:nt

        # rk4_step!(qh,dt,rhs!,ν, r, τ, f0,kx, ky,params)

        t+= dt

        qh = rk4_BT(qh, t, dt)

        if mod(n, output_every) == 0      # output a message
            ψh = @. -(k2D + Ld^-2).^-1 * qh

            ψ = real.(ifft(ψh))

            if isnan(ψ[2,2])
                error("Psi is NaN")
            else
                u, v = u_from_psi_fft(ψh, kx, ky)

                cfl = dt * maximum([maximum(u) / dx, maximum(v) / dy])

                elapsed_time = time() - start_time

                # modified from GophysicalFlows.jl ex
                log = @sprintf("step: %04d, t: %.1f, cfl: %.2f, KE1 avg.: %.4e, ens1: %.4e, walltime: %.2f min",
                n, (t0+n*dt)/3600/24, cfl, mean(u.^2 .+ v.^2), sum(real.(ifft(-k2D .* ψh)).^2), elapsed_time/60)

                println(log)

            end
        end

        # if mod(n, save_every) == 0          # save streamfunction fields
        #     if save_bool==true # && n >= start_saving*nt
        #         ψ1, ψ2 = invert_qg_pv(q1, q2, ψ1_bg, ψ2_bg, inversion_ops, dx, dy) # (q1, q2)

        #         save_streamfunction(save_path, ψ1[:,save_ind_start:save_ind_end], ψ2[:,save_ind_start:save_ind_end], t0+n*dt, params)
        #         cnt+=1

        #     end
        # end

        if mod(n, plot_every) == 0          # plot whatever is in save_basic_anim_panel() function
            ψh = @. -(k2D + Ld^-2).^-1 * qh

            # ψ_f = wavemaker(ψ0, x, y, x0, y0, δx, δy, t, τ)
            # q_f = real.(ifft((k2D .+ Ld^-2) .* fft(ψ_f)))        

            save_BT_anim_panel(fig_path, ell, real.(ifft(qh)), real.(ifft(ψh)))

            ell+=1
        end

    end

    if save_last==true
        ψh = @. -(k2D + Ld^-2).^-1 * qh

        save_streamfunction(save_path, real.(ifft(ψh)), t0+nt*dt, params)
    end

    # To turn the PyPlot GUI back on
    PyPlot.pygui(true)
end



function run_BT_model_FD(q, t, params; timestepper="RK4", output_every=500)
    start_time = time()
    # To turn off the PyPlot GUI
    PyPlot.pygui(false)

    # calculating y indices for prescribed meridional width of save domain
    save_ind_start = floor(Int, Ny * y_width / 2)
    save_ind_end   = floor(Int, Ny * (1 - y_width / 2))

    # defining a couple of counters
    cnt=1
    ell=1

    for n = 1:nt

        t+= dt

        q = rk4_BT_FD(q, t, dt)

        if mod(n, output_every) == 0      # output a message
            Q = q .- f0
            Q .-= mean(Q)           # ensure zero mean
            Q_vec = reshape(Q, Nx*Ny)
            ψ = reshape(L_re_FD_inv \ Q_vec, Nx, Ny)

            if isnan(ψ[2,2])
                error("Psi is NaN")
            else
                u, v = u_from_psi_fft(fft(ψ), kx, ky)

                cfl = dt * maximum([maximum(u) / dx, maximum(v) / dy])

                elapsed_time = time() - start_time

                # modified from GophysicalFlows.jl ex
                log = @sprintf("step: %04d, t: %.1f, cfl: %.2f, KE1 avg.: %.4e, ens1: %.4e, walltime: %.2f min",
                n, (t0+n*dt)/3600/24, cfl, mean(u.^2 .+ v.^2), sum(reshape(L_re_FD * reshape(ψ, Nx*Ny), Nx, Ny).^2), elapsed_time/60)

                println(log)

            end
        end

        # if mod(n, save_every) == 0          # save streamfunction fields
        #     if save_bool==true # && n >= start_saving*nt
        #         ψ1, ψ2 = invert_qg_pv(q1, q2, ψ1_bg, ψ2_bg, inversion_ops, dx, dy) # (q1, q2)

        #         save_streamfunction(save_path, ψ1[:,save_ind_start:save_ind_end], ψ2[:,save_ind_start:save_ind_end], t0+n*dt, params)
        #         cnt+=1

        #     end
        # end

        if mod(n, plot_every) == 0          # plot whatever is in save_basic_anim_panel() function
            Q = q .- f0
            Q .-= mean(Q)           # ensure zero mean
            Q_vec = reshape(Q, Nx*Ny)
            ψ = reshape(L_re_FD_inv \ Q_vec, Nx, Ny)

            # ψ_f = wavemaker(ψ0, x, y, x0, y0, δx, δy, t, τ)
            # q_f = real.(ifft((k2D .+ Ld^-2) .* fft(ψ_f)))        

            save_BT_anim_panel(fig_path, ell, q, ψ)

            ell+=1
        end

    end

    if save_last==true
        Q = q .- f0
        Q .-= mean(Q)           # ensure zero mean
        Q_vec = reshape(Q, Nx*Ny)
        ψ = reshape(L_re_FD_inv \ Q_vec, Nx, Ny)

        save_streamfunction(save_path, ψ, t0+nt*dt, params)
    end

    # To turn the PyPlot GUI back on
    PyPlot.pygui(true)
end



function init_params_BT_FD(Nx, Ny, dx, dy, r, f0, L_re_FD_inv, L_re_FD, Lap_op; backend=CPU())
    ip = [i == Nx ? 1 : i + 1 for i in 1:Nx]
    im = [i == 1  ? Nx : i - 1 for i in 1:Nx]
    jp = [j == Ny ? 1 : j + 1 for j in 1:Ny]
    jm = [j == 1  ? Ny : j - 1 for j in 1:Ny]

    L_fac = lu(L_re_FD_inv)    # huge speedup

    return BTParamsOpt(
        Nx, Ny, dx, dy, r, f0,
        L_fac, L_re_FD, Lap_op,
        ip, im, jp, jm,
        backend
    )
end

function run_BT_model_FD_KA(q, params; nt=10_000, dt, output_every=500)

    Nx, Ny = params.Nx, params.Ny

    # Preallocate all arrays once
    Q      = zeros(Nx, Ny)
    Q_vec  = zeros(Nx*Ny)
    ψ_v    = zeros(Nx, Ny)
    ψ_v_vec= zeros(Nx*Ny)
    ψ_f    = zeros(Nx, Ny)
    ψ_f_vec= zeros(Nx*Ny)
    q_f    = zeros(Nx, Ny)
    dqdt   = zeros(Nx, Ny)
    ψ_sum = similar(ψ_v)
    q_sum = similar(q)

    η = stationary_topo(η0, x, y, -x0, y0, δx, δy)

    lap_q_vec = zeros(Nx*Ny)
    lap2_q_vec = zeros(Nx*Ny)

    workspace = (ψ_v, ψ_v_vec, Q, Q_vec, ψ_f, ψ_f_vec, q_f, ψ_sum, q_sum, lap_q_vec, lap2_q_vec, η)

    # RK4 storage
    k1 = similar(q)
    k2 = similar(q)
    k3 = similar(q)
    k4 = similar(q)
    qnew = similar(q)

    t = 0.0

    start_time = time()
    # To turn off the PyPlot GUI
    PyPlot.pygui(false)

    # calculating y indices for prescribed meridional width of save domain
    save_ind_start = floor(Int, Ny * y_width / 2)
    save_ind_end   = floor(Int, Ny * (1 - y_width / 2))

    # defining a couple of counters
    cnt=1
    ell=1

    for n = 1:nt

        rk4_BT_FD!(q, qnew, k1, k2, k3, k4, workspace, params, dt, t)
        t += dt

        if mod(n, output_every) == 0      # output a message
            Q = q .- f0
            Q .-= mean(Q)           # ensure zero mean
            Q_vec = reshape(Q, Nx*Ny)
            ψ = reshape(L_re_FD_inv \ Q_vec, Nx, Ny)

            if isnan(ψ[2,2])
                error("Psi is NaN")
            else
                u, v = u_from_psi_fft(fft(ψ), kx, ky)

                cfl = dt * maximum([maximum(u) / dx, maximum(v) / dy])

                elapsed_time = time() - start_time

                # modified from GophysicalFlows.jl ex
                log = @sprintf("step: %04d, t: %.1f, cfl: %.2f, KE1 avg.: %.4e, ens1: %.4e, walltime: %.2f min",
                n, (t0+n*dt)/3600/24, cfl, mean(u.^2 .+ v.^2), sum(reshape(L_re_FD * reshape(ψ, Nx*Ny), Nx, Ny).^2), elapsed_time/60)

                println(log)

            end
        end

        # if mod(n, save_every) == 0          # save streamfunction fields
        #     if save_bool==true # && n >= start_saving*nt
        #         ψ1, ψ2 = invert_qg_pv(q1, q2, ψ1_bg, ψ2_bg, inversion_ops, dx, dy) # (q1, q2)

        #         save_streamfunction(save_path, ψ1[:,save_ind_start:save_ind_end], ψ2[:,save_ind_start:save_ind_end], t0+n*dt, params)
        #         cnt+=1

        #     end
        # end

        if mod(n, plot_every) == 0          # plot whatever is in save_basic_anim_panel() function
            Q = q .- f0
            Q .-= mean(Q)           # ensure zero mean
            Q_vec = reshape(Q, Nx*Ny)
            ψ = reshape(L_re_FD_inv \ Q_vec, Nx, Ny)

            # ψ_f = wavemaker(ψ0, x, y, x0, y0, δx, δy, t, τ)
            # q_f = real.(ifft((k2D .+ Ld^-2) .* fft(ψ_f)))        

            save_BT_anim_panel(fig_path, ell, q, ψ)

            ell+=1
        end

    end

    if save_last==true
        Q = q .- f0
        Q .-= mean(Q)           # ensure zero mean
        Q_vec = reshape(Q, Nx*Ny)
        ψ = reshape(L_re_FD_inv \ Q_vec, Nx, Ny)

        save_streamfunction(save_path, ψ, t0+nt*dt, params)
    end

    # To turn the PyPlot GUI back on
    PyPlot.pygui(true)
end
