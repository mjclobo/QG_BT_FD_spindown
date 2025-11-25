
################################################################################
# RK4 time steppers
################################################################################

# function rk4_step!(qh, t, dt,
#     rhs!, ν, r, τ, f0, kx, ky, params)

#     K1 = similar(qh)
#     K2 = similar(qh)
#     K3 = similar(qh)
#     K4 = similar(qh)
#     qtmp = similar(qh)

#     rhs!(K1, qh, t, ν, r, τ, f0, kx, ky, params)

#     @. qtmp = qh + 0.5*dt*K1
#     rhs!(K2, qtmp, t + 0.5*dt, ν, r, τ, f0, kx, ky, params)

#     @. qtmp = qh + 0.5*dt*K2
#     rhs!(K3, qtmp, t + 0.5*dt, ν, r, τ, f0, kx, ky, params)

#     @. qtmp = qh + dt*K3
#     rhs!(K4, qtmp, t + dt, ν, r, τ, f0, kx, ky, params)

#     @. qh += dt * (K1 + 2K2 + 2K3 + K4)/6

#     return nothing
# end


function rk4_BT(qh, t, dt)
    k1qh = rhs(qh, t)
    k2qh = rhs(qh .+ 0.5dt .* k1qh, t+0.5*dt)
    k3qh = rhs(qh .+ 0.5dt .* k2qh, t+0.5*dt)
    k4qh = rhs(qh .+ dt .* k3qh, t+dt)

    qh_new = qh .+ dt/6 .* (k1qh .+ 2k2qh .+ 2k3qh .+ k4qh)
    return qh_new
end

function rk4_BT_FD(q, t, dt)
    k1qh = rhs_BT_FD(q, t)
    k2qh = rhs_BT_FD(q .+ 0.5dt .* k1qh, t+0.5*dt)
    k3qh = rhs_BT_FD(q .+ 0.5dt .* k2qh, t+0.5*dt)
    k4qh = rhs_BT_FD(q .+ dt .* k3qh, t+dt)

    q_new = q .+ dt/6 .* (k1qh .+ 2k2qh .+ 2k3qh .+ k4qh)
    return q_new
end

function rk4_BT_FD!(q, qnew, k1, k2, k3, k4, workspace, params, dt, t)
    # workspace is a NamedTuple of (ψ_v, ψ_v_vec, Q, …)
    k1 .= rhs_BT_FD!(k1, q, workspace..., params, t)

    @. qnew = q + 0.5*dt*k1
    k2 .= rhs_BT_FD!(k2, qnew, workspace..., params, t + 0.5dt)

    @. qnew = q + 0.5*dt*k2
    k3 .= rhs_BT_FD!(k3, qnew, workspace..., params, t + 0.5dt)

    @. qnew = q + dt*k3
    k4 .= rhs_BT_FD!(k4, qnew, workspace..., params, t + dt)

    @. q = q + dt/6*(k1 + 2k2 + 2k3 + k4)
    return q
end


function rk4(q1, q2, dt)
    k1q1, k1q2 = rhs(q1, q2)
    k2q1, k2q2 = rhs(q1 .+ 0.5dt .* k1q1, q2 .+ 0.5dt .* k1q2)
    k3q1, k3q2 = rhs(q1 .+ 0.5dt .* k2q1, q2 .+ 0.5dt .* k2q2)
    k4q1, k4q2 = rhs(q1 .+ dt .* k3q1, q2 .+ dt .* k3q2)

    q1_new = q1 .+ dt/6 .* (k1q1 .+ 2k2q1 .+ 2k3q1 .+ k4q1)
    q2_new = q2 .+ dt/6 .* (k1q2 .+ 2k2q2 .+ 2k3q2 .+ k4q2)
    return q1_new, q2_new
end


