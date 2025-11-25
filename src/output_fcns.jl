
################################################################################
# Data write functions
################################################################################

function struct_to_string(s::T) where T
    field_strings = String[]
    for field_name in fieldnames(T)
        field_value = getfield(s, field_name)
        if startswith(string(field_name), "Δρ")
            push!(field_strings, "drho$field_value")
        elseif startswith(string(field_name), "H") && length(field_value)>1
            push!(field_strings, "H1"*string(field_value[1]))
            push!(field_strings, "H2"*string(field_value[2]))
        else
            push!(field_strings, "$field_name"*string(round(field_value, digits=3)))
        end
    end
    return join(field_strings, "_")
end

function save_streamfunction(dir, ψ, t, params)

    file_name =  "BT_test_file_t$t.jld"
    # Save variables to JLD file
    jld_data = Dict("ψ" => Array(ψ), "t" => t)
    jldsave(dir * file_name; jld_data)

    println("Saved streamfunction to $file_name")

end


function save_BT_anim_panel(fig_path, ell, q, ψ)
    plotname = "snapshots"
    # ψ1, ψ2 = invert_qg_pv(q1, q2, ψ1_bg, ψ2_bg, inversion_ops, dx, dy) # (q1, q2)

    fig, ax = plt.subplots(1,1, figsize=(8, 6.5))

    # lim = maximum(abs.(ψ .- mean(ψ, dims=1)))
    # pc = ax.pcolormesh(x, y, (ψ .- mean(ψ, dims=1))', cmap=PyPlot.cm.bwr, vmin=-lim, vmax=lim)
    # ax.set_title("ψ")

    # lim = maximum(abs.(q .- mean(q, dims=1)))
    pc = ax.pcolormesh(x ./ Ld, y ./ Ld, (q./f0)', cmap=PyPlot.cm.RdYlBu , vmin=0.975, vmax=1.45)

    ax.set_title(L"q / f_0")

    plt.colorbar(pc)

    ax.set_xlabel("x / Ld")
    ax.set_ylabel("y / Ld")

    ###
    local savename = @sprintf("%s_%04d.png", joinpath(fig_path, plotname), ell)
    PyPlot.savefig(savename)

    PyPlot.close()
end


