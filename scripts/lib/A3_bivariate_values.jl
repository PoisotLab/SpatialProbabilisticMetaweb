function bivariatelayer(l1, l2; n_stops=3, rscale=0.0:0.01:1.0)
    # Rescale layers
    r1 = convert(SimpleSDMResponse, rescale(l1, rscale));
    r2 = convert(SimpleSDMResponse, rescale(l2, rscale));

    # Get classes
    r1[keys(l1)] = round.((n_stops - 1) .* values(r1); digits=0) .+ 1;
    r2[keys(l2)] = round.((n_stops - 1) .* values(r2); digits=0) .+ 1;

    # Convert as integers
    r1 = convert(Int64, r1)
    r2 = convert(Int64, r2)

    # Get bivariate combinations
    function bivariator(n1, n2)
        function bv(v1, v2)
            return n2 * (v2 - 1) + v1
        end
        return bv
    end
    b = bivariator(n_stops, n_stops).(r1, r2)

    return b
end

function _get_bivariate_colormap(;
    n_stops=3,
    p0 = colorant"#e8e8e8ff",
    p1 = colorant"#5ac8c8ff",
    p2 = colorant"#be64acff",
    rev = false
)
    # Generate colormap
    cm1 = LinRange(p0, p1, n_stops)
    cm2 = LinRange(p0, p2, n_stops)
    if rev
        cmat = ColorBlendModes.BlendMultiply.(cm2, cm1')
    else
        cmat = ColorBlendModes.BlendMultiply.(cm1, cm2')
    end
    cmap = vec(cmat)
end

function bivariateheatmap!(
    f,
    l1,
    l2;
    n_stops=3,
    rscale=0.0:0.01:1.0,
    p0 = colorant"#e8e8e8ff",
    p1 = colorant"#5ac8c8ff",
    p2 = colorant"#be64acff",
    shading=false,
    rev=false,
    cmap=nothing
)
    # Generate colormap
    if isnothing(cmap)
        cmap = _get_bivariate_colormap(n_stops=3, p0=p0, p1=p1, p2=p2, rev=rev)
    end

    # Get bivariate layer
    b = bivariatelayer(l1, l2; n_stops=n_stops, rscale=rscale)

    # Create bivariate figure
    # m_biv = Axis(f[1, 1:3]; aspect = DataAspect())
    heatmap!(f, b; colormap = cmap, colorrange = (1, n_stops * n_stops))

    return f
end

function bivariatesurface!(
    f,
    l1,
    l2;
    n_stops=3,
    rscale=0.0:0.01:1.0,
    p0 = colorant"#e8e8e8ff",
    p1 = colorant"#5ac8c8ff",
    p2 = colorant"#be64acff",
    shading=false,
    rev=false,
    cmap=nothing
)
    # Generate colormap
    if isnothing(cmap)
        cmap = _get_bivariate_colormap(n_stops=n_stops, p0=p0, p1=p1, p2=p2, rev=rev)
    end

    # Get bivariate layer
    b = bivariatelayer(l1, l2; n_stops=n_stops, rscale=rscale)

    # Create bivariate figure
    # m_biv = Axis(f[1, 1:3]; aspect = DataAspect())
    surface!(f, b; colormap = cmap, colorrange = (1, n_stops * n_stops), shading=shading)

    return f
end

function bivariatelegend!(
    ax, l1, l2;
    n_stops=3,
    p0 = colorant"#e8e8e8ff",
    p1 = colorant"#5ac8c8ff",
    p2 = colorant"#be64acff",
    rev = false,
    cmap = nothing
)
    # Generate colormap
    if isnothing(cmap)
        cmap = _get_bivariate_colormap(n_stops=n_stops, p0=p0, p1=p1, p2=p2, rev=rev)
    end

    # Add bivariate legend
    # m_leg = Axis(f[1, 4]; aspect = 1, xlabel = "Temperature", ylabel = "Precipitation")
    x = LinRange(minimum(l1), maximum(l1), n_stops)
    y = LinRange(minimum(l2), maximum(l2), n_stops)
    heatmap!(ax, x, y, reshape(1:(n_stops * n_stops), (n_stops, n_stops)); colormap = cmap)
end
