"""
    Extract the values from a bivariate map. Mostly follows the internal codes
    from `bivariate()`. Returns a layer with the color category number for each
    pixel and a vector with the specific colors.
"""
function get_bivariate_values(
    l1,
    l2;
    classes=3,
    p0=colorant"#e8e8e8ff",
    p1=colorant"#64acbeff",
    p2=colorant"#c85a5aff",
    quantiles=true,
)
    # Verifications
    SimpleSDMLayers._layers_are_compatible(l1, l2)
    # Color ranges
    c1 = LinRange(p0, p1, classes)
    c2 = LinRange(p0, p2, classes)
    # Quantile breaks
    breakpoints = LinRange(0.0, 1.0, classes + 1)
    # Rescale values as quantiles (with 10 times the number of classes)
    if quantiles
        q1 = rescale(l1, collect(LinRange(0.0, 1.0, 10classes)))
        q2 = rescale(l2, collect(LinRange(0.0, 1.0, 10classes)))
    else
        q1 = rescale(l1, (0.0, 1.0))
        q2 = rescale(l2, (0.0, 1.0))
    end
    # Layer for the classifies values
    classified = similar(l1, Int)
    # Empty array for the colors (requires ColorBlendModes)
    cols = typeof(p0)[]
    # Get the values for all classes
    for i in 1:classes
        # Anonymous function to check in which class a value belongs (between which breakpoints)
        # We need an exception for the last class
        if isequal(classes)(i)
            fi = (v) -> breakpoints[i] < v <= breakpoints[i + 1]
        else
            fi = (v) -> breakpoints[i] <= v < breakpoints[i + 1]
        end
        # Check if the values belong is this class
        m1 = broadcast(fi, q1)
        # Now repeat for the 2nd layer
        for j in 1:classes
            if isequal(classes)(j)
                fj = (v) -> breakpoints[j] < v <= breakpoints[j + 1]
            else
                fj = (v) -> breakpoints[j] <= v < breakpoints[j + 1]
            end
            # Assign 2nd class
            m2 = broadcast(fj, q2)
            # Get the corresponding color (requires ColorBlendModes)
            push!(cols, ColorBlendModes.BlendMultiply(c1[i], c2[j]))
            # Get the locations where both classes match
            m = reduce(*, [m1, m2])
            replace!(m, false => nothing)
            # Assign value corresponding to the class-combination number
            if length(m) > 0
                classified[keys(m)] = fill(length(cols), length(m))
            end
        end
    end
    # Assign unclassified values to 1st class (if any)
    replace!(classified, 0 => 1)

    return (layer=classified, colors=cols)
end
