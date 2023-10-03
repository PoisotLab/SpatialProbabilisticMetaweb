#### Bivariate color palettes ###

# Prepare colors
p0 = colorant"#e8e8e8"
bv_pal_1 = (p0=p0, p1=colorant"#64acbe", p2=colorant"#c85a5a")
bv_pal_2 = (p0=p0, p1=colorant"#73ae80", p2=colorant"#6c83b5")
bv_pal_3 = (p0=p0, p1=colorant"#9972af", p2=colorant"#c8b35a")
bv_pal_4 = (p0=p0, p1=colorant"#be64ac", p2=colorant"#5ac8c8")

# Colors from https://jakubnowosad.com/posts/2020-08-25-cbc-bp2/#summary
cmat2 = [
    colorant"#73ae80" colorant"#5a9178" colorant"#2a5a5b"
    colorant"#b8d6be" colorant"#90b2b3" colorant"#567994"
    colorant"#e8e8e8" colorant"#b5c0da" colorant"#6c83b5"
]
cmap2 = vec(cmat2[3:-1:1,:])
cmat5 = [
    colorant"#f3b300" colorant"#b36600" colorant"#000000"
    colorant"#f3e6b3" colorant"#b3b3b3" colorant"#376387"
    colorant"#f3f3f3" colorant"#b4d3e1" colorant"#509dc2"
]
cmap5 = vec(cmat5[3:-1:1,:])