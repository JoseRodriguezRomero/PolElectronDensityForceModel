using Plots, SpecialFunctions, LaTeXStrings, PGFPlotsX

include("GetPotentialFromAngles.jl")

λ1 = 2.45;
λ2 = 1.65;

λ = (λ1*λ2)/(λ1+λ2);
l = collect((10.0/300.0).*(0:300));
aux_exp = exp.(-λ.*(l.^2.0));

num_pts = length(l);

l = zeros(Float64,num_pts,22) .+ l;
l = l .^ Float64.(collect(1:22)');

λ = zeros(Float64,22) .+ λ;
λ = λ .^ Float64.(collect(1:22));

l0 = zeros(Float64,20);

x0 = zeros(Float64,num_pts);
x1 = zeros(Float64,num_pts);
x2 = zeros(Float64,num_pts);
x3 = zeros(Float64,num_pts);
x4 = zeros(Float64,num_pts);
x5 = zeros(Float64,num_pts);
x5 = zeros(Float64,num_pts);
x6 = zeros(Float64,num_pts);
x7 = zeros(Float64,num_pts);
x8 = zeros(Float64,num_pts);
x9 = zeros(Float64,num_pts);
x10 = zeros(Float64,num_pts);
x11 = zeros(Float64,num_pts);

for i in 1:num_pts
    x0[i] = XCOrder0(λ,l[i,:])./XCOrder0(λ,l0);
    x1[i] = (aux_exp[i].*XCOrder1(λ))./XCOrder1(λ);
    x2[i] = (aux_exp[i].*XCOrder2(λ,l[i,:]))./XCOrder2(λ,l0);
    x3[i] = (aux_exp[i].*XCOrder3(λ,l[i,:]))./XCOrder3(λ,l0);
    x4[i] = (aux_exp[i].*XCOrder4(λ,l[i,:]))./XCOrder4(λ,l0);
    x5[i] = (aux_exp[i].*XCOrder5(λ,l[i,:]))./XCOrder5(λ,l0);
    x6[i] = (aux_exp[i].*XCOrder6(λ,l[i,:]))./XCOrder6(λ,l0);
    x7[i] = (aux_exp[i].*XCOrder7(λ,l[i,:]))./XCOrder7(λ,l0);
    x8[i] = (aux_exp[i].*XCOrder8(λ,l[i,:]))./XCOrder8(λ,l0);
    x9[i] = (aux_exp[i].*XCOrder9(λ,l[i,:]))./XCOrder9(λ,l0);
    x10[i] = (aux_exp[i].*XCOrder10(λ,l[i,:]))./XCOrder10(λ,l0);
    x11[i] = (aux_exp[i].*XCOrder11(λ,l[i,:]))./XCOrder11(λ,l0);
end

p1 = plot(l[:,1],x0,label=L"n = 0");
plot!(l[:,1],x1,label=L"n = 1");
plot!(l[:,1],x2,label=L"n = 2");
plot!(l[:,1],x3,label=L"n = 3");
plot!(l[:,1],x4,label=L"n = 4");
plot!(l[:,1],x5,label=L"n = 5");
plot!(l[:,1],x6,label=L"n = 6");
plot!(l[:,1],x7,label=L"n = 7");
plot!(l[:,1],x8,label=L"n = 8");
plot!(l[:,1],x9,label=L"n = 9");
plot!(l[:,1],x10,label=L"n = 10");
plot!(l[:,1],x11,label=L"n = 11");
plot!(xlims=(0,5));
plot!(ylims=(-0.25,1));

plot!(xticks=(0:5,[]));
plot!(ylabel=L"\overline{\mathrm{XC}^{\mathrm{EE} \left( n \right)}_\mathrm{AB}}");
plot!(legend=:outertop,legendcolumns=6);
plot!(size=(525,425));

PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.DEFAULT_PREAMBLE; "\\usepackage[bitstream-charter]{mathdesign}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{amsmath,accents}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{accents}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[utf8]{inputenc}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[T1]{fontenc}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{fontspec}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[fontsize=10pt]{fontsize}"];

fp1 = @pgf Axis(
    {
        legend_style =
        {
            at = Coordinate(0.5,1.3),
            anchor = "center",
            legend_columns = 6
        },
        xmajorgrids,
        ymajorgrids,
        xmin = 0,
        xmax = 5,
        ymin = -0.25,
        ymax = 1.0,
        ytick = "-0.25,0,0.25,0.5,0.75,1",
        yticklabels = L"$-0.25$,$0.00$,$0.25$,$0.50$,$0.75$,$1.00$",
        xtick = "-0.25,0,1,2,3,4,5",
        xticklabels = "\\empty",
        # xlabel = L"$d$",
        ylabel = L"\overline{\mathrm{XC}^{\mathrm{EE} \left( k \right)}_\mathrm{AB}}",
        width = 300, 
        height = 135,
        "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
    },
    Plot(
        {
            color = theme_palette(:auto)[1],
            style = {"thick"},
            no_marks,
        },
        Coordinates(l[:,1],x0)
    ),
    Plot(
        {
            color = theme_palette(:auto)[2],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x1),
    ),
    Plot(
        {
            color = theme_palette(:auto)[3],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x2),
    ),
    Plot(
        {
            color = theme_palette(:auto)[4],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x3),
    ),
    Plot(
        {
            color = theme_palette(:auto)[5],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x4),
    ),
    Plot(
        {
            color = theme_palette(:auto)[6],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x5),
    ),
    Plot(
        {
            color = theme_palette(:auto)[7],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x6),
    ),
    Plot(
        {
            color = theme_palette(:auto)[8],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x7),
    ),
    Plot(
        {
            color = theme_palette(:auto)[9],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x8),
    ),
    Plot(
        {
            color = theme_palette(:auto)[10],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x9),
    ),
    Plot(
        {
            color = theme_palette(:auto)[11],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x10),
    ),
    Plot(
        {
            color = theme_palette(:auto)[12],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x11),
    ),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $k = 0$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $k = 1$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $k = 2$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $k = 3$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $k = 4$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $k = 5$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $k = 6$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $k = 7$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $k = 8$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $k = 9$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $k = 10$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $k = 11$"),
);

function aux_mult(λ::Real,n::Int)
    aux_ret = gamma(5/2)/gamma(5/2+n-1);
    aux_ret *= (λ^(2*n+1))/(2^(2*n+1));
    return ((-1)^(n+1))*aux_ret;
end

for i in 1:num_pts
    x0[i] = aux_mult(λ[1],0).*(XCOrderD00(λ,l[i,:]).*aux_exp[i]);
    x0[i] += aux_mult(λ[1],0).*XCOrderD01(λ,l[i,:]);

    x1[i] = aux_mult(λ[1],1).*(XCOrderD1(λ,l[i,:]).*aux_exp[i]);
    x2[i] = aux_mult(λ[1],2).*(XCOrderD2(λ,l[i,:]).*aux_exp[i]);
    x3[i] = aux_mult(λ[1],3).*(XCOrderD3(λ,l[i,:]).*aux_exp[i]);
    x4[i] = aux_mult(λ[1],4).*(XCOrderD4(λ,l[i,:]).*aux_exp[i]);
    x5[i] = aux_mult(λ[1],5).*(XCOrderD5(λ,l[i,:]).*aux_exp[i]);
    x6[i] = aux_mult(λ[1],6).*(XCOrderD6(λ,l[i,:]).*aux_exp[i]);
    x7[i] = aux_mult(λ[1],7).*(XCOrderD7(λ,l[i,:]).*aux_exp[i]);
    x8[i] = aux_mult(λ[1],8).*(XCOrderD8(λ,l[i,:]).*aux_exp[i]);
    x9[i] = aux_mult(λ[1],9).*(XCOrderD9(λ,l[i,:]).*aux_exp[i]);
    x10[i] = aux_mult(λ[1],10).*(XCOrderD10(λ,l[i,:]).*aux_exp[i]);
    x11[i] = aux_mult(λ[1],11).*(XCOrderD11(λ,l[i,:]).*aux_exp[i]);
end

p2 = plot(l[:,1],x0,label=L"n = 0");
plot!(l[:,1],x1,label=L"n = 1");
plot!(l[:,1],x2,label=L"n = 2");
plot!(l[:,1],x3,label=L"n = 3");
plot!(l[:,1],x4,label=L"n = 4");
plot!(l[:,1],x5,label=L"n = 5");
plot!(l[:,1],x6,label=L"n = 6");
plot!(l[:,1],x7,label=L"n = 7");
plot!(l[:,1],x8,label=L"n = 8");
plot!(l[:,1],x9,label=L"n = 9");
plot!(l[:,1],x10,label=L"n = 10");
plot!(l[:,1],x11,label=L"n = 11");
plot!(xlims=(0,5));
plot!(ylims=(-0.05,0.4));

plot!(xlabel=L"d")
plot!(legend=false);
plot!(size=(525,50));
plot!(ylabel=L"\overline{\mathrm{XC}^{\mathrm{ED} \left( n \right)}_\mathrm{AB}}");

plot(p1,p2,layout=grid(2,1,heights=(4.5/8,3.5/8)));
# plot!(left_margin=3Plots.mm)

fp2 = @pgf Axis(
    {
        "legend pos = north east",
        xmajorgrids,
        ymajorgrids,
        xmin = 0,
        xmax = 5,
        ymin = -0.05,
        ymax = 0.4,
        ytick = "0,0.1,0.2,0.3,0.4",
        yticklabels = L"$0.0$,$0.1$,$0.2$,$0.3$,$0.4$",
        xtick = "0,1,2,3,4,5",
        xlabel = L"$d$",
        ylabel = L"\overline{\mathrm{XC}^{\mathrm{ED} \left( k \right)}_\mathrm{AB}}",
        width = 300, 
        height = 135,
        "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
    },
    Plot(
        {
            color = theme_palette(:auto)[1],
            style = {"thick"},
            no_marks,
        },
        Coordinates(l[:,1],x0)
    ),
    Plot(
        {
            color = theme_palette(:auto)[2],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x1),
    ),
    Plot(
        {
            color = theme_palette(:auto)[3],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x2),
    ),
    Plot(
        {
            color = theme_palette(:auto)[4],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x3),
    ),
    Plot(
        {
            color = theme_palette(:auto)[5],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x4),
    ),
    Plot(
        {
            color = theme_palette(:auto)[6],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x5),
    ),
    Plot(
        {
            color = theme_palette(:auto)[7],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x6),
    ),
    Plot(
        {
            color = theme_palette(:auto)[8],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x7),
    ),
    Plot(
        {
            color = theme_palette(:auto)[9],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x8),
    ),
    Plot(
        {
            color = theme_palette(:auto)[10],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x9),
    ),
    Plot(
        {
            color = theme_palette(:auto)[11],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x10),
    ),
    Plot(
        {
            color = theme_palette(:auto)[12],
            no_marks,
            style = {"thick"},
        },
        Coordinates(l[:,1],x11),
    )
);

@pgf GroupPlot({group_style = { group_size = "1 by 2"}},fp1,fp2)