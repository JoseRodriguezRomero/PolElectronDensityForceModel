using LaTeXStrings, Plots, PGFPlotsX

PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.DEFAULT_PREAMBLE; "\\usepackage[bitstream-charter]{mathdesign}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[T1]{fontenc}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{fontspec}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[fontsize=10pt]{fontsize}"];

function ReadFile(file_name::String)
    data = zeros(Float64,0);
    raw_data = readlines(file_name);

    for value in split(raw_data[1])
        push!(data,parse(Float64,value));
    end

    return data;
end

f1_data = ReadFile("f1_data.txt");
f2_data = ReadFile("f2_data.txt");

plot((0:(length(f1_data)-1)).*(3.0/(length(f1_data)-1)),f1_data,
    label=L"$K = 100$",linestyle=:solid,linewidth=2);
plot!((0:(length(f2_data)-1)).*(3.0/(length(f2_data)-1)),f2_data,
    label=L"$K = 10$",linestyle=:dash,linewidth=2)

plot!(xlims=(0,3),ylims=(-10,40));
plot!(xlabel=L"$d$",ylabel=L"$\sum^K_{k = 1} \frac{"*
    L"\left( -1 \right)^k \mathrm{XC}_\mathrm{AB}^{\mathrm{EN} \left( k \right)}}{\left(2 k \right)! \mathrm{exp} \left( - \lambda d^2 \right)}$");

plot!(left_margin=8Plots.mm);
plot!(right_margin=2Plots.mm);
plot!(bottom_margin=5Plots.mm);
plot!(top_margin=2Plots.mm);
plot!(size=(600,240))

# savefig("inf_lim.pdf")

PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.DEFAULT_PREAMBLE; "\\usepackage{amsmath,accents}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{accents}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[utf8]{inputenc}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[T1]{fontenc}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{fontspec}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{mathptmx}"];
PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[fontsize=10pt]{fontsize}"];

# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.DEFAULT_PREAMBLE; "\\usepackage[bitstream-charter]{mathdesign}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{amsmath,accents}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{accents}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[utf8]{inputenc}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[T1]{fontenc}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage{fontspec}"];
# PGFPlotsX.CUSTOM_PREAMBLE = [PGFPlotsX.CUSTOM_PREAMBLE; "\\usepackage[fontsize=10pt]{fontsize}"];

@pgf Axis(
    {
        "legend pos = north east",
        xmajorgrids,
        ymajorgrids,
        xmin = 0,
        xmax = 3,
        ymin = -10,
        ymax = 40,
        xtick = "{0,1,2,3}",
        ytick = "{-10,0,10,20,30,40}",
        xlabel = L"$d$",
        ylabel = L"$\displaystyle \sum^K_{k = 1} \frac{\left( - 1 \right)^k \mathrm{XC}_\mathrm{AB}^{\mathrm{EN} \left( k \right)}}{\left( "*
            L"2 k \right)! \exp \left( - \lambda d^2 \right)}$",
        width = 300, 
        height = 150,
        "grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!25}",
    },
    Plot(
        {
            color = theme_palette(:auto)[1],
            style = {"thick"},
            no_marks,
        },
        Coordinates((0:(length(f1_data)-1)).*(3.0/(length(f1_data)-1)),f1_data)
    ),
    Plot(
        {
            color = theme_palette(:auto)[2],
            no_marks,
            style = {"thick","dashed"},
        },
        Coordinates((0:(length(f1_data)-1)).*(3.0/(length(f1_data)-1)),f2_data),
    ),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $K = 100$"),
    LegendEntry({anchor = "west"},L"\fontsize{8pt}{8pt}\selectfont $K = 10$")
)
