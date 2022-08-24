# Gnuplot script for plotting the band structure as bitmap
# You need to generate *.bin file with ubs_bmp.m prior to plotting

#set terminal postscript eps color enhanced font "Whitney-Book,16" size 20cm,12.0cm
set terminal postscript eps color enhanced font "Helvetica,16" size 20cm,12cm


set output "f2b-band-structure.eps"
set encoding utf8 # nice minus
set minussign # nice minus
set multiplot layout 2,3

# PLOT 1
#set origin 0.0,0.0
#set size 1.0,0.33
set xlabel "Wave vector"
set ylabel "Energy (eV)"
set format y "%.1f" # 1 digit after period
set border linewidth 1.0
#set datafile separator "," # for *.csv file
set yrange [-2.0:0.5]
set xrange [0:1.43966]
#set cbrange [0:4]
set xtics ("G" 0, "X" 0.42118, "S" 0.84236, "G" 1.43966)
set border lw 2 lc rgb "white"
set xtics tc rgb "black"
set ytics tc rgb "black"
set cbtics tc rgb "black"
set tics front
set label 1 "SrIrO3 (2x2x2 -> 1x1x1)" at graph 0.05,0.95 front textcolor "white"
plot "case.f2b.bin" binary matrix with image title "",\
     0 title "" lc "white" lw 2 lt 2 dt 2

# PLOT 2
#set label 1 "GaP-wz (not opt.)" at graph 0.05,0.95 front textcolor "white"
#set xrange [0:0.34095]
#set xtics ("L" 0, "G" 0.15824, "X" 0.34095)
#plot "GaP-wz-supecell-notOptimized.bin" binary matrix with image title "",\
#     0 title "" lc "white" lw 2 lt 2 dt 2

# PLOT 3
#set label 1 "GaP-wz" at graph 0.05,0.95 front textcolor "white"
#set xrange [0:0.33995]
#set xtics ("L" 0, "G" 0.15725, "X" 0.33995)
#plot "GaP-wz-supecell-optimized.bin" binary matrix with image title "",\
#     0 title "" lc "white" lw 2 lt 2 dt 2

# PLOT 3
#set label 1 "GaN-wz" at graph 0.05,0.95 front textcolor "white"
#set yrange [-2.0:5.0]
#set xrange [0:0.41545]
#set xtics ("L" 0, "G" 0.19297, "X" 0.41545)
#plot "GaN-wz-supecell.bin" binary matrix with image title "",\
#     0 title "" lc "white" lw 2 lt 2 dt 2
