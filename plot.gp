# Set the terminal to output a PNG image (optional, remove to view interactively)
# set terminal pngcairo size 800,600
# set output 'plot.png'

# Add titles and labels
set title "Plot of tbtest.dat"
set xlabel "X-axis (Column 1)"
set ylabel "Y-axis (Columns 2 and 3)"

# Enable a grid for easier reading
set grid

# Adjust the axis ranges slightly so the points don't touch the very edges
set xrange [-2:2]
set yrange [-1.5:1.5]

# Plot the data
# using 1:2 plots column 2 against column 1
# using 1:3 plots column 3 against column 1
plot 'tbtest.dat' using 1:2 with linespoints linewidth 2 pointtype 7 title 'Column 2', \
     'tbtest.dat' using 1:3 with linespoints linewidth 2 pointtype 5 title 'Column 3'
