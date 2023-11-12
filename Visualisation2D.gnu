# ------------------------------------------------------------------------- #
# - Gnuplot two dimensional visualisation program for our system's masses - #
# ------------------------------------------------------------------------- #
#                                                                           #
#   To execute this code, write the following in the terminal:              #
#                                                                           #
#            gnuplot -persist -e "N=2" Visualisation2D.gnu                  #
#                                                                           #
#   Make sur to specify the number of masses you want to show on the plot.  #
#   To choose this number, replace "N=2" by the chosen number of masses in  #
#   the command written above. For X masses, write "N=X".                   #
#                                                                           #
# ------------------------------------------------------------------------- #

# Number of masses shown. Defaults to 2.
Masses = N

# Generic plot-tweaking
set object 1 rectangle from graph 0,0 to graph 1,1 behind fc rgb "#00008B" fs solid
set termoption font "Arial, 12"
set title sprintf("Visualisation 2D des trajectoires des %d masses", N)
set xlabel "Position en x"
set ylabel "Position en y"
# set xrange [-100:100]
# set xrange [-1000:1000]
set key box opaque


# Initialize an empty string to accumulate plot commands
plot_commands = ""

# Create a plot command in string format for each mass using their respective x and y coordinates
do for [i = 1:N] {
    
    xcoord = int(2. * i - 1)
    ycoord = int(2. * i)
    
    # Append the plot command to the string
    plot_commands = sprintf('%s, "bodies_movement.dat" using %d:%d with l t "Trajectoire masse %d"', plot_commands, xcoord, ycoord, i)
}

# Remove the leading comma and plot the accumulated commands
eval("plot".plot_commands[2:])

