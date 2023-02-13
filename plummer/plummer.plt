reset
set term gif animate
set output "plummer.2.gif"

set datafile separator ","

set xrange [-3:3]
set yrange [-3:3]
set zrange [-3:3]
set view equal xyz
set ticslevel 0

num = 1000
last_chunk_index =101
skip = 1

do for [i=0:last_chunk_index:skip]{
    start = i*(num+3)+1
    stop  = (num*3)*(i+1) - 2

    time = system(sprintf('awk "NR == %d" %s',start,"plummer.csv"))

    set title sprintf('t = %s',time)

    splot "plummer.csv" every ::1 index i using 1:2:3 title "point" with points pt 7 lc 6 pointsize 0.3
}

if(last_chunk_index%skip != 0){
    i     = last_chunk_index
    start = i*(num+3)+1
    stop  = (num*3)*(i+1) - 2

    time = system(sprintf('awk "NR == %d" %s',start,"plummer.csv"))

    set title sprintf('t = %s',time)

    splot "plummer.csv" every ::1 index i using 1:2:3 title "point" with points pt 7 lc 6 pointsize 0.1
}

print("")
