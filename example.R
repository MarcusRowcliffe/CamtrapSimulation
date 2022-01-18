devtools::source_url("https://raw.githubusercontent.com/MarcusRowcliffe/CamtrapSimulation/master/CamtrapSimulation.R")

#Create a correlated random walk movement path
path <- pathgen(5e3, kTurn=2, kCor=TRUE, pTurn=1, 
                logspeed=-2, speedSD=1, speedCor=0, 
                xlim=c(0,10), wrap=TRUE)

#Create a camera detection zone
dz <- data.frame(x=5, y=2, r=6, th=1, dir=0)

#Visualise
plot_wrap(path, lineargs = list(col="grey"))
plot_dzone(dz, border=2)

#Create position data for sequences falling within the detection zone
posdat <- sequence_data(path, dz)
points(posdat$x, posdat$y, col=2, pch=16, cex=0.5)

#Create speed data summarised for each sequence
seqdat <- calc_speed(posdat)
