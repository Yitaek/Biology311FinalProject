gedata = read.csv("GEdata.csv", header=T, row.names=1)
gedata = read.csv("GEdata.csv", header=T, row.names=1)
parent = t(gedata[1:8])
RosRDel = t(gedata[9:16])
parent.t<-t(parent)
RosRDel.t<-t(RosRDel)


parent.level<-10^parent.t
RosR.level<-10^RosRDel.t
dim(parent.level)
dim(RosR.level)
diff.level<-RosR.level-parent.level
parent.level.max<-apply(parent.level.max,1,max)
parent.level.max<-apply(parent.level,1,max)
RosR.level.max<-apply(RosR.level,1,max)
diff.level.max<-RosR.level.max-parent.level.max
hist(diff.level.max)
length(which(diff.level.max>(mean(diff.level.max)+sd(diff.level.max))|diff.level.max<(mean(diff.level.max)-sd(diff.level.max))))
length(which(diff.level.max>(mean(diff.level.max)+2*sd(diff.level.max))|diff.level.max<(mean(diff.level.max)-2*sd(diff.level.max))))
which(diff.level.max>(mean(diff.level.max)+sd(diff.level.max))|diff.level.max<(mean(diff.level.max)-sd(diff.level.max)))