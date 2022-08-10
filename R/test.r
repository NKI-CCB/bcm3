

setwd("D:/git/bcm3/bin/Release")
dyn.load("D:/git/bcm3/bin/Release/bcmrbridge.dll")


a <- "test"
b <- "longbufferlongbuffer"
res <- .C("testfun", a, b, 10, 5.0)

res2 <- .C("testfun2", res[[2]], c)

dyn.unload("D:/git/bcm3/bin/Release/bcmrbridge.dll")
