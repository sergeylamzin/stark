    # Copyright (C) 2014  Sergey Lamzin, https://github.com/sergeylamzin/stark

    # This file is part of the StarK genome assembler.

    # StarK is free software: you can redistribute it and/or modify it
    # under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 3 of the License, or
    # (at your option) any later version.

    # StarK is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.

    # You should have received a copy of the GNU General Public License
    # along with Foobar.  If not, see <http://www.gnu.org/licenses/>.


#data in matrix format provided by variable z
#x <- x-axis
#y <- y axis

cmap <- function(nbcol, alpha=1) {
jet.colors <- colorRampPalette( c("white", "green", "yellow", "blue", "red") )
jet.colors(nbcol)
}

# x = seq(1, nrow(z))
# y = seq(1, ncol(z))

ccl <- c(1, 2, 5)
covlab <- c(0, ccl, 10*ccl, 100*ccl, 1000*ccl, 10000*ccl, 100000*ccl, 1000000*ccl, 10000000*ccl, 100000000*ccl, 1000000000*ccl)

filled.contour(
x
,log(y)
,z
,color = cmap 
,plot.axes = { axis(1,x); axis(2,log(y), labels=y)}
,xlab="k"
,ylab="coverage"
,main="k* Coverage Histogram"
, key.axes={axis(4,at=log(covlab, 10), labels=covlab)}
,cex=1.5
,cex.lab=1.5
,cex.axis=1.5
,cex.main=1.5
,cex.sub=1.5
)


