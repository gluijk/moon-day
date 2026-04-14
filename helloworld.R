# Distance from Earth at which the Artemis II 'Hello, World' picture was taken
# www.overfitting.net
# https://www.overfitting.net/2026/04/a-que-distancia-de-la-tierra-se-hizo-la.html


############################################
# 1. DRAW ESTIMATED DISTANCE TO EARTH

library(Cairo)


# Radius
R <- 6371
d <- 6371 + 9889  # 16260
    
CairoPNG("orion512.png", width=512, height=800)
    par(cex = 1,
        cex.lab = 1,
        cex.axis = 1,
        cex.main = 1.1)
    
    # Prepare plot
    plot(0, 0, type = "n",
         xlim = c(-R, d*1.2),
         ylim = c(-R*1.1, R*1.1),
         asp = 1,
         xlab = "km", ylab = "km",
         main = "Situación de la Orion respecto a la Tierra al hacer la foto 'Hello, World'")
    
    # Circumference
    theta <- seq(0, 2*pi, length.out = 500)
    lines(R * cos(theta), R * sin(theta), col = "blue", lwd = 2)
    
    # Centre-Orion line
    abline(h=0, v=0)
    
    L <- d * 2
    # ---- LINES ±44.482° ----
    angulo1 <- 44.482 * pi / 180
    m1 <- tan(angulo1)
    m2 <- tan(-angulo1)
    
    segments(d - L, 0 - m1*L,
             d + L, 0 + m1*L,
             col = "darkorange", lwd = 1)
    segments(d - L, 0 - m2*L,
             d + L, 0 + m2*L,
             col = "darkorange", lwd = 1)
    
    # ---- LINES ±23.068° ----
    angulo2 <- 23.068 * pi / 180
    m3 <- tan(angulo2)
    m4 <- tan(-angulo2)
    
    segments(d - L, 0 - m3*L,
             d + L, 0 + m3*L,
             col = "darkorange", lwd = 1, lty = 2)
    segments(d - L, 0 - m4*L,
             d + L, 0 + m4*L,
             col = "darkorange", lwd = 1, lty = 2)
    
    # Orion
    points(d, 0, pch = 16, cex = 1.5)
    text(d, 0, labels = "Orion", pos = 3, offset = 1.5)

dev.off()



############################################
# 2. CONFIRM ESTIMATED DISTANCE OVER EARTH MAP

library(ggmap)  # map_data()
library(data.table)


# BITMAP DRAWING FUNCTIONS

NewBitmap = function(dimx, dimy, val=0) {
    # Crea bitmap de dimensiones dimx y dimy
    return(array(val,c(dimx,dimy)))
}

DrawCircle = function(img, x0, y0, r, inc=TRUE, val=1, fill=FALSE, thick=1) {
    # Dibuja círculo de centro (x0,y0) y radio r
    # Por defecto método no destructivo, con valor=1 y sin relleno
    # Puede elegirse el grosor si no se rellena
    img=DrawEllip(img, x0, y0, r, r, inc, val, fill, thick)
    
    return(img)
}

DrawEllip = function(img, x0, y0, a, b, inc=TRUE, val=1, fill=FALSE, thick=1) {
    # Dibuja elipse de centro (x0,y0) y radios a y b
    # Por defecto método no destructivo, con valor=1 y sin relleno
    # Puede elegirse el grosor si no se rellena
    # Aquí no redondeamos para tener más precisión en la división
    if (fill) {
        indices=which( ((row(img)-x0)/a)^2 + ((col(img)-y0)/b)^2 < 1 )
    } else {
        indices=which( ((row(img)-x0)/(a+thick/2))^2 + ((col(img)-y0)/(b+thick/2))^2 <  1 &
                           ((row(img)-x0)/(a-thick/2))^2 + ((col(img)-y0)/(b-thick/2))^2 >= 1 )
    }
    if (inc) img[indices]=img[indices]+val
    else img[indices]=val
    
    return(img)
}

SaveBitmap = function(img, name, trunc=TRUE, gamma=1) {
    # Guarda bitmap en formato PNG
    # Solo si trunc=FALSE y la imagen excede de 1 se reescala a 1
    library(png)
    img[img<0]=0
    if (trunc) img[img>1]=1
    if (tolower(substr(name, nchar(name)-3, nchar(name))) != ".png") name=paste0(name,".png")
    writePNG(t(img[,ncol(img):1] / max(max(img),1))^(1/gamma), name)
}


# COORDINATES CONVERSION

polar2x = function(R, phi, theta) {R*cos(theta)*sin(phi)}
polar2y = function(R, phi, theta) {R*sin(theta)}
polar2z = function(R, phi, theta) {-R*cos(theta)*cos(phi)}


# DRAW EARTH

# Physical parameters
Rearth=6371  # Earth average radius (km)
Oriond=9889  # Orion altitude (km)

# Obtain and process coordinates
DT=data.table(map_data("world"))  # long/lat pairs for all countries
DT=DT[, .(num=.N), by=.(long, lat)]  # summarize to deduplicate points

# deg to rad conversion
DT$phi=DT$long*pi/180  # longitude
DT$theta=DT$lat*pi/180  # latitude

# polar to XYZ coordinates conversion
# NOTE: x and y don't depend on d, but z depends on d so z is nested
DT$x=polar2x(Rearth, DT$phi, DT$theta)
DT$y=polar2y(Rearth, DT$phi, DT$theta)

IMAGESIZE=512
TH=0.9  # allow border

d=Oriond
DT$z=polar2z(Rearth, DT$phi, DT$theta) + Rearth+d  # Earth along Z axis
# Distance from each map point to observation point (0,0,0)
DT$dist=(DT$x^2+DT$y^2+DT$z^2)^0.5

distmax=( (d+Rearth)^2 - Rearth^2 )^0.5  # max distance to visible points
DTplot=DT[DT$dist<=distmax]  # keep only visible points

# Draw globe map
img=NewBitmap(IMAGESIZE, IMAGESIZE)
NCOLDIV2=ncol(img)/2
NROWDIV2=nrow(img)/2

# Calculate focal length to fit Earth in the final image and FOV
thetamax=acos(Rearth/(d+Rearth))
f=min(NCOLDIV2,NROWDIV2)*
    (d+Rearth-Rearth*cos(thetamax))/(Rearth*sin(thetamax))*TH
FOV=(pi-2*thetamax)*180/pi  # FOV in deg
print(paste0("distance d=", round(d), "km, FOV=",
             round(FOV,1), "º, ", nrow(DTplot), " points"))

img=DrawCircle(img, NCOLDIV2, NROWDIV2, min(NCOLDIV2,NROWDIV2)*TH,
               fill=TRUE, val=0.25)
DTplot$factor=f/DTplot$z
DTplot$xp=DTplot$x*DTplot$factor + NCOLDIV2
DTplot$yp=DTplot$y*DTplot$factor + NROWDIV2
img[round(cbind(DTplot$xp, DTplot$yp))]=1  # draw points

SaveBitmap(img, paste0("earthfromorion"))

