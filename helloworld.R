# Distance from Earth at which the Artemis II 'Hello, World' picture was taken
# www.overfitting.net
# https://www.overfitting.net/2026/04/a-que-distancia-de-la-tierra-se-hizo-la.html


#########################################
# 1. PHOTOGRAPHIC METHOD: FOCAL LENGTH, EARTH DIAMETER RATIO,...

library(Cairo)

    # Radios
    R <- 6371
    d <- 6371 + 9889  # 16260
    
    CairoPNG("orion512.png", width=512, height=800)
    par(cex = 1,
        cex.lab = 1,
        cex.axis = 1,
        cex.main = 1.1)
    
    # Preparar ventana
    plot(0, 0, type = "n",
         xlim = c(-R, d*1.2),
         ylim = c(-R*1.1, R*1.1),
         asp = 1,
         xlab = "km", ylab = "km",
         main = "Situación de la Orion respecto a la Tierra al hacer la foto 'Hello, World'")
    
    # Circunferencia
    theta <- seq(0, 2*pi, length.out = 500)
    lines(R * cos(theta), R * sin(theta), col = "blue", lwd = 2)

    # Punto Orion
    points(d, 0, pch = 16, cex = 1.5)
    
    # Línea centro-Orion
    abline(h=0, v=0)
    
    # ---- RECTAS ±44.482° ----
    angulo1 <- 44.482 * pi / 180
    m1 <- tan(angulo1)
    m2 <- tan(-angulo1)
    
    L <- d * 2
    segments(d - L, 0 - m1*L,
             d + L, 0 + m1*L,
             col = "darkorange", lwd = 1)
    segments(d - L, 0 - m2*L,
             d + L, 0 + m2*L,
             col = "darkorange", lwd = 1)
    
    # ---- RECTAS ±23.068° ----
    angulo2 <- 23.068 * pi / 180
    m3 <- tan(angulo2)
    m4 <- tan(-angulo2)
    
    segments(d - L, 0 - m3*L,
             d + L, 0 + m3*L,
             col = "darkorange", lwd = 1, lty = 2)
    segments(d - L, 0 - m4*L,
             d + L, 0 + m4*L,
             col = "darkorange", lwd = 1, lty = 2)
    
    # Etiqueta
    text(d, 0, labels = "Orion", pos = 3, offset = 1.5)

dev.off()
