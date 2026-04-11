# False 3D animation of a whole Moon day
# www.overfitting.net
# https://www.overfitting.net/2026/04/un-dia-en-la-luna.html

library(tiff)  # save 16-bit TIFF's
library(Rcpp)


# Equirectangular to orthographic projection conversion
# 3. Two radius + bilinear interpolation
cppFunction('
    Rcpp::NumericVector equirect_to_orthographic_V2_bilinear(
            Rcpp::NumericMatrix imgR,
            Rcpp::NumericMatrix imgG,
            Rcpp::NumericMatrix imgB,
            int outW,
            int outH,
            double Rx,
            double Ry)
    {
        int texH = imgR.nrow();
        int texW = imgR.ncol();
    
        Rcpp::NumericVector out(outH * outW * 3);
        double cx = (outW - 1) * 0.5;
        double cy = (outH - 1) * 0.5;
    
        for (int j = 0; j < outH; j++) {
            for (int i = 0; i < outW; i++) {
                int base = j + i * outH;
    
                double X = (i - cx) / Rx;
                double Y = (cy - j) / Ry;
    
                if (X*X + Y*Y > 1.0) {
                    out[base] = 0.0;
                    out[base + outH*outW] = 0.0;
                    out[base + 2*outH*outW] = 0.0;
                    continue;
                }
    
                double Z = sqrt(std::max(0.0, 1.0 - X*X - Y*Y));
    
                double lat = asin(Y);
                double lon = atan2(X, Z);
    
                double u = (lon + M_PI/2) / M_PI;
                double v = (M_PI/2 - lat) / M_PI;
    
                double px = u * (texW - 1);
                double py = v * (texH - 1);
    
                // Bilinear interpolation indices
                int x0 = std::floor(px);
                int y0 = std::floor(py);
                int x1 = std::min(x0 + 1, texW - 1);
                int y1 = std::min(y0 + 1, texH - 1);
    
                double dx = px - x0;
                double dy = py - y0;
    
                // Interpolation weights
                double w00 = (1 - dx) * (1 - dy);
                double w10 = dx * (1 - dy);
                double w01 = (1 - dx) * dy;
                double w11 = dx * dy;
    
                // R channel
                double R00 = imgR(y0, x0);
                double R10 = imgR(y0, x1);
                double R01 = imgR(y1, x0);
                double R11 = imgR(y1, x1);
                out[base] =
                    R00 * w00 + R10 * w10 + R01 * w01 + R11 * w11;
    
                // G channel
                double G00 = imgG(y0, x0);
                double G10 = imgG(y0, x1);
                double G01 = imgG(y1, x0);
                double G11 = imgG(y1, x1);
                out[base + outH*outW] =
                    G00 * w00 + G10 * w10 + G01 * w01 + G11 * w11;
    
                // B channel
                double B00 = imgB(y0, x0);
                double B10 = imgB(y0, x1);
                double B01 = imgB(y1, x0);
                double B11 = imgB(y1, x1);
                out[base + 2*outH*outW] =
                    B00 * w00 + B10 * w10 + B01 * w01 + B11 * w11;
            }
        }
        return out;
    }
')


# Hillshade calculation
hillshademap=function(DEM, dx=25, dlight=c(0, 2, 3), gamma=1) {
    # hillshademap() inputs DEM data and outputs a hillshade matrix
    #
    # DEM: digital elevation map matrix
    # dx: DEM resolution/cell size (same units as elevation values)
    # dlight: lighting direction. It can be defined in three ways:
    #   a) 1D value indicating light source Azimuth in degrees (0-360)
    #      (0=North, 90=East, 180=South, 270=West)
    #   b) 2D vector indicating light source (X,Y) coordinates
    #   c) 3D vector indicating light source (X,Y,Z) coordinates:
    #      (X=South, Y=East, Z=Up)
    #      dlight=c(0, 2, 3)  # sunrise
    #      dlight=c(0, 0, 1)  # midday
    #      dlight=c(0,-2, 3)  # sunset
    #   NOTE: both in a) and b) a 45º Elevation angle is applied
    # gamma: optional output gamma lift
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    # If array turn its first dimension into a genuine matrix
    if (!is.matrix(DEM)) {
        print("WARNING: input DEM is not a matrix but an array. First dimension is used")
        DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    }
    
    # Deal with lighting direction
    if (length(dlight)==1) dlight=c(-cos(dlight*pi/180), sin(dlight*pi/180))
    if (length(dlight)==2) dlight=c(dlight, (dlight[1]^2+dlight[2]^2)^0.5)
    dlightM=sum(dlight^2)^0.5
    
    # Vectorial product to calculate n (orthogonal vector)
    nx = 2*dx*(DEM[1:(DIMY-2), 2:(DIMX-1)] - DEM[3:DIMY,     2:(DIMX-1)])
    ny = 2*dx*(DEM[2:(DIMY-1), 1:(DIMX-2)] - DEM[2:(DIMY-1), 3:DIMX])
    nz = 4*dx^2
    nM = (nx^2 + ny^2 + nz^2)^0.5
    
    # Dot product to calculate cos(theta)
    dn = dlight[1]*nx + dlight[2]*ny + dlight[3]*nz  # (DIMY-2)x(DIMX-2) matrix
    
    # Reflectance (=cos(theta))
    hillshadepre=dn/(dlightM*nM)
    hillshadepre[hillshadepre<0]=0  # clip negative values
    
    # Add 1-pix 'lost' borders
    hillshademap=matrix(0, nrow=DIMY, ncol=DIMX)
    hillshademap[2:(DIMY-1), 2:(DIMX-1)]=hillshadepre
    rm(hillshadepre)
    hillshademap[c(1,DIMY),]=hillshademap[c(2,DIMY-1),]
    hillshademap[,c(1,DIMX)]=hillshademap[,c(2,DIMX-1)]
    
    return(hillshademap^(1/gamma))
}


###########################################################
# MOON DAY


# Process DEM to obtain HILLSHADE

DEMWHOLE=readTIFF("ldem_16_uint.tif")
hist(DEMWHOLE, breaks=500, xlim=c(0,1))
DEMWHOLE=DEMWHOLE-min(DEMWHOLE)
DEMWHOLE=DEMWHOLE/max(DEMWHOLE)

writeTIFF(DEMWHOLE, "ldem_16_uint_norm.tif")
# MAXHEIGHT=19.5*1000  # 19.5km between Moon's deepest basin and highest point

# dlight(X=South, Y=East, Z=Up)
hillshade=hillshademap(DEMWHOLE, dx=1/10, dlight=c(0, 1, 1))

# Save hillshade
writeTIFF(hillshade, "hillshade.tif", bits.per.sample=16, compression="LZW")


# Build "moon_final_extended.tif" in Photoshop as a combination of colour and hillshade

# Build Full HD animation frames by reprojecting sections of "moon_final_extended.tif"
# from equirectangular to ortographic
ANCHO=1000
ALTO=ANCHO
img=readTIFF("moon_final_extended.tif")
width=nrow(img)
mask=readTIFF("mask.tif")  # solar mask (dark area)
background=readTIFF("backgroundfinal.tif")  # full HD title and label
alpha=0.05  # weighted mix

NFRAMES=640  # total number of frames
INC=9  # 9px offset from frame to frame
for (f in 1:NFRAMES) {
    print(paste0("Processing frame ", f, "/", NFRAMES, "..."))
    
    # Crop
    DEM=img[, ((f-1)*INC+1):((f-1)*INC+width), ]
    # Sun mask
    DEM=DEM*(alpha+(1-alpha)*replicate(3, mask))
    
    # Projection conversion
    # imgR, imgG, imgB: full planet texture from -90..+90 longitude and latitude
    out <- equirect_to_orthographic_V2_bilinear(imgR=DEM[,,1], imgG=DEM[,,2], imgB=DEM[,,3],
                                                outW=ANCHO, outH=ALTO, 
                                                Rx=ANCHO/2, Ry=ANCHO/2)
    dim(out) <- c(ALTO, ANCHO, 3)
    out2=background
    out2[(1080/2-ANCHO/2+1):(1080/2+ANCHO/2), (1920/2-ANCHO/2+1):(1920/2+ANCHO/2),]=out
    writeTIFF(out2, sprintf("moon_%05d.tif", f), bits.per.sample=8)
}



# MP4 Video (MPEG-4 AVC/H.264):
# ffmpeg -framerate 24 -i moon_%05d.tif -i 26seconds.wav -c:v libx264 -crf 18 -pix_fmt yuv420p -vf reverse moonday.mp4