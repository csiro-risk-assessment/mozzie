library(data.table)
species = c("Ac", "Ag")
genotypes = c("ww", "wc", "wr", "cc", "cr", "rr")
sexes = c("male", "female")
years = 2022:2023
monthdays = c(31,28,31,30,31,30,31,31,30,31,30,31) # No leap years

mozzies = array(NA, c(2,6,2,2,12,31,100,100))


for (s in 1:2)
{
  S = species[s]
  print(S)
  for (g in 1:6)
  {
    G = genotypes[g]
    print(G)
    {
      for (e in 1:2)
      {
        E = sexes[e]
        print(E)
        {
          for (y in 1:2)
          {
            Y = years[y]
            print(Y)
            for (m in 1:12)
            {
              for (d in 1:monthdays[m])
              {
                X = fread(sprintf("%s_%s_%s_%d_%02d_%02d.csv",
                                  S, G, E, Y, m, d), 
                          header = FALSE, skip = 2, sep = ',')
                mozzies[s,g,e,y,m,d,,] = as.matrix(X)
              }
            }
          }
        }
      }
    }
  }
}

save(mozzies, file = "mozzies.Rdata")

library(png)

min.Ac = min(
  apply(mozzies[1,,,,,,,], 3:7, sum), 
  na.rm = TRUE)
max.Ac = max(
  apply(mozzies[1,,,,,,,], 3:7, sum), 
  na.rm = TRUE)
min.Ag = min(
  apply(mozzies[2,,,,,,,], 3:7, sum), 
  na.rm = TRUE)
max.Ag = max(
  apply(mozzies[2,,,,,,,], 3:7, sum), 
  na.rm = TRUE)
max.notww = max(
  apply(mozzies[,-1,,,,,,], 4:8, sum), 
  na.rm = TRUE)


n = 1
for (y in 1:2)
{
  Y = years[y]
  print(Y)
  for (m in 1:12)
  {
    print(m)
    for (d in 1:monthdays[m])
    {
      red = apply(mozzies[1,,,y,m,d,,], 3:4, sum)
      red = (red - min.Ac)/(max.Ac - min.Ac)
      green = apply(mozzies[2,,,y,m,d,,], 3:4, sum)
      green = (green - min.Ag)/(max.Ag - min.Ag)
      blue = apply(mozzies[,-1,,y,m,d,,], 4:5, sum)
      blue = blue / max.notww
      image = array(NA, c(100,100,3))
      image[,,1] = red
      image[,,2] = green
      image[,,3] = blue
      writePNG(image, target = paste0(n, ".png"))
      n = n + 1
    }
  }
}



library(magick)
img_list = lapply(paste0(1:730, ".png"), image_read)
## join the images together
img_joined <- image_join(img_list)
## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 20)
## view animated image
#img_animated
## save to disk
image_write(image = img_animated,
            path = "anim.gif")
#for (i in 1:730) file.remove(paste0(i,".png"))

gm = apply(mozzies[,-1,,,,,,], 4:8, sum)
gm = array(
  aperm(gm, c(3:1, 4, 5)),
        c(2*12*31, 100, 100))
gm = apply(gm, 2:3, function(x) min(which(x > 1), na.rm = TRUE))
gm[gm == Inf] = NA
library(raster)

africa.albers = '+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 
+units=m +no_defs'
gm = raster(as.matrix(gm)[nrow(gm):1,], 
       xmn=-1799134.0,ymn=-1299134.0,xmx=-1799134.0+5000*100,ymx=-1299134.0+5000*100, 
       crs = africa.albers)

pdf("figure.pdf", 7, 7)
plot(gm)
contour(gm, add = TRUE, drawlabels = FALSE, lwd = 1)
dev.off()

