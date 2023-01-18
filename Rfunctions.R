path.procorr = '/Users/do/Dropbox/Code/Fortran/ProCorr/ProCorr'

example <- function() {
  file.test = paste0(path.procorr,'/cdm_redshift0')
  compile(clean = T)
  procorr(file.test,'-input 4 -output pdexl')
  delta = read.box(paste0(file.test,'_d.bin'))
  show.field(delta, title = paste0('Density field (',file.test,')'))
  epsil = read.box(paste0(file.test,'_e.bin'))
  show.field(epsil, title = paste0('Epsilon field (',file.test,')'), contrast = 0.7, gamma = 1.3)
  show.powerspectrum(file.test,xlim=c(0.1,5),ylim=c(20,10000))
  show.autocorrelation(file.test)
  show.linecorrelation(file.test)
}
  
compile <- function(clean = T) {
  if (clean) {
    system(paste0('cd ',path.procorr,'; make clean'),intern=T)
  }
  system(paste0('cd ',path.procorr,'; make'),intern=F)
}

procorr <- function(filename,options) {
  system(paste0('cd ',path.procorr,'; ./procorr ',filename,' ',options),intern=FALSE)
}

format <- function(filename) {
  # returns the number of bytes per cell in a 3D-regular grid stored as binary file
  fs = file.info(filename)$size
  n = round((fs/4)^(1/3))
  if (n^3*4==fs) {
    return(list(ncells=n,bytes=4))
  } else {
    n = round((fs/8)^(1/3))
    if (n^3*8==fs) {
      return(list(ncells=n,bytes=8))
    } else {
      stop('File format unknown.')
    }
  }
}

read.box <- function(filename) {
  fm = format(filename)
  dat = readBin(filename,numeric(),size=fm$bytes,n=fm$ncells^3)
  return(array(dat,rep(fm$ncells,3)))
}

show.field <- function(field, title = NULL,
                       dimension = 1, # dimension perpendicular to the screen
                       thickness = c(0.4,0.6), # thickness of slice in coordinates normalized to (0,1)
                       auto.rescale = T,
                       contrast = 1.0, gamma = 0.5, offset = 0.0, show.positive.negative = T,
                       export_filename = NA,
                       show.figure = T) {
  n = dim(field)[1]
  range = seq(max(1,round(thickness[1]*n)),min(n,round(thickness[2]*n)))
  if (dimension==1) {
    imgp = apply(field[range,,],c(2,3),mean)
  } else if (dimension==2) {
    imgp = apply(field[,range,],c(3,1),mean)
  } else if (dimension==3) {
    imgp = apply(field[,,range],c(1,2),mean)
  }
  if (auto.rescale) {
    imgp = (imgp-mean(field))*contrast/2/sd(img)+offset
  } else {
    imgp = imgp*contrast+offset
  }
  
  imgn = -imgp
  imgp[imgp<0] = 0
  imgp[imgp>1] = 1
  imgn[imgn<0] = 0
  imgn[imgn>1] = 1
  imgp = imgp^gamma
  imgn = imgn^gamma
  rgb = array(0,c(n,n,3))
  if (show.positive.negative) {
    rgb[,,2] = imgp
    rgb[,,1] = imgn
  } else {
    for (i in seq(3)) {
      rgb[,,i] = imgp
    }
  }
  if (show.figure) {
    plot(0,0,type='n',xlim=c(0,1),ylim=c(0,1),xlab='',ylab='',xaxs='i',yaxs='i',xaxt='n',yaxt='n',asp=1,bty='n')
    title(title)
    rasterImage(rgb,0,0,1,1)
  }
  if (!is.na(export_filename)) {
    require(png)
    writePNG(rgb,export_filename)
  }
  invisible(rgb)
}

custom.plot <- function(filename, xlim = NULL, ylim = NULL, add = F, col = 'black', yscale = 1, xlab = '', ylab = '', title = '', output = FALSE) {
  dat = read.table(filename,skip=7,sep=',')
  x = dat[,1]
  y = dat[,2]*yscale
  if (is.null(xlim)) {
    xlim = c(0.9*min(x),1.1*max(x))
  }
  if (is.null(ylim)) {
    ylim = c(0.9*min(y),1.1*max(y))
    if (ylim[1]<=0) {ylim[1]=0.01*max(y)}
  }
  if (!add) {plot(0.1,0.1,type='n',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=title,xaxs='i',yaxs='i',log='xy')}
  list = x>0 & y>0
  lines(x[list],y[list],col=col)
  points(x[list],y[list],pch=20,col=col)
  if (output) return(list(x=x,y=y))
}

show.powerspectrum <- function(filename, xlim = NULL, ylim = NULL, add = F, col = 'black', yscale = 1,
                               xlab = 'Wavevector k [(simulation units)'^'-1'~']',
                               ylab = 'Power spectrum p(k)',
                               output = FALSE) {
  out = custom.plot(paste0(filename,'_p.txt'), xlim, ylim, add, col, yscale = yscale, title = 'Power spectrum',
                     xlab=xlab, ylab=ylab, output)
  if (output) return(out)
}

show.autocorrelation <- function(filename, xlim = NULL, ylim = NULL, add = F, col = 'black', yscale = 1,
                                 xlab = 'Correlation scale r [simulation units]',
                                 ylab = '2-point correlation xi(r)',
                                 output = FALSE) {
  out = custom.plot(paste0(filename,'_x.txt'), xlim, ylim, add, col, yscale = yscale, title = '2-point correlation',
                     xlab=xlab, ylab=ylab, output)
  if (output) return(out)
}

show.linecorrelation <- function(filename, xlim = NULL, ylim = NULL, add = F, col = 'black', yscale = 1,
                                 xlab = 'Correlation scale r [simulation units]',
                                 ylab = 'Line correlation l(r)',
                                 output = FALSE) {
  out = custom.plot(paste0(filename,'_l.txt'), xlim, ylim, add, col, yscale = yscale, title = 'Line correlation',
                     xlab=xlab, ylab=ylab, output)
  if (output) return(out)
}