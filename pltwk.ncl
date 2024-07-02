; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"

  cdf_path="/glade/derecho/scratch/jpan/archive/b.e23.BMOM.f09_sx0.66av1.aqua.production.0606dlyout/atm/hist/"
  ;cdf_path="/glade/derecho/scratch/zarzycki/archive/QPC5-ne30np4-aquap10-seed3x3/atm/hist/"
  ;cdf_path="/glade/derecho/scratch/jpan/QPC5-ne30np4-aquap10-seed3x3_regrid/"
  ; opening up 6 hourly files 
  ncfiles=systemfunc("ls "+cdf_path+"*h1*.nc")
  ;ncfiles=systemfunc("ls "+cdf_path+"*h1*2002-0*00.nc")
  nfiles=dimsizes(ncfiles)
  f=addfiles(ncfiles,"r")
  print(ncfiles)

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  case     = "0606dlyout"  
  ;case     = "seed3x3"
  diro     = "./wkout"        ; output dir: location of plots 

  latN     = 15
  latS     = -latN       ; symmetric about the equator

  nDayWin  = 96      ; Wheeler-Kiladis [WK] temporal window length (days)
  nDaySkip = -65      ; time (days) between temporal windows [segments]
                     ; negative means there will be overlapping 
                     ; temporal segments

  opt      = True
  opt@debug= True
  opt@pltType     = "png"   ; send graphics to PNG file
  opt@cnLinesOn   = False   ; turn off contour lines
                     ; specify the plot intervals

  spd      = 4       ; samples per day (6 hour sampling rate)
  var      = "OMEGA500"  ; user specified 

; diri     = "/project/cas/shea/WHEELER/"
  diri     = "./"    ; new input directory

  ;fili     = "OMEGA.2001-2005.nc"    


  ;f        = addfile (diri+fili, "r") 
  ;omega    = f->omega(:,0,{latS:latN},:)     ; (time,lat,lon)

  u200 = f[:]->OMEGA500
  print(f[:]->time)

                     ; calculation and create plots
  wkSpaceTime (u200, diro, case, var             \
              ,latN, spd, nDayWin, nDaySkip, opt  )
