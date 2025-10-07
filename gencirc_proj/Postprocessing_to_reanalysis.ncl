;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                 ;;
;; Postproccessing.ncl                                             ;;
;;                                                                 ;;
;; AUTHOR                                                          ;;
;; -> Allen Mewhinney ;; Based on Mingyu's Postproccessing script  ;;
;;                                                                 ;;
;; DESCRIPTION                                                     ;;
;; -> 2021.01.15 Process GFDL model output experiments             ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load constants.ncl and energetics.ncl ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; load "/glade/work/jpan/597_gencirc_gfdl/MODEL/constants.ncl"
; load "/glade/work/jpan/597_gencirc_gfdl/MODEL/energetics.ncl"
load "./constants.ncl"
load "./energetics.ncl"

begin  ; begin main script

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; Define script specific variables  ;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	ti           = 0                ; initial time step
	tf           = 19               ; Final time step
	levout       = LEV37            ; output pressure levels
	nlev         = dimsizes(levout)
	nlat         = 145  ; convert to this number of latitudes
	nlon         = 288 ; convert to this number of longitudes
	trunc        = 42  ; truncate with at this wavenumber
	CenterLon    = 0

	model_run    = "Pacific_trend_central_pacific_nodiab"
	climate      = "Pacific_trend_central_pacific_diabatic"
	
	mdl_path     = "/glade/work/jpan/597_gencirc_gfdl/ModelCopy/out/" + model_run + "/ensemble0/atmos_daily.nc"
    heating_path = "/glade/work/jpan/597_gencirc_gfdl/ModelCopy/out/" + model_run + "/ensemble0/convective_f.nc"
	clm_path     = "/glade/work/jpan/597_gencirc_gfdl/ModelCopy/out/" + climate   + "/ensemble0/atmos_daily.nc"
	shum_clm_path = "/glade/work/jpan/597_gencirc_gfdl/ModelCopy/input/ERA5.Annual.1980to2020.T42.nc"
	outpath      = "/glade/work/jpan/597_gencirc_gfdl/model_results/postprocessed_" + model_run + ".nc"
	
	; read in model data
	mdl_file     = addfile(mdl_path,"r")
    heating_file = addfile(heating_path,"r")
	clm_file     = addfile(clm_path,"r")
    shum_clm_file = addfile(shum_clm_path,"r")

	print("Loading data.")
	time         = mdl_file->time(ti:tf)
	ucomp_mdl    = mdl_file->ucomp(ti:tf,:,:,:)
	vcomp_mdl    = mdl_file->vcomp(ti:tf,:,:,:)
	temp_mdl     = mdl_file->temp(ti:tf,:,:,:)
	height_mdl   = mdl_file->height(ti:tf,:,:,:)
	ps_mdl       = mdl_file->ps(ti:tf,:,:)
	tracer_mdl   = mdl_file->tracer1(ti:tf,:,:,:)
	omega_mdl    = mdl_file->omega(ti:tf,:,:,:)

    heating      = dble2flt(heating_file->dteq1(0,:,:,:))
    copy_VarCoords(ucomp_mdl(0,:,:,:),heating)

	ucomp_clm    = clm_file->ucomp(ti:tf,:,:,:)
	vcomp_clm    = clm_file->vcomp(ti:tf,:,:,:)
	temp_clm     = clm_file->temp(ti:tf,:,:,:)
	height_clm   = clm_file->height(ti:tf,:,:,:)
	ps_clm       = clm_file->ps(ti:tf,:,:)
	tracer_clm   = shum_clm_file->shum(:,:,:)
	omega_clm    = clm_file->omega(ti:tf,:,:,:)

    
	print("Converting data to a fixed grid and horizontally interpolating")
	ucompF_mdl   = g2fsh_Wrap(ucomp_mdl , (/nlat,nlon/))
	vcompF_mdl   = g2fsh_Wrap(vcomp_mdl , (/nlat,nlon/))
	tempF_mdl    = g2fsh_Wrap(temp_mdl  , (/nlat,nlon/))
	heightF_mdl  = g2fsh_Wrap(height_mdl, (/nlat,nlon/))
	psF_mdl      = g2fsh_Wrap(ps_mdl    , (/nlat,nlon/))
	tracerF_mdl  = g2fsh_Wrap(tracer_mdl, (/nlat,nlon/))
	omegaF_mdl   = g2fsh_Wrap(omega_mdl , (/nlat,nlon/))

    heatingF     = g2fsh_Wrap(heating, (/nlat,nlon/))

	ucompF_clm   = g2fsh_Wrap(ucomp_clm , (/nlat,nlon/))
	vcompF_clm   = g2fsh_Wrap(vcomp_clm , (/nlat,nlon/))
	tempF_clm    = g2fsh_Wrap(temp_clm  , (/nlat,nlon/))
	heightF_clm  = g2fsh_Wrap(height_clm, (/nlat,nlon/))
	psF_clm      = g2fsh_Wrap(ps_clm    , (/nlat,nlon/))
	tracerF_clm  = g2fsh_Wrap(tracer_clm, (/nlat,nlon/))
	omegaF_clm   = g2fsh_Wrap(omega_clm , (/nlat,nlon/))

	latitude     = dble2flt(ucompF_mdl&lat)
	longitude    = dble2flt(ucompF_mdl&lon)

	; Hybrid coordinate definitions
	hyB          = mdl_file->pfull / 1000.
	nhyB         = dimsizes(hyB)
	hyA          = hyB
	hyA          = 0.0 ; Change all values to 0 in hyA
	tempF_bot    = tempF_mdl(:,nhyB-1,:,:) ; lowest level temperature
	ZF_bot       = G * heightF_mdl(:,nhyB-1,:,:) ; lowest level geopotential

	print("Vertically interpolating to ECMWF pressure levels.")
	ucompFp_mdl  = vinth2p_ecmwf(ucompF_mdl , hyA, hyB, levout, psF_mdl, 1, 1000., 1, True, 0, tempF_bot, ZF_bot)
	vcompFp_mdl  = vinth2p_ecmwf(vcompF_mdl , hyA, hyB, levout, psF_mdl, 1, 1000., 1, True, 0, tempF_bot, ZF_bot)
	tempFp_mdl   = vinth2p_ecmwf(tempF_mdl  , hyA, hyB, levout, psF_mdl, 1, 1000., 1, True, 0, tempF_bot, ZF_bot)
	heightFp_mdl = vinth2p_ecmwf(heightF_mdl, hyA, hyB, levout, psF_mdl, 1, 1000., 1, True, 0, tempF_bot, ZF_bot)
	tracerFp_mdl = vinth2p_ecmwf(tracerF_mdl, hyA, hyB, levout, psF_mdl, 1, 1000., 1, True, 0, tempF_bot, ZF_bot)
	omegaFp_mdl  = vinth2p_ecmwf(omegaF_mdl , hyA, hyB, levout, psF_mdl, 1, 1000., 1, True, 0, tempF_bot, ZF_bot)

    heatingFp = vinth2p_ecmwf(heatingF, hyA, hyB, levout, psF_mdl(0,:,:),1, 1000., 1, True, 0, tempF_bot(0,:,:), ZF_bot(0,:,:))
    
	; Do same for climate file
	hyB          = clm_file->pfull / 1000.
	nhyB         = dimsizes(hyB)
	hyA          = hyB
	hyA          = 0.0 ; Change all values to 0 in hyA
	tempF_bot    = tempF_clm(:,nhyB-1,:,:) ; lowest level temperature
	ZF_bot       = G * heightF_clm(:,nhyB-1,:,:) ; lowest level geopotential

	print("Vertically interpolating to ECMWF pressure levels.")
	ucompFp_clm  = vinth2p_ecmwf(ucompF_clm , hyA, hyB, levout, psF_clm, 1, 1000., 1, True, 0, tempF_bot, ZF_bot)
	vcompFp_clm  = vinth2p_ecmwf(vcompF_clm , hyA, hyB, levout, psF_clm, 1, 1000., 1, True, 0, tempF_bot, ZF_bot)
	tempFp_clm   = vinth2p_ecmwf(tempF_clm  , hyA, hyB, levout, psF_clm, 1, 1000., 1, True, 0, tempF_bot, ZF_bot)
	heightFp_clm = vinth2p_ecmwf(heightF_clm, hyA, hyB, levout, psF_clm, 1, 1000., 1, True, 0, tempF_bot, ZF_bot)
	tracerFp_clm = vinth2p_ecmwf(tracerF_clm, hyA, hyB, levout, psF_clm(0,:,:), 1, 1000., 1, True, 0, tempF_bot(0,:,:), ZF_bot(0,:,:))
	omegaFp_clm  = vinth2p_ecmwf(omegaF_clm , hyA, hyB, levout, psF_clm, 1, 1000., 1, True, 0, tempF_bot, ZF_bot)
	
	; Rename coordinates
	ucompFp_mdl!1         = "level"
	ucompFp_mdl&level  = levout
	ucompFp_mdl!2         = "latitude"
	ucompFp_mdl&latitude  = latitude
	ucompFp_mdl!3         = "longitude"
	ucompFp_mdl&longitude = longitude

    ; Rename coordinates for trace clm
	tracerFp_clm!0         = "level"
	tracerFp_clm&level  = levout
	tracerFp_clm!1         = "latitude"
	tracerFp_clm&latitude  = latitude
	tracerFp_clm!2         = "longitude"
	tracerFp_clm&longitude = longitude

	; use copr_VarCoords to update the other variables too
	copy_VarCoords(ucompFp_mdl, vcompFp_mdl) 
	copy_VarCoords(ucompFp_mdl, tempFp_mdl)
	copy_VarCoords(ucompFp_mdl, heightFp_mdl)
	copy_VarCoords(ucompFp_mdl, tracerFp_mdl)
	copy_VarCoords(ucompFp_mdl, omegaFp_mdl)
	
    ; Rename coordinates for heating
	heatingFp!0         = "level"
	heatingFp&level  = levout
	heatingFp!1         = "latitude"
	heatingFp&latitude  = latitude
	heatingFp!2         = "longitude"
	heatingFp&longitude = longitude

    ; Compute differences 
    ucompFp_anom = ucompFp_mdl - ucompFp_clm
    vcompFp_anom = vcompFp_mdl - vcompFp_clm
    tempFp_anom = tempFp_mdl - tempFp_clm
    heightFp_anom = heightFp_mdl - heightFp_clm
    tracerFp_anom = tracerFp_mdl - conform(tracerFp_mdl,tracerFp_clm,(/1,2,3/))
    omegaFp_anom = omegaFp_mdl - omegaFp_clm

	copy_VarCoords(ucompFp_mdl, ucompFp_anom)
	copy_VarCoords(ucompFp_mdl, vcompFp_anom)
	copy_VarCoords(ucompFp_mdl, tempFp_anom)
	copy_VarCoords(ucompFp_mdl, heightFp_anom)
	copy_VarCoords(ucompFp_mdl, tracerFp_anom)
	copy_VarCoords(ucompFp_mdl, omegaFp_anom)

    printVarSummary(ucompFp_anom)
    
	; Find Streamfuntion
	sfvp_mdl              = uv2sfvpF(ucompFp_mdl,vcompFp_mdl) 
	psi_mdl               = sfvp_mdl(0,:,:,:,:)
	psi_mdl!0             = "time"
	psi_mdl&time          = ispan(ti+1,tf+1,1)
	psi_mdl!1             = "level"
	psi_mdl&level         = levout
	psi_mdl!2             = "latitude"
	psi_mdl&latitude      = latitude
	psi_mdl!3             = "longitude"
	psi_mdl&longitude     = longitude
	psi_mdl@units         = "m**2 s**-1" 
	; find climate Streamfunction
	sfvp_clm              = uv2sfvpF(ucompFp_clm,vcompFp_clm)
	psi_clm               = sfvp_clm(0,:,:,:,:)
	copy_VarCoords(psi_mdl,psi_clm)
	; Find Streamfunction anomaly
	psi_anom              = psi_mdl - psi_clm
	copy_VarCoords(psi_mdl,psi_anom)
	
	
	; Save data to file
	print("Saving " + outpath + ".")
	system("rm " + outpath)                 ; remove file if present
	outfile = addfile(outpath, "c")	        ; create file

	; create global attributes
	foutAtt                  = True ; assign file attributes
	foutAtt@path             = outpath
	foutAtt@filename         = model_run + ".nc"
	foutAtt@title            = "GFDL Model Experiment: " + model_run
	foutAtt@source_file      = mdl_file
	foutAtt@creation_date    = systemfunc("date")
	foutAtt@Convections      = "None"	
	;fileattdef(outfile, foutAtt) ; copy file attributes to file

	;Save data values to file
	outfile->u_mdl  = ucompFp_mdl
	outfile->v_mdl  = vcompFp_mdl
	outfile->temp_mdl   = tempFp_mdl
	outfile->height_mdl = heightFp_mdl
	outfile->tracer_mdl = tracerFp_mdl
	outfile->omega_mdl  = omegaFp_mdl
	outfile->psi_mdl    = psi_mdl
	outfile->u_anom  = ucompFp_anom
	outfile->v_anom  = vcompFp_anom
	outfile->temp_anom   = tempFp_anom
	outfile->height_anom = heightFp_anom
	outfile->tracer_anom = tracerFp_anom
	outfile->omega_anom  = omegaFp_anom
	outfile->psi_anom    = psi_anom
    outfile->heating     = heatingFp

end ; end main script
