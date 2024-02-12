pro reconrates
  ;index 0 is positive rate, index 1 is negative rate
  stdrecs=make_array(481,2)
  meanrecs=make_array(481,2)
  medianrecs=make_array(481,2)
  maximprecs=make_array(481,2)
  maxrecrates=make_array(481,2)
  recflx=make_array(481,2)
  file = '/Users/coletamburri/Desktop/Impulsiveness_Paper/imp_dev/best_r2corr.txt'
  OPENR, lun, file, /GET_LUN
  ; Read one line at a time, saving the result into array
  flaredirectories = ''
  line = ''
  WHILE NOT EOF(lun) DO BEGIN & $
    READF, lun, line & $
    flaredirectories= [flaredirectories, line] & $
  ENDWHILE
  
  file2 = '/Users/coletamburri/Desktop/Impulsiveness_Paper/imp_dev/best_r2corr_inds.txt'
  OPENR, lun2, file2, /GET_LUN
  ; Read one line at a time, saving the result into array
  flareinds = ''
  line2 = ''
  WHILE NOT EOF(lun2) DO BEGIN & $
    READF, lun2, line2 & $
    flareinds= [flareinds, line2] & $
  ENDWHILE
; Close the file and free the file unit
FREE_LUN, lun
FREE_LUN, lun2
for i= 1,481 do begin
  print,i
  avgposrec = []
  avgnegrec = []
  stdposrec = []
  stdnegrec = []
  medianposrec = []
  mediannegrec = []

  flarename =flaredirectories[i]
  iflare = flareinds[i]
  plot_rbn4cole,cut=8,maxposrecrate=maxposrecrate,maxnegrecrate=maxnegrecrate,maximpposrec=maximpposrec,maximpnegrec=maximpnegrec,avgposrec=avgposrec,avgnegrec=avgnegrec,stdposrec=stdposrec,stdnegrec=stdnegrec,medianposrec=medianposrec,mediannegrec=mediannegrec,flaredir0=flarename,flx=flx,dir='/Users/coletamburri/Desktop/Old Mac/CU_Research/REC_RATES/recfiles_rad/',iflare=iflare
  stdrecs[i-1,0]=stdposrec
  stdrecs[i-1,1]=stdnegrec
  meanrecs[i-1,0]=avgposrec
  meanrecs[i-1,1]=avgnegrec
  medianrecs[i-1,0]=medianposrec
  medianrecs[i-1,1]=mediannegrec
  maximprecs[i-1,0]=maximpposrec
  maximprecs[i-1,1]=maximpnegrec
  maxrecrates[i-1,0]=maxposrecrate
  maxrecrates[i-1,1]=maxnegrecrate
  recflx[i-1,0]=max(flx[*,2])
  recflx[i-1,1]=min(flx[*,0])
  
endfor
SAVE,maximprecs,maxrecrates,stdrecs,meanrecs,medianrecs,flaredirectories,flareinds,recflx, FILENAME = '/Users/coletamburri/Desktop/Impulsiveness_Paper/imp_dev/recratesidl_r2corr.sav'

end