;pro moses2_RichardsonLucyDeconvol

  ;Load MOSES2 level one dataset into memory
  restore, 'IDLWorkspace/MOSES2_LevelOne_Dataset/mosesDarkSub.sav'
  
  ; Select single image to deconvolve
  img = cube_zero[*,*,17]
  
  iterations = 5   ; total number of blind interations to perform
  m = 5  ; Number of Richardson-Lucy interations per blind iteration
  
  ; Construct our guess for the PSF
;  winX = 40   ; Size of the Hanning window in x
;  winY = 40   ; Size of Hanning window in Y
;  atrous_tr, img, 20, decomp=imgd   ; Use atrous filter on the image
;  ac_atrous = convol_fft(imgd[*,*,6], temp, /auto_correlation)  ; take the autocorrelation of the filtered image
;  ;window = make_array(Nx, Ny, value=0.0)  ; Allocate space for window
;  window = gaussian_function([4,4], width=100, maximum=1)  ; Construct hanning window
;  psf_guess = window * ac_atrous
  
  atv, psf_guess
  
  ; Allocate intital guesses
  f_k = make_array(Nx, Ny, value=1.0)   ; intial guess for the object
  g_k = make_array(Nx, Ny, value=0.0) + gaussian_function([6,6], width=2048, maximum=1)   ;intial guess for the PSF
  
  atv, g_k
  
  for n = 0, iterations - 1 do begin
    
    ; Compute the next iteration of the PSF
    for i = 0, m - 1 do begin
      
      g_k = convol_fft((img / convol_fft(g_k, f_k)), reverse(reverse(f_k, 2))) * g_k
      
    endfor
    
    ; Compute the next iteration of the image
    for i = 0, m - 1 do begin
      
      f_k = convol_fft( (img / convol_fft(f_k, g_k)), reverse(reverse(g_k, 2)) ) * f_k
      
    endfor
    
    
  endfor

end