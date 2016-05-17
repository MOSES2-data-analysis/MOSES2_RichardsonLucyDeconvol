pro atrous_tr,  image, nsc, decomposition = decomp, filter = filter, $
             n_scales = n_scales, check = check
;+
; NAME:
;   ATROUS
; PURPOSE:
;   Perform a 2-D "a trous" wavelet decomposition on an image.
;
; CALLING SEQUENCE:
;   ATROUS, image [, decomposition = decomposition, $
;                 filter = filter, n_scales = n_scales, /check]
;
; INPUTS:
;   IMAGE -- A 2-D Image to Filter
;
; KEYWORD PARAMETERS:
;   FILTER -- A 1D (!!) array representing the filter to be used.
;             Defaults to a B_3 spline filter (see Stark & Murtaugh
;             "Astronomical Image and Data Analysis", Spinger-Verlag,
;             2002, Appendix A)
;   N_SCALES -- Set to name of variable to receive the number of
;               scales performed by the decomposition.
;   CHECK -- Check number of scales to be performed and return.
; OUTPUTS: 
;   DECOMPOSITION -- A 3-D array with scale running along the 3rd axis
;                    (large scales -> small scales). The first plane
;                    of the array is the smoothed image. To recover
;                    the image at any scale, just total the array
;                    along the 3rd dimension up to the scale desired.
;                  
;
; RESTRICTIONS:
;   Uses FFT convolutions which edge-wrap instead of mirroring the
;   edges as suggested by Stark & Mutaugh.  Wait for it.
;
; MODIFICATION HISTORY:
;
;       Mon Oct 6 11:49:50 2003, Erik Rosolowsky <eros@cosmic>
;		Written.
;       July 22, 2012 Thomas Rust
;            Implemented mirroring for FFT's, which removes edge effects
;-


; Start with simple filter
;  filter = [0.25, 0.5, 0.25]


  if n_elements(filter) eq 0 then filter = 1./[16, 4, 8/3., 4, 16]
  fmat = filter#transpose(filter)
  sz = size(image)

  n_scales = floor(alog((sz[1] < sz[2])/n_elements(filter))/alog(2))
  if keyword_set(check) then return
  decomp = fltarr(sz[1], sz[2], n_scales+1)

  if nsc gt n_scales then nsc = n_scales

  im = image
  for k = 0, nsc-1 do begin
; Smooth image with a convolution by a filter
    
    im2 = fltarr(2*sz(1),2*sz(2))
    im2(0:sz(1)-1,0:sz(2)-1)=im
    im2(0:sz(1)-1,sz(2):*)=reverse(im,2)
    im2(sz(1):*,0:sz(2)-1)=reverse(im,1)
    im2(sz(1):*,sz(2):*)=reverse(reverse(im,1),2)

    fmat2 = fltarr(2*sz(1),2*sz(2))
    sf = size(fmat)
    fmat2(sz(1)-(sf(1)-1)/2:sz(1)+(sf(1)-1)/2,$
          sz(2)-(sf(2)-1)/2:sz(2)+(sf(2)-1)/2) = fmat
    fmat2 = shift(fmat2,-sz(1),-sz(2))
    
    imq = fft(im2,-1)
    fmat2q = fft(fmat2,-1)
    smoothed = real_part(fft(imq*fmat2q,/inverse))*(2*float(sz(1))*2*float(sz(2)))
    smoothed = smoothed(0:sz(1)-1,0:sz(2)-1)
   
    decomp(*,*,n_scales-k) = im-smoothed
    
    im = smoothed;(smoothed)(0:sz(1)-1,0:sz(2)-1)

; Generate new filter
    newfilter = fltarr(2*n_elements(filter)-1) 
    newfilter[2*findgen(n_elements(filter))] = filter
; note filter is padded with zeros between the images
    fmat = newfilter#transpose(newfilter)
    filter = newfilter
  endfor

; Stick last smootheded image into end of array
  decomp[*, *, 0] = smoothed(0:sz(1)-1,0:sz(2)-1)

  return
end
