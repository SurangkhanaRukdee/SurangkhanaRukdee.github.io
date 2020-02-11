var exptime = 600;
var expnumber = 1;
var minimum = 10;
var maximum = 18;
var signoise = 100;

var pwv = 1;
var seeing = 1;

var waveminY1 = 0.85;
var wavemaxY1 = 1.12;

var waveminY2 = 1.12;
var wavemaxY2 = 1.4;

var jmag = 12;
var temperature = 3000;

var restwave = 1.0;
var doppler = 0;
var lineflux = 1;
var linewidth = 100;

//<!--Constants-->
var C = 299792500.0           // [m s-1]
var h = 6.6260755E-34         // [J s]

// Parameters: Basic
var D_telescope   = 6.5     // 6.5for the full aperture. 
var T_ambient     = 275     // [K]
var PI = Math.PI
var W_slit = 1.0
var n_slit = 3.66

var R = [60000, 60000]   // [Y band]

var Tau_atmosphere_band = [0.85, 1.12]

// Read noise
var n_read = 5               // [electron, Read Noise, IR Array parameters]
// dark current
var n_dark = 0.02            // [electron sec-1]



var wave_band = [1.0, 1.25];
var delta_wave = [0.27, 0.32];
var band_min = [0.85, 1.12];
var band_max = [1.12, 1.4];

// [ph s-1 m-2 arcsec-2], Moon Light
var S_ZM = [2.03e-23,   1.53e-23]
var phi_zod = [55, 41]
// [ph s-1 m-2 arcsec-2], Moon Light
var phi_moon_60_half = [54, 30]
var phi_moon_60_full = [595, 277]
// [ph s-1 m-2 arcsec-2], Use H-band data for J and K-bands.
var phi_OH_band = [4486,   4486] 

// Parameters: in Cal_emissivity_total.py
var Tau_M_1                  = 0.95
var Tau_M_1_e                = [0.55, 0.55]                 // Tau_M_1 * np.square(D_telescope_e / D_telescope)
var Tau_M_2                  = [1.00, 1.00]
var Tau_M_3                  = [1.00, 1.00]
var Tau_Window               = [0.95, 0.95]                 // (and AO)
//Tau_slit_loss            = [0.64, 0.64]	        // only for point source
//Pupil_Stop = [0.90, 0.90]
var Doublet    = [0.96, 0.96]
var Collimator = [0.58, 0.58]
var Echelle    = [0.88, 0.88]
var X_disp     = [0.80, 0.75]
var camera     = [0.89, 0.89]
var detector   = [0.80, 0.80]

var Tau_Optics = [0,0]

for (i = 0; i < 2; i++) {
    Tau_Optics[i] = Tau_M_1_e[i]*Tau_Window[i] * Doublet[i]*Collimator[i]*Echelle[i]*X_disp[i]*camera[i]*detector[i];
	//alert(Tau_Optics[i])
}

// This function calculated background emission.
// This function is called by "n_thermal()", "Mode_SN_band()", "Plot()", "Plot_wavelegnth()", "Signal_to_Noise()", "Magnitude()" functions.
function Cal_Background(wave_um, Temperature){
	y = 1.4745E-50 * Math.pow((C / (wave_um * 0.000001)),3)/ (Math.exp(0.0000000000479922 * C / (wave_um * 0.000001)/ Temperature)-1)
   return y
}

function Cal_Template_Mag_band(J_template, T_template){
	mag = [0,0];
	Omega_template = S_ZM[1]/Cal_Background(wave_band[1], T_template)
   for (BAND = 0; BAND < 2; BAND++){
        y_flux = Omega_template * Cal_Background(wave_band[BAND], T_template)
        y_mag  = J_template - 2.5*Math.log10(y_flux/S_ZM[BAND])
        mag[BAND] = Math.round(y_mag*1e3)/ 1e3
	}	    
    return mag
}

function Y_template_calc() {
    var J_template = parseInt(document.getElementById('J_mag').value); 
    var T_template = parseInt(document.getElementById('T_bb').value);
    var Y_template = Cal_Template_Mag_band(J_template, T_template);
    
    document.getElementById("modal-content").innerHTML="Y[mag] = "+Y_template[0]
}

function Doppler(restwave, vshift){
    d = ((restwave * vshift) / (C * Math.pow(10, -3))) + restwave
    return d
}

function GaussFunction(wave_um, linewidth, restwave){
    sigma = linewidth / 2.35482
    A = (wave_um - restwave) / sigma
    g = (1/(sigma*Math.sqrt(2.0*PI)))* Math.exp(-0.5 * A * A)
    return g
}

function Cal_LSignal(lineflux, wave_um, linewidth, restwave, vshift){
	 line = lineflux * GaussFunction(wave_um, linewidth, Doppler(restwave, vshift))
    return line
}

// This is Signal function.
function Cal_Signal(mag, Tau_total, t_exp, n_exp, i){
    x = (t_exp * n_exp) *(PI * Math.pow((D_telescope/2),2))* Tau_total * S_ZM[i] * Math.pow(10,(-0.4 * mag))/ (h * R[i])   
    return x
}

// This is Noise function.
function Cal_Noise(t_n_zod, t_thermal, t_signal, t_exp, n_exp){
    x = Math.sqrt(2 * Math.pow(n_slit,2) * n_exp * (t_exp * (t_n_zod + t_thermal + n_dark)+ Math.pow(n_read,2))+ t_signal)
    return x
}

// This function calculate the photo-electrons from background emissions
// This function is called by "Signal_to_Noise()", "Magnitude()" functions.
function Cal_n_zod_OH(wave_band, Tau_total, phi_OH_band, phi_moon, seeing_type, BAND){
    // A Omega calculation
    // [arcsec], We are planning on around 85 mas for all wavelengths.
    // [pixel], 4 pixels in the spectral and spatial resolution elements
    // We assume that the image in the cross-dispersion direction has a FWHM equal to the slit width.
    A_Omega = PI * Math.pow(D_telescope/2,2) * Math.pow((W_slit/n_slit)/206265,2)  // [m^2 sr]
    Tau_Slitloss = Cal_Tau_slitloss(seeing_type)
    Tau_total = Tau_total / Tau_Slitloss // Tau for extended source
    if (phi_OH_band > 4000){
        y = (wave_band / delta_wave[BAND] / R[BAND]) * Tau_total * (A_Omega * Math.pow(206265,2)) * (phi_zod[BAND] + phi_moon + 0.1 * phi_OH_band)
    }
    else{
        y_zod = (wave_band / delta_wave[BAND] / R[BAND]) * Tau_total * (A_Omega * Math.pow(206265,2)) * (phi_zod[BAND] + phi_moon)
        y_OH  = Tau_total * (A_Omega * Math.pow(206265,2)) * (phi_OH_band)
        y = y_zod + y_OH
	 }    
    return y
}

// This function calculated background emission.
// This function is called by "n_thermal()", "Mode_SN_band()", "Plot()", "Plot_wavelegnth()", "Signal_to_Noise()", "Magnitude()" functions.
function Cal_Background(wave_um, Temperature){
    // Background Emission Calculation
    // Units in bands, [um], from <gmtnirs_AO_strawman_070706b.doc>
    // Background for Zodiacal and OH emissions
    y = 1.4745E-50 * Math.pow((C / (wave_um * 0.000001)),3)/ (Math.exp(0.0000000000479922 * C / (wave_um * 0.000001)/ Temperature)-1)
    return y
}

// This function calculate the thermal-electrons.
// This function is called by "Signal_to_Noise()", "Mode_SN_band()", "Plot()", "Plot_wavelegnth()", "Magnitude()" functions.
function Cal_n_thermal(emissivity, background, BAND){
    A_Omega = PI * Math.pow(D_telescope/2,2) * Math.pow((W_slit/n_slit)/206265,2)  // [m^2 sr]
    y = (A_Omega * emissivity * background) / (h * R[BAND])
    return y
}

// This function calculate tau from optical parameters. 
// This function is called by "Signal_to_Noise()", "Magnitude()" functions.
// Tau_Optics[] values are calculated seperately (see Cal_Tau_Optics())
function Cal_Tau_total(Tau_atmosphere, seeing_type, BAND){
    //Total transmission for point source
    Tau_Slitloss = Cal_Tau_slitloss(seeing_type)
    y = Tau_atmosphere * Tau_Optics[BAND] * Tau_Slitloss
    return y
}

function erf(x) {
  // save the sign of x
  var sign = (x >= 0) ? 1 : -1;
  x = Math.abs(x);

  // constants
  var a1 =  0.254829592;
  var a2 = -0.284496736;
  var a3 =  1.421413741;
  var a4 = -1.453152027;
  var a5 =  1.061405429;
  var p  =  0.3275911;

  // A&S formula 7.1.26
  var t = 1.0/(1.0 + p*x);
  var y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
  return sign * y; // erf(-x) = -erf(x);
}



// This function calculated emissivity from optical parameters
// This function is called by "Cal_n_thermal()", "Signal_to_Noise()", "Magnitude()" functions. 
// The emissivities at different wavelengths are different enough that this should be folded into the real calculator
function Cal_Tau_slitloss(seeing_type){
    // Calculate Tau slit loss for point source.
    sigma = seeing_type / 2.35482
    Y = 0.5 // half of slit width
    para = Y/(Math.sqrt(2)*(sigma))
    y = erf(para)
//    invexp = lambda t: Math.exp(-1*t*t)
//    y = (2 / Math.sqrt(PI)) * integrate.quad(invexp, 0, para)[0]    
    return y
}

function Cal_emissivity_total(Tau_atmosphere, Tau_total, seeing_type, BAND){
    // Total emissivity        
    var Eta_M1 = 0.25
    var Eta_M2 = 0.00
    var Eta_M3 = 0.00 
    var Eta_Window     = 0.05
    var Eta_atmosphere = 0.05  
    var Tau_Slitloss = Cal_Tau_slitloss(seeing_type)
    y = (Eta_atmosphere + (Eta_M1 + (Eta_M2 + (Eta_M3 + Eta_Window / Tau_Window[BAND]) / Tau_M_3[BAND]) / Tau_M_2[BAND]) / Tau_M_1_e[BAND]) / Tau_atmosphere * (Tau_total / Tau_Slitloss)
    return y 
}

// This is Cal_Signal-to-Noise function.
// This function is called "Mode_SN_band()", "Plot()", Plot_wavelength()" functions.
function Cal_SN(Tau_atmosphere, wave_um, phi_OH_band, phi_moon, mag, t_exp , n_exp, linesignal, seeing_type, BAND){
    Tau_band = Cal_Tau_total(Tau_atmosphere, seeing_type, BAND)
    Emissivity_band = Cal_emissivity_total(Tau_atmosphere, Tau_band, seeing_type, BAND)
    n_zod_OH_band = Cal_n_zod_OH(wave_um, Tau_band, phi_OH_band, phi_moon, seeing_type, BAND)
    B_band = Cal_Background(wave_um, T_ambient)
    n_thermal_band = Cal_n_thermal(Emissivity_band, B_band, BAND)
    y_signal = Cal_Signal(mag, Tau_band, t_exp, n_exp, BAND)
    
    y_line = (linesignal *  (n_exp * t_exp * (PI * Math.pow((D_telescope/2),2)) * Tau_band * Math.pow(wave_um,2) *  Math.pow(10, -6))) / (h * C * R[1])
    y_signal_c = y_signal + y_line
    y_noise  = Cal_Noise(n_zod_OH_band, n_thermal_band, y_signal_c, t_exp, n_exp)
    y = y_signal_c / y_noise
    return y
}


function Mode_SN_band(t_exp, n_exp, J_template, T_template, lineflux, linewidth, restwave, vshift, seeing_type, moon_type){
	 if (moon_type == 'No Moon'){
        phi_moon = [0, 0]}
    else if (moon_type == 'Half Moon'){
        phi_moon = phi_moon_60_half}
    else{
        phi_moon = phi_moon_60_full}
//    print 'moon = ', phi_moon
    mag = Cal_Template_Mag_band(J_template, T_template)
    SN = [0,0]
    for (BAND = 0; BAND < 2; BAND++){
        signal = Cal_SN(Tau_atmosphere_band[BAND], wave_band[BAND], phi_OH_band[BAND], phi_moon[BAND], mag[BAND], t_exp, n_exp, Cal_LSignal(lineflux, wave_band[BAND], linewidth, restwave, vshift), seeing_type, BAND)
        SN[BAND] = Math.round(signal * 100) / 100
	 }
    return SN
}

function Cal_MSignal(Tau_atmosphere, wave_um, mag, t_exp , n_exp, seeing_type, BAND){
    Tau_band = Cal_Tau_total(Tau_atmosphere, seeing_type, BAND)
    //print 'tau', Tau_band
    y_signal = Cal_Signal(mag, Tau_band, t_exp, n_exp, BAND) * Math.pow(10, 6)
    y_m      = (h * C * R[BAND]) / ((n_exp * t_exp * (PI * Math.pow((D_telescope/2),2)) * Tau_band * Math.pow(wave_um,2)))
    y = y_signal * y_m 
    return y
}
// Test with calibrated sginal
function Mode_MS_band(t_exp, n_exp, J_template, T_template, seeing_type){
    mag = Cal_Template_Mag_band(J_template, T_template)
    SN = [0,0]
    for (BAND = 0; BAND < 2; BAND++){
        //y = Cal_SN(Tau_atmosphere_band[BAND], wave_band[BAND], phi_OH_band[BAND], mag[BAND], t_exp, n_exp, BAND)
        y = Cal_MSignal(Tau_atmosphere_band[BAND], wave_band[BAND], mag[BAND], t_exp , n_exp, seeing_type, BAND)
        SN[BAND] = y 
    }
    return SN
}

function arange(start, stop, step) {
    if (stop == null) {
      stop = start || 0;
      start = 0;
    }
    if (!step) {
      step = stop < start ? -1 : 1;
    }
    var length = Math.max(Math.ceil((stop - start) / step), 0);
    var range = Array(length);
    for (var idx = 0; idx < length; idx++, start += step) {
      range[idx] = start;
    }
    return range;
}

//  This function calculate Signal-to-Noise vs. Magnitude values.
function Mode_SN_Mag_Plot(t_exp, n_exp, mag_max, mag_min, lineflux, linewidth, restwave, vshift, seeing_type, moon_type){
    if (moon_type == 'No Moon'){
        phi_moon = [0, 0]}
    else if (moon_type == 'Half Moon'){
        phi_moon = phi_moon_60_half}
    else{
        phi_moon = phi_moon_60_full}

	 var series = []    
    var mag = arange(mag_min, mag_max + 0.5, 0.5)
    //var mag = [10,11,12,13,14,15,16,17,18]
    //S_N = [arange(v_mag_min, v_mag_max + 0.5, 0.5), arange(v_mag_min, v_mag_max + 0.5, 0.5)]
	 
    
    for (BAND = 0; BAND < 2; BAND++){
    	  var band_series = [];
        for (i = 0; i < mag.length; i++){
			    var v = Cal_SN(Tau_atmosphere_band[BAND], wave_band[BAND], phi_OH_band[BAND], phi_moon[BAND], mag[i], t_exp , n_exp, Cal_LSignal(lineflux, wave_band[BAND], linewidth, restwave, vshift), seeing_type, BAND) 
             band_series.push([mag[i], Math.log10(v)]);
        }   
         series.push(band_series);
    }
    return series;
}

// This function calculates atmosphere transmission vs. wavelength.
// This function is called by "Plot_wavelength" function.
function Get_Tau_atmo(atmo, pwv){
//    N_data = len(wavelength)
    y = []
    if (pwv >= 1 && pwv <= 2){
        for (i = 0; i < (atmo.length); i++){
	        y[i] = atmo[i][1] + (pwv - 1) * (atmo[i][2] - atmo[i][1])
        }
     }
     return y
}

// This function calculate Magnitude values vs wavelength.
// This function is called by "Mode_SN_band()", "Plot_wavelegnth()" function.
function Cal_Template_Mag_wave(x, J_template, T_template, i){
    Omega_template = S_ZM[1] / Cal_Background(wave_band[1], T_template)
    flux = Omega_template * Cal_Background(x, T_template)
    magnitude = J_template - 2.5*Math.log10(flux/S_ZM[i])
    return magnitude
}

function Mode_SN_wave_Plot(t_exp, n_exp, pwv_type, min_x, max_x, J_template, T_template, lineflux, linewidth, restwave, vshift, seeing_type, moon_type){
    if (moon_type == 'No Moon'){
        phi_moon = [0, 0]}
    else if (moon_type == 'Half Moon'){
        phi_moon = phi_moon_60_half}
    else{
        phi_moon = phi_moon_60_full}

	 series = [];
    if (min_x >= band_min[0] && max_x <= band_max[0]){
        Tau_atmosphere = Get_Tau_atmo(ATMO_TAO_Y, pwv_type);
       
        for (i = 0; i < (ATMO_TAO_Y.length); i++){
            BAND = 0
            Wavelength= ATMO_TAO_Y[i][0];
            if (Wavelength >= min_x && Wavelength <= max_x){
            S_N = Cal_SN(Tau_atmosphere[i], ATMO_TAO_Y[i][0], OH_Y[i][1], phi_moon[BAND], Cal_Template_Mag_wave(ATMO_TAO_Y[i][0], J_template, T_template, BAND), t_exp , n_exp, Cal_LSignal(lineflux, ATMO_TAO_Y[i][0], linewidth, restwave, vshift), seeing_type, BAND)
            series[i] = [Wavelength, S_N];
        }}
        
        //print S_N
    }
    return [series]
    //Plot
    //if (min_x >= band_min[0] && max_x <= band_max[0]){
        //Title_para = ' (t=' + '%d' %t_exp + 'sec, N=' + '%d'%n_exp + ', PWV=' + '%d' % pwv_type + ')'
        //Plot_SN_Wave(Wavelength, S_N, min_x, max_x) 
    //}	
}    


function Plot_SN_Mag(){
	 var v_exptime = parseFloat(document.getElementById('exptime').value); 
    var v_expnumber = parseFloat(document.getElementById('expnumber').value);
	 var v_seeing = parseFloat(document.getElementById('seeingnumber').value); 
    var moon_get = (document.getElementById('moon').value);
    var v_restwave = parseFloat(document.getElementById('rest_wavelength').value); 
    var v_lineflux = parseFloat(document.getElementById('line_flux').value)*1e-18;
    var linewidth = parseFloat(document.getElementById('line_width').value);
    var v_doppler = parseFloat(document.getElementById('doppler_shift').value);
    
 //Change unit to Km s-1
    var v_linewidth = (v_restwave * (linewidth)) / (299792500 * Math.pow(10, -3)); // Unit in um
    var o_doppler = (v_restwave * (v_doppler)) / (299792500 * Math.pow(10, -3)); // Unit in um

 // Derive for calibrate line signal peak and line equivalent width based on input of line width and line flux
    var v_peaksignal = (v_lineflux * 2.35482) / (Math.sqrt(2.0*PI)* v_linewidth);
    //var v_ewidth = (Math.sqrt(2.0*PI)* v_linewidth) / 2.35482;    
    
    var v_mag_min = parseFloat(document.getElementById("mag_min").value)
    var v_mag_max = parseFloat(document.getElementById("mag_max").value) 
	 series = Mode_SN_Mag_Plot(v_exptime, v_expnumber, v_mag_max, v_mag_min, v_lineflux, v_linewidth, v_restwave, v_doppler, v_seeing, moon_get)    

	 var options = {	 	
      /* yaxis: {
	    ticks: [[0, "10<sup>0</sup>"], [1, "10<sup>1</sup>"], [2, "10<sup>2</sup>"], [3, "10<sup>3</sup>"]]
		}, */
		xaxis: {
			   show: true,
            axisLabel: "Magnitude",
        },
      yaxis: {
      	       show: true,
                ticks: [[0, "10<sup>0</sup>"], [1, "10<sup>1</sup>"], [2, "10<sup>2</sup>"], [3, "10<sup>3</sup>"]],
                axisLabel: "Signal to Noise",
            },
		grid: {
			show: true
		},
	 }
    $.plot("#placeholder", [series[0]], options);
}


function Plot_SN_Wave(){
	 var v_exptime = parseFloat(document.getElementById('exptime').value); 
    var v_expnumber = parseFloat(document.getElementById('expnumber').value);
	 var v_seeing = parseFloat(document.getElementById('seeingnumber').value);
	 var pwv = parseFloat(document.getElementById('pwv').value);  
    var moon_get = (document.getElementById('moon').value);
    var v_jmag = parseFloat(document.getElementById('J_mag').value); 
    var v_temperature = parseFloat(document.getElementById('T_bb').value);
    var v_restwave = parseFloat(document.getElementById('rest_wavelength').value); 
    var v_lineflux = parseFloat(document.getElementById('line_flux').value)*1e-18;
    var linewidth = parseFloat(document.getElementById('line_width').value);
    var v_doppler = parseFloat(document.getElementById('doppler_shift').value);
    
 //Change unit to Km s-1
    var v_linewidth = (v_restwave * (linewidth)) / (299792500 * Math.pow(10, -3)); // Unit in um
    var o_doppler = (v_restwave * (v_doppler)) / (299792500 * Math.pow(10, -3)); // Unit in um

 // Derive for calibrate line signal peak and line equivalent width based on input of line width and line flux
    var v_peaksignal = (v_lineflux * 2.35482) / (Math.sqrt(2.0*PI)* v_linewidth);
    //var v_ewidth = (Math.sqrt(2.0*PI)* v_linewidth) / 2.35482;    
    
    var min_x = parseFloat(document.getElementById("Wave_min").value);
    var max_x = parseFloat(document.getElementById("Wave_max").value); 
	 series = Mode_SN_wave_Plot(v_exptime, v_expnumber, pwv, min_x, max_x, v_jmag, v_temperature, v_lineflux, linewidth, v_restwave, v_doppler, v_seeing, moon_get)    

	 var options = {	 	
      /* yaxis: {
	    ticks: [[0, "10<sup>0</sup>"], [1, "10<sup>1</sup>"], [2, "10<sup>2</sup>"], [3, "10<sup>3</sup>"]]
		}, */
		xaxis: {
			   show: true,
            axisLabel: "Wavelength [um]",
        },
      yaxis: {
      	     show: true,
              axisLabel: "S/N per Resolution Element",
            },
		grid: {
			show: true
		},
	 }
    $.plot("#placeholder", [series[0]], options);
}
    	         
function Calc_SN(){
	 var v_exptime = parseFloat(document.getElementById('exptime').value); 
    var v_expnumber = parseFloat(document.getElementById('expnumber').value);
	 var v_seeing = parseFloat(document.getElementById('seeingnumber').value); 
    var moon_get = (document.getElementById('moon').value);
	 var v_jmag = parseFloat(document.getElementById('J_mag').value); 
    var v_temperature = parseFloat(document.getElementById('T_bb').value);
    var v_restwave = parseFloat(document.getElementById('rest_wavelength').value); 
    var v_lineflux = parseFloat(document.getElementById('line_flux').value)*1e-18;
    var linewidth = parseFloat(document.getElementById('line_width').value);
    var v_doppler = parseFloat(document.getElementById('doppler_shift').value);
    
 //Change unit to Km s-1
    var v_linewidth = (v_restwave * (linewidth)) / (299792500 * Math.pow(10, -3)); // Unit in um
    var o_doppler = (v_restwave * (v_doppler)) / (299792500 * Math.pow(10, -3)); // Unit in um

 // Derive for calibrate line signal peak and line equivalent width based on input of line width and line flux
    var v_peaksignal = (v_lineflux * 2.35482) / (Math.sqrt(2.0*PI)* v_linewidth);
    //var v_ewidth = (Math.sqrt(2.0*PI)* v_linewidth) / 2.35482;    
    
    var v_snr = Mode_SN_band(v_exptime, v_expnumber, v_jmag, v_temperature, v_lineflux, v_linewidth, v_restwave, v_doppler, v_seeing, moon_get);
    var v_mag = Cal_Template_Mag_band(v_jmag, v_temperature);
    document.getElementById("Out_SN").value= v_snr[0];
    document.getElementById("Out_SN").style.backgroundColor = "lightblue";
}


