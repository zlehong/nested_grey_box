
# ClimaPlus validation study, copy made at March 17 2021
# Alpha Yacob Arsano
########################

import os, json, zipfile
from ClimaBox import CPlus_Climate as cc
from ClimaBox import climabox_ng as cb
import math
# import data_climates as dc

import pvlib
import pandas as pd
import numpy as np

# from Climate_aa_sam import RadFacade

from matplotlib import pyplot as plt


################## PART 1: EWP TO JSON
#######copied from data_climates file
#######Location of json folder, all EPW files are stored as json files. Place all EPW files for multiple simulations.

#dstDirJ = "B:/ClimateOneBuildingOrg/ashrae_drycold" ###Folder where json files are stored
#dstDir = 'B:/ClimateOneBuildingOrg/ashrae_drycold' ###Folder where ewp files are stored

dstDirJ = "C:/Users/zoele/Dropbox (MIT)/ML Tests/WeatherFiles/JSON"
dstDir = 'C:/Users/zoele/Dropbox (MIT)/ML Tests/WeatherFiles'

Qh_total, Qc_total = [],[]
hour_cc, tempout_cc = [], []
peopleLoad_cc,heatLoad_cc,coolLoad_cc,tempIn_cc,eqpLoad_cc,lgtLoad_cc,radLoad_cc,infillossLoad_cc,infilgainLoad_cc=[],[],[],[],[],[],[],[],[]
AOI_S, alt_cc, azi_cc, hourang_cc,dirradHS_cc,dirradHN_cc, difradHS_cc=[],[],[],[],[],[],[] #angle of incidence on the south window
solarDirGL_cc, solarDifGL_cc = [],[]
SHGC_L_dif_full = [[17.0, 0.93], [41.0, 2.03], [55.0, 1.87], [67.0, 1.6], [78.0, 1.36], [90, 1.1]]
SHGC_L_dir_full =[[34.0, 0.94], [42.0, 0.93], [51.0, 0.89], [61.0, 0.82], [70.0, 0.67], [90, 0.49]]

def epwtojson(dstDir):
	epws = os.listdir(dstDir)
	for epw in epws:
		if epw.endswith(".epw"):
			epwf = dstDir + "/"+epw
			w = cc.wthImp(epwf)
			wi = w.dataWth()
			ww = cc.wrangleWth(wi)
			wthdata = json.dumps(ww.WTH())
			jsnfile_ = epw.replace(".epw",".json")
			jsnfile = dstDirJ + "/" + jsnfile_
			with open(jsnfile,'w') as f:
			    f.write(wthdata)


##############
################## PART 2: RUN CLIMABOX
def initialize_cb():
	dstDirJ = "C:/Users/zoele/Dropbox (MIT)/ML Tests/WeatherFiles/JSON"
	dstDir = 'C:/Users/zoele/Dropbox (MIT)/ML Tests/WeatherFiles'

	Qh_total, Qc_total = [],[]
	hour_cc, tempout_cc = [], []
	peopleLoad_cc,heatLoad_cc,coolLoad_cc,tempIn_cc,eqpLoad_cc,lgtLoad_cc,radLoad_cc,infillossLoad_cc,infilgainLoad_cc=[],[],[],[],[],[],[],[],[]
	AOI_S, alt_cc, azi_cc, hourang_cc,dirradHS_cc,dirradHN_cc, difradHS_cc=[],[],[],[],[],[],[] #angle of incidence on the south window
	solarDirGL_cc, solarDifGL_cc = [],[]
	SHGC_L_dif_full = [[17.0, 0.93], [41.0, 2.03], [55.0, 1.87], [67.0, 1.6], [78.0, 1.36], [90, 1.1]]
	SHGC_L_dir_full =[[34.0, 0.94], [42.0, 0.93], [51.0, 0.89], [61.0, 0.82], [70.0, 0.67], [90, 0.49]]

def runMultipleClimabox():
	jsons = os.listdir(dstDirJ)
	for js in jsons:
		jsonfile = dstDirJ+"/"+js
		print (jsonfile)
		with open(jsonfile) as wthj:
			wth = json.load(wthj)

		dd = wth['db']
		tempout_cc.append(wth['db'])
		hour_cc.append(wth['hour'])
		wind_d = wth['wd']
		wind_v = wth['wv']
		LPD, EPD = 4,6 
		ground_temp = 10 



		# L,W,H,infilt_ach = 8,6,2.7,0.001
		L,W,H,infilt_ach = 10,10,3,0.001


		#Coefficents of SHGC based Kmeans clustering study using Denver's climate
		SHGC_L_dif_full = [[17.0, 0.93], [41.0, 2.03], [55.0, 1.87], [67.0, 1.6], [78.0, 1.36], [90, 1.1]] #[[18.0, 0.83], [43.0, 1.65], [57.0, 1.58], [68.0, 1.4], [79.0, 1.23], [90, 0.98]]
		SHGC_L_dir_full = [[34.0, 0.94], [42.0, 0.93], [51.0, 0.89], [61.0, 0.82], [70.0, 0.67], [90, 0.49]]

		DoubleClr,DoubleLoEe2 = [2.67,0.703],[1.493,0.373] ##first:U-val, second:SHGC Source CS, Eplu,RSMEANS,Citherlet et al
		DoubleClr_CBbase = [3,0.8] #DoublePaneClrCB in idf file

		################# OLD HERE #################
		# 'sp_lower':18, 'sp_upper':27,"nv_setpoint":20,
		################# OLD HERE #################

		zone_input = {'L':L,'W':W,'H':H,'WWRs':0.4, 'WWRw':0, 'WWRn':0, 'WWRe':0,
		                'U_wall':1.2,'U_roof':0.3,'U_floor':"adiabatic",'U_glazing':DoubleClr_CBbase[0],
		                'V_dot':(infilt_ach*L*W*H/3600),
		                'A_TM_factor':1, 'T_TM': 0.12, 'cp_TM':1130, 'density_TM':1000, 'density_TM0':1400,
		                'h_convection':5, 'shgc_ldif':SHGC_L_dif_full, 'shgc_ldir':SHGC_L_dir_full, 'shgc':DoubleClr_CBbase[1],'shadingFac':0.8, #'shgc':0.703, 'shadingFac':0.75,
		                'sp_lower':20, 'sp_upper':24,"nv_setpoint":18,
		                'Hcop':1, 'Ccop':3, 'NVach':10, #NVach has to be >1 for NV to be used
		                'coef_CO2elec':0.275, 'coef_CO2gas':0.18, 'coef_Costgas':0.04, 'coef_Costelec':0.20}

		lat, lon, tz = round((float(wth['latitude'])),0), round((float(wth["longitude"]))*-1,0), round((float(wth["timezone"]))*-15,0)
		rf = cc.radFacade(lat,lon,tz,wth['radn'], wth['radd'],wth['radg'])
		# azi_cc.append([math.degrees(x) for x in rf.azimuth])
		azi_cc.append( rf.azimuth_full)
		# alt_cc.append([math.radians(x) for x in rf.altitude])
		
		alt_cc.append(rf.altitude_full)
		hourang_cc.append(rf.hourangle)
		# print (rf.hourangle[0:48])
		# TODO: ZOE EDIT HERE: 
		dirradHS_cc.append(rf.dirradS)
		dirradHN_cc.append(rf.dirradN)
		difradHS_cc.append(rf.difradS)
		################# OLD HERE #################
		# dirradHS_cc.append(rf.dirradHS)
		# dirradHN_cc.append(rf.dirradHN)
		# difradHS_cc.append(rf.difradHS)
		################# OLD HERE #################
		# print(wth['radd'][0:8760])
		print ()
		print ()
		# print(rf.difradHS[0:8760])
		# print(rf.declination[0:48])
		# print("hour angle    ",rf.hourangle[0:48])
		# print (wth['radg'][600:900])

		alt_surf, azi_surf = math.radians(90), math.radians(0)
		for i in range(0,len(rf.azimuth)):
			# _am = math.radians(90-alt_surf)
			# _as = (rf.altitude[i])
			# _AM = azi_surf
			# _AS= (rf.azimuth[i])
			# _aoi = math.cos(_am)*math.cos(_as)*math.cos(_AM-_AS) + math.sin(_am)*math.sin(_as)
			_aoi_a = -1*math.cos(math.radians(lat))*math.sin(rf.declination[i])*math.cos(azi_surf)
			_aoi_b = math.sin(math.radians(lat))*math.cos(rf.declination[i])*math.cos(math.radians(rf.hourangle[i]))*math.cos(azi_surf) 
			_aoi_c = math.cos(rf.declination[i])*math.sin(math.radians(rf.hourangle[i]))*math.sin(azi_surf)   #for vertical surface, source: solar enery engineering
			_AOI = round(math.degrees(math.acos(_aoi_a+_aoi_b+_aoi_c)),2)
			################# OLD HERE #################
			# if -90 < _AOI < 90 and rf.dirradHS[i]>0.5: AOI_S.append(_AOI)
			################# OLD HERE #################
			if -90 < _AOI < 90 and rf.dirradS[i]>0.5: AOI_S.append(_AOI)
			else:AOI_S.append(0)
		
		
		
		######
		#### def intGain(zone_input,dirqradHS_hr,dirqradHW_hr,dirqradHN_hr,dirqradHE_hr,difqradHS_hr,ground_temp,zone_temp,AOI_hr,varSHGCbool)
		### if varSHGCbool==0, static SHGC of selected glazing is used. Else if varSHGCbool is 1, then variable SHGC based on solar angle is used.

		rads = cc.RadFacadePvlib(jsonfile)

		run_ = cb.runRCClass(dd,wind_d,wind_v,rads,zone_input,LPD,EPD,ground_temp,AOI_S,1)
		#########
		
		EUI_H_ = run_.EUI_H
		Qh_total.append(EUI_H_*L*W/1000)
		EUI_C_ = run_.EUI_H
		Qc_total.append(EUI_C_*L*W/1000)

		infilgainLoad_cc.append(run_.QinfilgainL)
		infillossLoad_cc.append(run_.QinfillossL)

		radLoad_cc.append(run_.solarGL) #[-6])
		solarDirGL_cc.append(run_.solarDirGL)
		# print (solarDirGL_cc)
		solarDifGL_cc.append(run_.solarDifGL)

		eqpLoad_cc.append(run_.q_eqptL) #[-5])
		lgtLoad_cc.append(run_.q_lightL) #[-4])
		peopleLoad_cc.append(run_.q_lightL)
		heatLoad_cc.append(run_.qsystemL_h_hr)#[-1])
		coolLoad_cc.append(run_.qsystemL_c_hr)#[-2])
		tempIn_cc.append(run_.t_inL_q) #[3])

		q_pp_cp = round(sum(peopleLoad_cc[0]))
		q_light_cp = round(sum(lgtLoad_cc[0]))
		q_eqpt_cp = round(sum(eqpLoad_cc[0]))
		q_heat_cp = round(sum(heatLoad_cc[0]))
		q_cool_cp = round(sum(coolLoad_cc[0]))
		q_rad_cp = round(sum(radLoad_cc[0]))
		q_raddir_cp = round(sum(solarDirGL_cc[0]))
		q_raddif_cp = round(sum(solarDifGL_cc[0]))
		
		

		###using pvlib (adopted from sam: def(prep_model(wth)))
		# handle a json-format weather file
		print(zone_input)
	# return [Qh_total, Qc_total, tempIn_cc, infilgainLoad_cc, radLoad_cc,  coolLoad_cc, eqpLoad_cc, lgtLoad_cc, peopleLoad_cc, heatLoad_cc, coolLoad_cc]
	return run_

def runClimabox(jsonfile):
	print (jsonfile)
	with open(jsonfile) as wthj:
		wth = json.load(wthj)

	dd = wth['db']
	tempout_cc.append(wth['db'])
	hour_cc.append(wth['hour'])
	wind_d = wth['wd']
	wind_v = wth['wv']
	LPD, EPD = 4,6 
	ground_temp = 10 



	# L,W,H,infilt_ach = 8,6,2.7,0.001
	L,W,H,infilt_ach = 10.5,3,3,0.001


	#Coefficents of SHGC based Kmeans clustering study using Denver's climate
	SHGC_L_dif_full = [[17.0, 0.93], [41.0, 2.03], [55.0, 1.87], [67.0, 1.6], [78.0, 1.36], [90, 1.1]] #[[18.0, 0.83], [43.0, 1.65], [57.0, 1.58], [68.0, 1.4], [79.0, 1.23], [90, 0.98]]
	SHGC_L_dir_full = [[34.0, 0.94], [42.0, 0.93], [51.0, 0.89], [61.0, 0.82], [70.0, 0.67], [90, 0.49]]

	DoubleClr,DoubleLoEe2 = [2.67,0.703],[1.493,0.373] ##first:U-val, second:SHGC Source CS, Eplu,RSMEANS,Citherlet et al
	DoubleClr_CBbase = [3,0.8] #DoublePaneClrCB in idf file

	################# OLD HERE #################
	# 'sp_lower':18, 'sp_upper':27,"nv_setpoint":20,
	################# OLD HERE #################

	zone_input = {'L':L,'W':W,'H':H,'WWRs':0.4, 'WWRw':0, 'WWRn':0, 'WWRe':0,
					'U_wall':1.2,'U_roof':0.3,'U_floor':"adiabatic",'U_glazing':DoubleClr_CBbase[0],
					'V_dot':(infilt_ach*L*W*H/3600),
					'A_TM_factor':1, 'T_TM': 0.12, 'cp_TM':1130, 'density_TM':1000, 'density_TM0':1400,
					'h_convection':5, 'shgc_ldif':SHGC_L_dif_full, 'shgc_ldir':SHGC_L_dir_full, 'shgc':DoubleClr_CBbase[1],'shadingFac':0.8, #'shgc':0.703, 'shadingFac':0.75,
					'sp_lower':20, 'sp_upper':24,"nv_setpoint":18,
					'Hcop':1, 'Ccop':3, 'NVach':10, #NVach has to be >1 for NV to be used
					'coef_CO2elec':0.275, 'coef_CO2gas':0.18, 'coef_Costgas':0.04, 'coef_Costelec':0.20}

	lat, lon, tz = round((float(wth['latitude'])),0), round((float(wth["longitude"]))*-1,0), round((float(wth["timezone"]))*-15,0)
	rf = cc.radFacade(lat,lon,tz,wth['radn'], wth['radd'],wth['radg'])
	# azi_cc.append([math.degrees(x) for x in rf.azimuth])
	azi_cc.append( rf.azimuth_full)
	# alt_cc.append([math.radians(x) for x in rf.altitude])
	
	alt_cc.append(rf.altitude_full)
	hourang_cc.append(rf.hourangle)
	# print (rf.hourangle[0:48])
	# TODO: ZOE EDIT HERE: 
	dirradHS_cc.append(rf.dirradS)
	dirradHN_cc.append(rf.dirradN)
	difradHS_cc.append(rf.difradS)
	################# OLD HERE #################
	# dirradHS_cc.append(rf.dirradHS)
	# dirradHN_cc.append(rf.dirradHN)
	# difradHS_cc.append(rf.difradHS)
	################# OLD HERE #################
	# print(wth['radd'][0:8760])
	print ()
	print ()
	# print(rf.difradHS[0:8760])
	# print(rf.declination[0:48])
	# print("hour angle    ",rf.hourangle[0:48])
	# print (wth['radg'][600:900])

	alt_surf, azi_surf = math.radians(90), math.radians(0)
	for i in range(0,len(rf.azimuth)):
		# _am = math.radians(90-alt_surf)
		# _as = (rf.altitude[i])
		# _AM = azi_surf
		# _AS= (rf.azimuth[i])
		# _aoi = math.cos(_am)*math.cos(_as)*math.cos(_AM-_AS) + math.sin(_am)*math.sin(_as)
		_aoi_a = -1*math.cos(math.radians(lat))*math.sin(rf.declination[i])*math.cos(azi_surf)
		_aoi_b = math.sin(math.radians(lat))*math.cos(rf.declination[i])*math.cos(math.radians(rf.hourangle[i]))*math.cos(azi_surf) 
		_aoi_c = math.cos(rf.declination[i])*math.sin(math.radians(rf.hourangle[i]))*math.sin(azi_surf)   #for vertical surface, source: solar enery engineering
		_AOI = round(math.degrees(math.acos(_aoi_a+_aoi_b+_aoi_c)),2)
		################# OLD HERE #################
		# if -90 < _AOI < 90 and rf.dirradHS[i]>0.5: AOI_S.append(_AOI)
		################# OLD HERE #################
		if -90 < _AOI < 90 and rf.dirradS[i]>0.5: AOI_S.append(_AOI)
		else:AOI_S.append(0)
	
	
	
	######
	#### def intGain(zone_input,dirqradHS_hr,dirqradHW_hr,dirqradHN_hr,dirqradHE_hr,difqradHS_hr,ground_temp,zone_temp,AOI_hr,varSHGCbool)
	### if varSHGCbool==0, static SHGC of selected glazing is used. Else if varSHGCbool is 1, then variable SHGC based on solar angle is used.

	rads = cc.RadFacadePvlib(jsonfile)

	run_ = cb.runRCClass(dd,wind_d,wind_v,rads,zone_input,LPD,EPD,ground_temp,AOI_S,1)
	#########
	
	EUI_H_ = run_.EUI_H
	Qh_total.append(EUI_H_*L*W/1000)
	EUI_C_ = run_.EUI_H
	Qc_total.append(EUI_C_*L*W/1000)

	infilgainLoad_cc.append(run_.QinfilgainL)
	infillossLoad_cc.append(run_.QinfillossL)

	radLoad_cc.append(run_.solarGL) #[-6])
	solarDirGL_cc.append(run_.solarDirGL)
	# print (solarDirGL_cc)
	solarDifGL_cc.append(run_.solarDifGL)

	eqpLoad_cc.append(run_.q_eqptL) #[-5])
	lgtLoad_cc.append(run_.q_lightL) #[-4])
	peopleLoad_cc.append(run_.q_lightL)
	heatLoad_cc.append(run_.qsystemL_h_hr)#[-1])
	coolLoad_cc.append(run_.qsystemL_c_hr)#[-2])
	tempIn_cc.append(run_.t_inL_q) #[3])

	q_pp_cp = round(sum(peopleLoad_cc[0]))
	q_light_cp = round(sum(lgtLoad_cc[0]))
	q_eqpt_cp = round(sum(eqpLoad_cc[0]))
	q_heat_cp = round(sum(heatLoad_cc[0]))
	q_cool_cp = round(sum(coolLoad_cc[0]))
	q_rad_cp = round(sum(radLoad_cc[0]))
	q_raddir_cp = round(sum(solarDirGL_cc[0]))
	q_raddif_cp = round(sum(solarDifGL_cc[0]))
	
	###using pvlib (adopted from sam: def(prep_model(wth)))
	# handle a json-format weather file
	print(zone_input)
	# return [Qh_total, Qc_total, tempIn_cc, infilgainLoad_cc, radLoad_cc,  coolLoad_cc, eqpLoad_cc, lgtLoad_cc, peopleLoad_cc, heatLoad_cc, coolLoad_cc]
	return run_

if __name__ == "__main__":
	initialize_cb()
	runMultipleClimabox()
	jsonfile = 'C:/Users/zoele/Dropbox (MIT)/ML Tests/WeatherFiles/JSON/USA_FL_Miami.Intl.AP.722020_TMY3.json'
	rads = cc.RadFacadePvlib(jsonfile)
	print(len(rads.dirradN))#dirradS

	# wth_pd = pd.read_json(jsonfile).rename(
	#     columns={"radg": "ghi", "radn": "dni", "radd": "dhi"}
	# )
	# city = wth_pd["city"][0]
	# lat = wth_pd["latitude"][0]
	# lon = wth_pd["longitude"][0]
	# tz = f"Etc/GMT{int(-wth_pd['timezone'][0]):+d}"
	# wth_pd.index = pd.date_range("2018-01-01", freq="1H", periods=len(wth_pd))
	# wth_pd = wth_pd.tz_localize(tz)
	# rf_pvlib = cc.RadFacade(
	#     latitude=lat,
	#     longitude=lon,
	#     name=city,
	#     altitude=wth_pd.get("alt", 10),
	#     timezone=tz,
	#     weather=wth_pd[["ghi", "dni", "dhi"]]
	# )

	# dirradN_pvlib = rf_pvlib.get_irradiance_for_surface(90, 0).loc[:, "poa_direct"].values
	# dirradE_pvlib = (
	#     rf_pvlib.get_irradiance_for_surface(90, 90).loc[:, "poa_direct"].values
	# )
	# dirradS_pvlib = (
	#     rf_pvlib.get_irradiance_for_surface(90, 180).loc[:, "poa_direct"].values
	# )
	# dirradW_pvlib, difradS_pvlib = (
	#     rf_pvlib.get_irradiance_for_surface(90, 270)
	#     .loc[:, ["poa_direct", "poa_diffuse"]]
	#     .values.T
	# )
	# print ((dirradW_pvlib))
	# print (len(dirradW_pvlib[5040:5088]))

	# plt.plot(difradS_pvlib[5040:5168])
	# plt.plot(difradHS_cc[0][5040:5168])
	# plt.plot(dirradS_pvlib[5041:5169])
	# plt.plot(dirradHS_cc[0][5040:5168])
	# plt.plot(dirradN_pvlib[5041:5169])
	# plt.plot(dirradHN_cc[0][5040:5168])


	plt.show()
	# print (hour_cc[0][5040:5088])
	# print (tempout_cc[0][5040:5088]) #Two days in first week of July



	# print (tempIn_cc[0][5040:5088])
	# print (radLoad_cc[0][5040:5088])
	# print (eqpLoad_cc[0][5040:5088])
	# print (peopleLoad_cc[0][5040:5088])
	# print (lgtLoad_cc[0][5040:5088])
	# print (coolLoad_cc[0][5040:5088])
	# print (infilgainLoad_cc[0][5040:5088])
	# print (infillossLoad_cc[0][5040:5088])


	# print (tempIn_cc)
	# print (tempout_cc)
	# print (coolLoad_cc)
	# print (heatLoad_cc)
	# print ((peopleLoad_cc))
	# print (Qh_total)
	# print (Qc_total)
	
	
