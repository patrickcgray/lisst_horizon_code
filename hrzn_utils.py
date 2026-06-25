import pandas as pd
import geopandas as gpd
import numpy as np
from datetime import date
import matplotlib.pyplot as plt
import glob
import scipy.io
from scipy import interpolate


import warnings
warnings.filterwarnings("ignore")

def load_in_hrzn_parameters():
    # Compute Angles for Forward and Side Scatter Detectors
    # this is based of Wayne Slade's raytracing work
    
    fn = 'anglesSsd.mat'
    mat = scipy.io.loadmat(fn)
    angles_ssd = mat['anglesSsd'][0]
    
    fn = 'angles_fwd.mat'
    mat = scipy.io.loadmat(fn)
    angles_fwd = mat['anglesFwd'][0]
    
    fn = 'anglesSsdWidth.mat'
    mat = scipy.io.loadmat(fn)
    widths_ssd = mat['anglesSsdWidth'][0]
    
    fn = 'ssdPathlengthMeters.mat'
    mat = scipy.io.loadmat(fn)
    pathlength_ssd = mat['pathlengthSsdMeters'][0]
    return(angles_ssd, angles_fwd, widths_ssd, pathlength_ssd)


def open_raw_files(csv_files):
    """
    inputs:
        csv_files: a list of strings denoting all LISST Horizon files to process 
    
    outputs:
        a pandas dataframe of these files 
    """
    df_list = []

    # Read and concatenate each CSV file
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        if len(df) == 0: # there are lots of empty files in there
            continue
        df_list.append(df)
    raw = pd.concat(df_list, ignore_index=True)

    # convert the sametime to a value that can be loaded directly with pandas to_datetime()
    raw.sampleTime = raw.sampleTime/86400
    
    return(raw)

def organize_sample_ids(raw):
    """
    inputs:
        raw dataframe
    
    outputs:
        sample_ids of all samples
    """
    # Determine the start and end indices of each sample
    newSampleType = np.concatenate([[0],np.diff(raw.type)!=0]) # TODO I think maybe the first sample should be 1??
     # TODO update this to python datetime
    newSampleByTime = np.concatenate([[0],np.diff(raw.sampleTime) > 3.4722e-04]) # note that datenum(0,0,0,0,0,30) = 3.4722e-04
    newSampleByPol = np.concatenate([[0],np.diff(raw.polarizationAngle)!=0])

    sampleStart = np.argwhere(newSampleType | newSampleByTime).flatten() # this use to also count pols are new samples but I'm not  | newSampleByPol
    sampleEnd = sampleStart-1
    sampleStart = np.concatenate([[0], sampleStart])
    
    sampleEnd = np.concatenate([sampleEnd, [len(raw.sampleTime)-1]])
    numSamples = len(sampleStart)

    # create the sample id numbers
    sample_ids = []
    sample_id = 0
    
    for k in range(numSamples):
        for subsample in range(sampleStart[k],sampleEnd[k]+1):
            sample_ids.append(sample_id)
        sample_id += 1     

    return(sample_ids)

def hrzn_load(raw):
    """
    inputs:
        raw: a pandas dataframe of raw files
    outputs:
        
    """
    print('we have: ' + str(len(raw)) + ' measurements to process')
    # Subtract Off Dark Measurements
    # reference and transmitted power
    LP = raw.RawTLight - raw.RawTDark # dark corrected
    LREF = raw.LRefLight - raw.LRefDark # dark corrected
    
    # forward scatter detectors
    scatFwd = raw.filter(regex='FSLight') - raw.filter(regex='FSDark').values # dark corrected
    
    # side scatter detectors
    tmpConv = raw.filter(regex='SSIntTime').values*20 # factor converts from msec back to CONV units; assumed 50usec/CONV
    tmpRange = raw.filter(regex='SSRange').values # range capacitor factor

    # print(tmpConv[raw.polarizationAngle==0])
    # print(tmpRange[raw.polarizationAngle==0])

    # print(tmpConv[raw.polarizationAngle==90])
    # print(tmpRange[raw.polarizationAngle==90])
    
    counterSsdLight = np.array([raw.filter(regex='SS1CountL'), raw.filter(regex='SS2CountL'), raw.filter(regex='SS3CountL')]).T[0] # need a shape of (samples, det boards) or (N, 3)
    counterSsdDark = np.array([raw.filter(regex='SS1CountD'), raw.filter(regex='SS2CountD'), raw.filter(regex='SS3CountD')]).T[0] # need a shape of (samples, det boards) or (N, 3)
    
    # unclear why these were put in a temp var before just being put in an array below but doing it consistently here
    # tmpSsdRawLight{1} = raw.filter(regex='SS1Light')
    # tmpSsdRawLight{2}  = raw.filter(regex='SS2Light')
    # tmpSsdRawLight{3}  = raw.filter(regex='SS3Light');
    
    # tmpSsdRawDark{1} = raw{:,startsWith(VarNames,'SS1Dark')};
    # tmpSsdRawDark{2} = raw{:,startsWith(VarNames,'SS2Dark')};
    # tmpSsdRawDark{3} = raw{:,startsWith(VarNames,'SS3Dark')};
    
    # this is not called so not implemented right now
    # % change ssb order, if specified
    # if ssbOrder~=[1 2 3]
    #     disp('** Reordering SideScat Boards **')
    #     disp(ssbOrder)
    #     tmpSsdRawLight = tmpSsdRawLight(ssbOrder);
    #     tmpSsdRawDark = tmpSsdRawDark(ssbOrder);
    #     tmpConv = tmpConv(:,ssbOrder);
    #     tmpRange = tmpRange(:,ssbOrder);
    #     counterSsdLight = counterSsdLight(:,ssbOrder);
    #     counterSsdDark = counterSsdDark(:,ssbOrder);
    # end
    
    scatSsdRawLight = np.array([raw.filter(regex='SS1Light'), raw.filter(regex='SS2Light'), raw.filter(regex='SS3Light')],dtype=np.float32)
    scatSsdRawDark = np.array([raw.filter(regex='SS1Dark'), raw.filter(regex='SS2Dark'), raw.filter(regex='SS3Dark')],dtype=np.float32)    
    
    # plt.plot(scatSsdRawLight.swapaxes(0,1)[raw.polarizationAngle==0].reshape(-1,24).T,c='k',label='SSDRawLight, pol=0')
    # plt.plot(scatSsdRawDark.swapaxes(0,1)[raw.polarizationAngle==0].reshape(-1,24).T,c='orange',label='SSDRawDark, pol=0')

    # plt.plot(scatSsdRawLight.swapaxes(0,1)[raw.polarizationAngle==90].reshape(-1,24).T,c='k',label='SSDRawLight, pol=90', ls='--')
    # plt.plot(scatSsdRawDark.swapaxes(0,1)[raw.polarizationAngle==90].reshape(-1,24).T,c='orange',label='SSDRawDark, pol=90', ls='--')
    # plt.yscale('log')
    # plt.ylim(1000,1300000)
    # plt.xlabel('Detector #')
    # plt.ylabel('Counts')
    # plt.show()
    
    
    # TODO should this be 0.01 or something? Or 10?
    badZeros = scatSsdRawLight-scatSsdRawDark<=10 # PG changed it here to only delete the individual detector not the whole sample or det array
    
    print('** Data rows removed due to light-dark error **')
    print('   ' + str(np.sum(badZeros)) + ' row(s) removed')
    
    scatSsdRawLight[badZeros] = np.nan # TODO don't remove whole row just make this a zero??
    scatSsdRawDark[badZeros] = np.nan
    
    scatSsdRawLight[np.any(scatSsdRawLight>(0.95*2**20),axis=2)] = np.nan # remove data near saturation
    scatSsdRawDark[np.any(scatSsdRawDark>(0.95*2**20),axis=2)] = np.nan # remove data near saturation
    
    print('** Data rows removed due to saturation **')
    print('   scatSsdRawLight, ' + str(np.sum(np.any(np.isnan(scatSsdRawLight),2))) + ' row(s) removed')
    print('   scatSsdRawDark, '  + str(np.sum(np.any(np.isnan(scatSsdRawDark),2)))  + ' row(s) removed')
    
    conv = np.repeat(tmpConv.T[:, :, np.newaxis],8,axis=2)
    rangeMat = np.repeat(tmpRange.T[:, :, np.newaxis],8,axis=2)
        
        # % maxWord = 2^20;
        # % fsRange = [50 100 150 200 250 300 350].*1e-12;  % in Coulombs
        # % convUnit = 55.1e-6; % converts to seconds
        # % sideScat1 = sideScatRaw1./maxWord .* repmat(fsRange(dat.range(:,1))',1,8);
        # % sideScat1 = sideScat1 ./ repmat(dat.conv(:,1)*convUnit,1,8);
    
    scatSsdLight = scatSsdRawLight*rangeMat/conv # TODO why aren't we dividing by 2^20 here?
    scatSsdDark = scatSsdRawDark*rangeMat/conv
    
    scatSsd = scatSsdLight - scatSsdDark

    # identify all the unique samples and create an id for each
    sample_ids = organize_sample_ids(raw)

    # now rename to the full names and add in the calculated values
    loaded = raw.copy()

    loaded = loaded.rename(columns={"type": "sampleType", 
        "Twater": "tempWaterCelcius",
        "Tcell": "tempCellCelcius",
        "Telec": "tempIntCelcius",
        "humidity": "humiIntPct",
        "VSup": "voltageSupply",
        "VHoldup": "voltageHoldup",
        "ManifoldPSI": "pressureManifoldPsi",
        "PumpPSI": "pressurePumpPsi",
        "PumpRPM": "speedPumpRpm",
        "MixerSpeed": "speedMixerRaw",
        "tank1": "tankLevelsPct1",
        "tank2": "tankLevelsPct2",
        "tank3": "tankLevelsPct3",
    })
    
    loaded['LP'] = LP
    loaded['LREF'] = LREF
    # this stores a len = 3 np.array in the tank column for each row
    # TODO why do we rename tanks and then set it to an array called tanks?
    loaded['tanks'] = list(raw.filter(regex='tank').values)
    loaded['scatFwd'] = list(scatFwd.values)
    loaded['scatSsd'] = list(scatSsd.swapaxes(0,1).reshape(-1,24))
    loaded['scatSsdDark'] = list(scatSsdDark.swapaxes(0,1).reshape(-1,24))
    loaded['counterSsdLight'] = list(counterSsdLight)
    loaded['counterSsdDark'] = list(counterSsdDark)
    loaded['conv'] = list(conv.swapaxes(0,1).reshape(-1,24))
    loaded['range'] = list(rangeMat.swapaxes(0,1).reshape(-1,24))

    loaded['sample_id'] = sample_ids

    # sample_dts = pd.to_datetime(loaded.sampleTime-719529, unit='D')
    sample_dts = pd.to_datetime(loaded.sampleTime, unit='D')
    loaded['datetime'] = sample_dts
    
    loaded['utc_dt'] = loaded['datetime']
    loaded  = loaded.set_index('utc_dt')

    return(loaded)

def hrzn_process(loaded, using_fscat=False, apply_atten_corr=True, subtract_zscat=True, empirical_corr_factor=True, specific_zscat_file=False, transmission_factor=1, nnf_factor=0.2,beamc_forced=None):
    # TODOS
    # this currently doesn't correct for the attenuation of the fwd scatter

    angles_ssd, angles_fwd, widths_ssd, pathlength_ssd = load_in_hrzn_parameters()
    
    dcal0pol = pd.read_csv('Hrzn_dcal_0DegPol.asc',header=None, delimiter='\t').values[0]
    dcal90pol = pd.read_csv('Hrzn_dcal_90DegPol.asc',header=None, delimiter=',').values[0]

    loaded['sample_id_groupby']= loaded.sample_id
    loaded_grouped = loaded[['sample_id_groupby', 'sample_id','sampleType', 'datetime']].groupby('sample_id_groupby').median()

    loaded_grouped['utc_dt'] = loaded_grouped['datetime']
    loaded_grouped  = loaded_grouped.set_index('utc_dt')

    loaded_grouped_samples = loaded_grouped[loaded_grouped.sampleType == 0]
    fscat_df = loaded_grouped[loaded_grouped.sampleType == 2]
    zscat_df = loaded_grouped[loaded_grouped.sampleType == 3]

    zscat_df = zscat_df.add_suffix('_zscat')
    fscat_df = fscat_df.add_suffix('_fscat')

    tol = pd.Timedelta('1440 minute') # currently requiring it be within a day to include the zscat
    # TODO maybe I want the nearest to the left??
    # not requiring a tolerance right now
    loaded_grouped_zscat = pd.merge_asof(left=loaded_grouped_samples,right=zscat_df,right_index=True,left_index=True,direction='nearest')#,tolerance=tol)
    loaded_grouped_fscat = pd.merge_asof(left=loaded_grouped_zscat,right=fscat_df,right_index=True,left_index=True,direction='nearest')#,tolerance=tol)
    
    sample_ids = loaded_grouped_fscat.sample_id.values
    matched_fscat_ids = loaded_grouped_fscat.sample_id_fscat.values
    matched_zscat_ids = loaded_grouped_fscat.sample_id_zscat.values
    
    fscat_lookup = dict(zip(sample_ids, matched_fscat_ids))
    zscat_lookup = dict(zip(sample_ids, matched_zscat_ids))

    ## Calculate Corrected Scatter (cscat)
    processed_samps = []

    monte_carlo_results_0 = []
    monte_carlo_results_90 = []
    cscatFwd_monte_carlo = []
    
    for sid in sample_ids:        
        samp_df = loaded[loaded.sample_id == sid]
        
        if len(samp_df.polarizationAngle.unique()) < 2: # this means it can't have both polarizations
            continue
    
        baseline_id_lookup = zscat_lookup
        if using_fscat:
            baseline_id_lookup = fscat_lookup
    
        zscat_full = loaded[loaded.sample_id == baseline_id_lookup[sid]]
        if len(specific_zscat_file)>0:
            zscat_full = specific_zscat_file
            print(zscat_full.datetime[0])
    
        if len(zscat_full) == 0:
            print(str(sid) + ' has no background measurement and thus skipping')
            continue

        # I need to use the same beamc or it can introduce a lot of error, here I calculate the median beamc
        median_beamc = []
        for pol_state in [0,90]:
            # select measurements with polarization
            zscat_pol = zscat_full[zscat_full.polarizationAngle==pol_state]
            samp_pol = samp_df[samp_df.polarizationAngle == pol_state]
            
            # TODO check ratio and if it changes
            # compute optical betransmission and beamc
            
            # PG if we multiply the perpendicular polarization ref (first meas) by 1.16
            # then it is much closer to matching the other pol on the ref and this
            # needs to be done on the zscat and on the pdat LREF
            # this is totally empirical
            corrFactor = 1
            if empirical_corr_factor & (pol_state == 90):
                corrFactor = 1.165
            LaserRatio = zscat_pol.LP.median()/(zscat_pol.LREF.median()*corrFactor)       # ratio of transmitted power / laser ref
            tau = samp_pol.LP.values/LaserRatio/(samp_pol.LREF.values*corrFactor) # compute optical transmission, taking drift in laser power into account
            # tau = (samp_pol.LP.values/zscat_pol.LP.median()) / (samp_pol.LREF.values/zscat_pol.LREF.median()) # compute optical transmission, taking drift in laser power into account
            # tau = samp_pol.LP.values / zscat_pol.LP.median()
            beamc = -(1/0.07)*np.log(tau)     # compute beamc for 70mm pathlength
            median_beamc.append(beamc)
        # print(median_beamc)
        
        median_beamc = (np.nanmedian(median_beamc[0]) + np.nanmedian(median_beamc[1]))/2
        if beamc_forced:
            median_beamc = beamc_forced
    
        for pol_state in [0,90]:
            # select measurements with polarization
            zscat_pol = zscat_full[zscat_full.polarizationAngle==pol_state]
            samp_pol = samp_df[samp_df.polarizationAngle == pol_state]
            
            # TODO check ratio and if it changes
            # compute optical betransmission and beamc
            
            # PG if we multiply the perpendicular polarization ref (first meas) by 1.16
            # then it is much closer to matching the other pol on the ref and this
            # needs to be done on the zscat and on the pdat LREF
            # this is totally empirical
            corrFactor = 1
            if empirical_corr_factor & (pol_state == 90):
                corrFactor = 1.165
    
            # TODO nanmedians?
            LaserRatio = zscat_pol.LP.median()/(zscat_pol.LREF.median()*corrFactor)       # ratio of transmitted power / laser ref
            tau = samp_pol.LP.values/LaserRatio/(samp_pol.LREF.values*corrFactor) # compute optical transmission, taking drift in laser power into account
            # tau = (samp_pol.LP.values/zscat_pol.LP.median()) / (samp_pol.LREF.values/zscat_pol.LREF.median()) # compute optical transmission, taking drift in laser power into account

            # tau = samp_pol.LP.values / zscat_pol.LP.median()
            beamc = -(1/0.07)*np.log(tau)     # compute beamc for 70mm pathlength
            # replace with median from both polarizations
            # beamc = median_beamc
            
            ###############################
            ######### FWD SCATTER #########
            ###############################
            
            # these are the parameters for the forward scattering
            rho = 1.18
            theta0air = (0.102/120)
            
            edge_angles = theta0air*rho**np.arange(0,37)
            edge_angles = np.arcsin(np.sin(edge_angles)/1.33)
            
            # find solid angle
            dOmega=np.cos(edge_angles[0:36])-np.cos(edge_angles[1:37]);  
            dOmega=dOmega*2*np.pi/6 # factor 6 takes care of rings covering only 1/6th circle
                        
            # correct forward scatter
            # todo I'm not currently accounting for changes in laser ref power because I don't trust them
            zscatFwd_scaled = np.nanmedian(np.vstack(zscat_pol.scatFwd),axis=0) * samp_pol.LREF.median() / zscat_pol.LREF.median() # correct background for changes in Lref between background and sample

            scatFwd = np.vstack(samp_pol.scatFwd)
            initial_scatFwd = scatFwd

            # correct particle scatter for attenuation
            scatFwd = scatFwd / (tau.reshape(-1,1)) 

            # subtract the background
            if subtract_zscat:
                scatFwd = scatFwd - zscatFwd_scaled 
            else:
                scatFwd = scatFwd # if we don't want to subtract the zscat

            scatFwd_tofilter = scatFwd

            # select dcal with polariztion that matches sample
            dcal_pol = None
            if pol_state == 0:
                dcal_pol = dcal0pol
            elif pol_state == 90:
                dcal_pol = dcal90pol

            # apply dcal
            cscatFwd = scatFwd * dcal_pol  
            # cscatFwd = scatFwd # don't apply dcal

            cscatFwd = cscatFwd / dOmega / 1000000
            if pol_state == 0:
                LREF_nominal = np.array([673])
            elif pol_state == 90:
                LREF_nominal = np.array([801])
            cscatFwd = cscatFwd * (LREF_nominal/samp_pol.LREF.values).reshape(-1,1)

            cscatFwd[scatFwd_tofilter<20] = np.nan

            ##########################################
            ###### now monte carlo the fwd scat ######
            ##########################################

            cscatFwd_mc_list = []

            # TODO maybe see if the zscat is > than cscat and then flag it out or something
            
            for i in [1]:# np.linspace(.98,1.02,3):
                # for rand in range(5):
                scatFwd = np.vstack(samp_pol.scatFwd)

                # correct particle scatter for attenuation
                scatFwd = scatFwd / (tau.reshape(-1,1)) 

                # subtract background
                scatFwd = scatFwd - zscatFwd_scaled * samp_pol.LREF.median() / zscat_pol.LREF.median() # correct background for changes in Lref between background and sample

                # add random values for sensitivity
                # scatFwd = scatFwd + np.random.uniform(low=-10, high=10, size=(scatFwd.shape))
    
                # select dcal with polariztion that matches sample
                dcal_pol = None
                if pol_state == 0:
                    dcal_pol = dcal0pol
                elif pol_state == 90:
                    dcal_pol = dcal90pol
    
                # apply dcal
                cscatFwd_mc = scatFwd * dcal_pol  
                # cscatFwd = scatFwd # don't apply dcal
    
                cscatFwd_mc = cscatFwd_mc / dOmega / 1000000
                cscatFwd_mc = cscatFwd_mc * (LREF_nominal/samp_pol.LREF.values).reshape(-1,1)

                cscatFwd_mc[scatFwd_tofilter<20] = np.nan
                cscatFwd_mc_list.append(cscatFwd_mc)

            

            ###############################
            ######### SSD SCATTER #########
            ###############################

            # parameters for ssd corrections
            
            # these are determined from 0.1 um beads on the perpendicular polarization (so that it is nearly isotropic)
            k_factors_perp = np.array([[
            1.61492324e-05, 1.83815677e-05, 2.05782186e-05, 2.32344816e-05,
            2.41606937e-05, 3.01886216e-05, 5.04243718e-05, 1.61539175e-04,
            4.88011978e-05, 3.37287778e-05, 2.88274429e-05, 2.96529573e-05,
            2.86983490e-05, 2.90227115e-05, 3.06281499e-05, 4.55926815e-05,
            1.32786626e-04, 4.66061385e-05, 2.73729258e-05, 2.22251387e-05,
            2.01312527e-05, 1.86630717e-05, 1.62262316e-05, 1.51668477e-05]])[0]

            k_factors_perp = (k_factors_perp + np.flip(k_factors_perp))/2

            # these are determine from 60nm beads on the perp polarization
            # k_factors_perp = np.array([1.47415135e-05, 1.67690011e-05, 1.89293080e-05, 2.15385948e-05,
            # 2.24640514e-05, 2.81337922e-05, 4.72764250e-05, 1.56067790e-04,
            # 4.58004889e-05, 3.13468901e-05, 2.68791865e-05, 2.76252171e-05,
            # 2.67599698e-05, 2.72017513e-05, 2.88672666e-05, 4.32221458e-05,
            # 1.29497157e-04, 4.48017343e-05, 2.63149051e-05, 2.12645341e-05,
            # 1.93149315e-05, 1.79160900e-05, 1.56123491e-05, 1.43750076e-05])

            # # these are with a minor ref correction
            # k_factors_perp = np.array([1.79404724e-05, 1.94247507e-05, 2.19028540e-05, 2.40085963e-05,
            #    2.57302296e-05, 3.18455847e-05, 5.36420529e-05, 1.61739567e-04,
            #    5.21033169e-05, 3.53802746e-05, 3.18635047e-05, 3.21084104e-05,
            #    3.21084104e-05, 3.18635047e-05, 3.53802746e-05, 5.21033169e-05,
            #    1.61739567e-04, 5.36420529e-05, 3.18455847e-05, 2.57302296e-05,
            #    2.40085963e-05, 2.19028540e-05, 1.94247507e-05, 1.79404724e-05])

            # k_factors_perp = np.array([2.03143673e-05, 2.23781330e-05, 2.45124115e-05, 2.77917370e-05,
            #    2.86784865e-05, 3.60363637e-05, 6.07089650e-05, 1.96200188e-04,
            #    5.58144775e-05, 3.86358260e-05, 3.27605043e-05, 3.42385070e-05,
            #    3.30489783e-05, 3.41401020e-05, 3.58082274e-05, 5.44004759e-05,
            #    1.59236064e-04, 5.57060720e-05, 3.29011747e-05, 2.70182479e-05,
            #    2.44220907e-05, 2.30333882e-05, 2.01645880e-05, 1.92409883e-05])

            # this option would basically take the average of the flipped version to make it symmetrical
            # k_factors_perp = np.mean([k_factors_perp,np.flip(k_factors_perp)],axis=0)
            # for now we set them equal on both polarizations
            # k_factors_perp = 1
            k_factors_para = k_factors_perp

            # this is the array mismatch correction which I've just created empirically
            ref_corr = np.zeros(24)
            ref_corr[7] = 1.0
            ref_corr[6] = 0.3
            ref_corr[5] = 0.15
            ref_corr[4] = 0.1
            
            ref_corr_back = np.zeros(24)
            ref_corr_back[22] = 0.025
            ref_corr_back[21] = 0.05
            ref_corr_back[20] = 0.1
            ref_corr_back[19] = 0.15
            ref_corr_back[18] = 0.2
            ref_corr_back[17] = 0.32
            ref_corr_back[16] = 1

            ref_corr = np.flip(ref_corr_back)

            # updates as of Jan 23rd 2026
           #  ref_corr = np.array([-0.84800071, -0.88489135, -0.80628523, -0.64113914, -0.42291661,
           # -0.2002593 ,  0.05965305,  0.83944931, -0.25015213, -0.16697565,
           # -0.12047388, -0.07965578,  0.        ,  0.        ,  0.        ,
           #  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
           #  0.        ,  0.        ,  0.        ,  0.        ])

           #  ref_corr_back = np.array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
           #  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
           #  0.        ,  0.        , -0.10122784, -0.09674394, -0.1420754 ,
           # -0.17309439,  0.93415145,  0.31278274,  0.06645415, -0.2084256 ,
           # -0.43949708, -0.5770984 , -0.63203992, -0.75381941])

            # these are pathlengths in meters for the light reflected off the recieve lens from Kirby
            # we then subtract them from the normal pathlengths to get just the additional pathlength which is basically 7cm every time since we're traveling the whole chamber and then reflecting
            mirror_attenuation_pathlengths = np.array([169,154, 143, 135, 128, 121, 113, 106, 157, 150, 145, 141, 138, 135, 131, 126, 148, 146, 144, 142, 140, 139, 137, 136])/1000 - np.flip(pathlength_ssd)

            # transmission_factor = 0.25*0.25 # from the ND filter because it is an OD of 0.6 and we pass through it twice
            # transmission_factor = 0 # because none comes back from the felt
            # transmission_factor = 1
            
            if pol_state == 0: # unsure if these should be 16% or 16.5%
                mirror_ref = (0.165)*transmission_factor # this is an assumption about the reflection factor off the recieve lens
            elif pol_state == 90:
                mirror_ref = (0.16)*transmission_factor # this is an assumption about the reflection factor off the recieve lens
            
            # do we think this is appropriate, really it depends on the near fwd VSF
            # this basically says if the detector is wider then we attenuate the light less because it is seeing more of the near fwd scattered light
            # atten_factor = widths_ssd
            atten_factor = 1

            zscatSsd_scaled = np.nanmedian(np.vstack(zscat_pol.scatSsd),axis=0)# * (samp_pol.LREF.median() / zscat_pol.LREF.median()) # correct background for changes in Lref between background and sample
            # TODO I should really change this to either the average value or not do it but I'm finding it is set to nan when sometimes the zscat is extremely clean so it is basically zero
            # zscatSsd_scaled[np.isnan(zscatSsd_scaled)] = 2
            if pol_state == 0:
                zscatSsd_scaled[np.isnan(zscatSsd_scaled)] = 6.5
            elif pol_state == 90:
                zscatSsd_scaled[np.isnan(zscatSsd_scaled)] = 1.6

            ###############################################
            ###### start processing the side scatter ######
            ###############################################

            scatSsd = np.vstack(samp_pol.scatSsd)

            plot = False
            
            # zscatSsd_scaled[np.isnan(zscatSsd_scaled)] = 2

            if plot:
                plt.plot(angles_ssd,scatSsd.T,c='k')
                plt.plot(angles_ssd,zscatSsd_scaled.T,c='grey')
                plt.yscale('log')
                plt.title('pol state = ' + str(pol_state)+'before zscat sub')
                plt.ylim(30,8000)
                plt.ylabel('uncalibrated scat value')
                plt.xlabel('scattering angle [°]')
                plt.show()
                
        
            # subtract the zscat            
            if subtract_zscat:
                scatSsd = scatSsd - zscatSsd_scaled        # subtract the background

            if plot:
                plt.plot(angles_ssd,scatSsd.T,c='k')
                plt.yscale('log')
                plt.title('pol state = ' + str(pol_state)+'after zscat sub')
                plt.ylim(30,8000)
                plt.ylabel('uncalibrated scat value')
                plt.xlabel('scattering angle [°]')
                plt.show()

            # now apply the k factors
            cscatSsdUncal = scatSsd*k_factors_perp
            cscatSsdUncal = cscatSsdUncal * (LREF_nominal/samp_pol.LREF.values).reshape(-1,1)

            if plot:
                print('with k factors then before and after mirror')
                plt.plot(angles_ssd,cscatSsdUncal.T,c='k')
                plt.yscale('log')
                plt.title('pol state = ' + str(pol_state) +'after k factors')
                plt.ylim(0.001,.1)
                plt.ylabel('calibrated scat value')
                plt.xlabel('scattering angle [°]')
                plt.show

            cscatSsdUncal = cscatSsdUncal - ref_corr*(cscatSsdUncal[:,7]-cscatSsdUncal[:,8]).reshape(-1,1) -  ref_corr_back*(cscatSsdUncal[:,16]-cscatSsdUncal[:,15]).reshape(-1,1)
            
            # now do the mirror reflection correction
            intermediate_mirror_correction = cscatSsdUncal-np.flip(cscatSsdUncal,axis=1) *(mirror_ref) /np.exp(mirror_attenuation_pathlengths*beamc.reshape(-1,1)/atten_factor)
            cscatSsdUncal = cscatSsdUncal-np.flip(intermediate_mirror_correction,axis=1) *(mirror_ref) /np.exp(mirror_attenuation_pathlengths*beamc.reshape(-1,1)/atten_factor)

            cscatSsdMirror = np.flip(intermediate_mirror_correction,axis=1) *(mirror_ref) /np.exp(mirror_attenuation_pathlengths*beamc.reshape(-1,1) / atten_factor)
            
            if plot:
                plt.plot(angles_ssd,cscatSsdUncal.T,c='k')
                plt.yscale('log')
                plt.title('pol state = ' + str(pol_state) +'after mirror sub')
                plt.ylim(0.001,.1)
                plt.ylabel('calibrated scat value')
                plt.xlabel('scattering angle [°]')
                plt.show()

            # now do the attenuation correction
            attenCorr = 1 # if we don't want to add any correction
            # TODO we know this value is often bad so consider not doing it?
            if apply_atten_corr:
                attenCorr = np.exp((beamc.reshape(-1,1)/atten_factor*pathlength_ssd.reshape(-1,24))) # calculate attenuation correction for side scatter detectors        

            cscatSsdUncal = cscatSsdUncal * attenCorr  # apply attenuation correction

            if plot:
                plt.plot(angles_ssd,cscatSsdUncal.T,c='k')
                plt.yscale('log')
                plt.title('pol state = ' + str(pol_state) +'after attenuation correction')
                plt.ylim(0.001,.1)
                plt.ylabel('calibrated scat value')
                plt.xlabel('scattering angle [°]')
                plt.show()

            # now correct for the array mismatch
            cscatSsdUncal = cscatSsdUncal - ref_corr*(cscatSsdUncal[:,7]-cscatSsdUncal[:,8]).reshape(-1,1) -  ref_corr_back*(cscatSsdUncal[:,16]-cscatSsdUncal[:,15]).reshape(-1,1)    

            if plot:
                plt.plot(angles_ssd,cscatSsdUncal.T,c='k')
                plt.yscale('log')
                plt.title('pol state = ' + str(pol_state) +'after array mismatch correction')
                plt.ylim(0.001,.1)
                plt.ylabel('calibrated scat value')
                plt.xlabel('scattering angle [°]')
                plt.show()

            # full equation = 
            #              ((scatSsd - zscatSsd_scaled)*k_factors_perp - np.flip((scatSsd - zscatSsd_scaled)*k_factors_perp) *(mirror_ref) /np.exp(mirror_attenuation_pathlengths*beamc.reshape(-1,1)/4)) * attenCorr

            ########################################
            ######### Finalize Data Format #########
            ########################################
            
            # add data to existing table
            samp_pol['tau'] = list(tau)
            samp_pol['initial_scatFwd'] = list(initial_scatFwd)
            samp_pol['scatFwd_tofilter'] = list(scatFwd_tofilter)
            samp_pol['beamc'] = list(beamc)
            samp_pol['scatFwd'] = list(scatFwd)
            samp_pol['cscatFwd'] = list(cscatFwd)
            samp_pol['cscatSsdUncal'] = list(cscatSsdUncal)
            samp_pol['cscatSsdMirror'] = list(cscatSsdMirror)
            samp_pol['zscatSsd_scaled'] = list(zscatSsd_scaled.reshape(1,-1).repeat(len(tau),axis=0))
            samp_pol['zscatFwd_scaled'] = list(zscatFwd_scaled.reshape(1,-1).repeat(len(tau),axis=0))
            
            ##########################################
            ## now we'll do the Monte Carlo version ##
            ##########################################

            num_in_dist = 3

            distribution = np.linspace(.75,1.25,num_in_dist)
            zscat_dist = np.array([i*distribution for i in zscatSsd_scaled]).T

            distribution = np.linspace(.85,1.15,num_in_dist)
            kfactors_dist = np.array([i*distribution for i in k_factors_perp]).T

            # here we're just allowing for a small amout of nonlinearity in the detector response in case it is close to maxing out
            # this leads to an underestimate of the true value and thus we subtract too little
            # nonlinearity_dist = np.linspace(1,1.01,num_in_dist) 
            nonlinearity_dist = [1]

            distribution = np.linspace(.95,1.05,num_in_dist)
            ref_corr_dist = np.array([i*distribution for i in ref_corr]).T
            ref_corr_dist = [ref_corr]

            distribution = np.linspace(.95,1.05,num_in_dist)
            ref_corr_back_dist = np.array([i*distribution for i in ref_corr_back]).T
            ref_corr_back_dist = [ref_corr_back]

            distribution =  np.linspace(.8,1.2,num_in_dist)
            # distribution = np.array([1])
            mirror_ref_dist = distribution*mirror_ref

            # distribution = np.linspace(1,2,num_in_dist)
            # distribution = np.logspace(.01,1,5)
            # distribution = np.logspace(-.5,.05,5)
            distribution = np.linspace(.5,2,5)
            atten_factor = widths_ssd
            # atten_factor = 1
            atten_factor_dist = np.array([i*distribution for i in atten_factor]).T
            # atten_factor_dist = np.array([atten_factor*1.5])
            atten_factor_dist = [1]

            # TODO revisit this
            # beamc = np.array([2.5]) # this is just if I need to specify the beamc
            beamc = median_beamc
            distribution = np.linspace(.85,1.15,num_in_dist)
            beamc_dist = np.array([i*distribution for i in beamc.reshape(-1,1)]).T

            modify_nonnearfwd = np.array([1])# np.linspace(.8,1.2,num_in_dist)

            full_cscat_dist = []

            first_ref_corr_idx = 0

            full_vsf_intg_list = []

            # TODO the nonnearforward should be based on the VSF or phase function, not on the polarization specific function

            for zscat in zscat_dist:
                # for kfactors in kfactors_dist:
                # for nonlinear_meas in nonlinearity_dist:
                for mod_nnf in modify_nonnearfwd:
                    # for first_ref_corr_idx in range(num_in_dist):
                        for mirror in mirror_ref_dist:
                            for atten_factor in atten_factor_dist:
                                for beamc in beamc_dist:
                                    scatSsd = np.vstack(samp_pol.scatSsd)
                                    
                                    if subtract_zscat:
                                        # subtract the background
                                        scatSsd = scatSsd - zscat 

                                    # now apply the k factors
                                    # cscatSsdUncal = scatSsd*kfactors
                                    cscatSsdUncal = scatSsd*k_factors_perp
                                    cscatSsdUncal = cscatSsdUncal * (LREF_nominal/samp_pol.LREF.values).reshape(-1,1)
                                    
                                    # do the array mismatch correction
                                    cscatSsdUncal = cscatSsdUncal - ref_corr_dist[first_ref_corr_idx]*(cscatSsdUncal[:,7]-cscatSsdUncal[:,8]).reshape(-1,1) -  ref_corr_back_dist[first_ref_corr_idx]*(cscatSsdUncal[:,16]-cscatSsdUncal[:,15]).reshape(-1,1)
                                    # TODO figure out some way if there are nans to save the remainder of the data
        
                                    # now do the mirror reflection correction
                                    
                                    intermediate_mirror_correction = cscatSsdUncal-np.flip(cscatSsdUncal,axis=1) *(mirror) /np.exp(mirror_attenuation_pathlengths*beamc.reshape(-1,1))
                                    first_guess_at_shape = cscatSsdUncal-np.flip(intermediate_mirror_correction,axis=1) *(mirror) /np.exp(mirror_attenuation_pathlengths*beamc.reshape(-1,1))
                                    nonnear_fwd_ratio = np.nan
                                    if mirror == 0:
                                        first_guess_at_shape = cscatSsdUncal
                                    ####
                                    try:
                                        ssd_meas = np.nanmedian(first_guess_at_shape,axis=0)
                                        vsf_fwd_med = np.nanmedian(np.array(cscatFwd_mc_list).reshape(-1,36),axis=0) / np.nanmedian(beamc)
                                        
                                        full_angles = np.hstack([angles_fwd,angles_ssd])
                                        full_vsf = np.hstack([vsf_fwd_med,ssd_meas])
                                        
                                        # unfortunately basically none of these detectors are reliable at angles less than 1.5 so we don't use them here unless they are over a threshold
                                        full_vsf[:22][full_vsf[:22] < 10] = np.nan
                                        full_vsf[:15] = np.nan
                                        
                                        # nearest neighbor interpolation for all nan values
                                        mask = np.isnan(full_vsf)
                                        full_vsf[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), full_vsf[~mask])
                                                                            
                                        # now extrapolate it to 0.01°
                                        f = interpolate.interp1d(np.log10(full_angles), np.log10(full_vsf), fill_value='extrapolate',kind='linear')
                                        interpolated_vsf = f(np.log10(np.linspace(0.01,150,1000)))
                                    
                                        # now do the integrations of the whole VSF and the amount of scattering at angles <3°
                                        full_vsf_intg = np.trapz(10**interpolated_vsf, x=np.deg2rad(np.linspace(0.01,150,1000)))
                                        vsf_3deg_intg = np.trapz((10**interpolated_vsf)[:21], x=np.deg2rad(np.linspace(0.01,150,1000))[:21])

                                        full_vsf_intg_list.append(full_vsf_intg/2/np.pi)
    
                                        # this is the best guess at the relevant beam-c correction factor given the detector widths
                                        nonnear_fwd_ratio = 1 - (vsf_3deg_intg / full_vsf_intg)
                                    except Exception:
                                        # Did not have any fwd data we assume assigning a fwd ratio of 0.2 for most natural samples
                                        nonnear_fwd_ratio = nnf_factor

                                    if np.isnan(nonnear_fwd_ratio):
                                        nonnear_fwd_ratio = nnf_factor
                                    # nonnear_fwd_ratio = 0.3
                                    # TODO this isn't quite right...
                                    if mirror == 0:
                                        pass
                                    else:
                                        intermediate_mirror_correction = cscatSsdUncal-np.flip(cscatSsdUncal,axis=1) *(mirror) /np.exp(mirror_attenuation_pathlengths*nonnear_fwd_ratio*mod_nnf*beamc.reshape(-1,1))
                                        cscatSsdUncal = cscatSsdUncal-np.flip(cscatSsdUncal,axis=1) *(mirror) /np.exp(mirror_attenuation_pathlengths*nonnear_fwd_ratio*mod_nnf*beamc.reshape(-1,1))

                                    ####
                    
                                    # now do the attenuation correction
                                    attenCorr = 1 # if we don't want to add any correction
                                    if apply_atten_corr:
                                        attenCorr = np.exp((beamc.reshape(-1,1)*nonnear_fwd_ratio*mod_nnf*pathlength_ssd.reshape(-1,24))) # calculate attenuation correction for side scatter detectors        
                                    # apply attenuation correction
                                    cscatSsdUncal = cscatSsdUncal * attenCorr  
                        
                                    # now correct for the array mismatch again
                                    cscatSsdUncal = cscatSsdUncal - ref_corr_dist[first_ref_corr_idx]*(cscatSsdUncal[:,7]-cscatSsdUncal[:,8]).reshape(-1,1) -  ref_corr_back_dist[first_ref_corr_idx]*(cscatSsdUncal[:,16]-cscatSsdUncal[:,15]).reshape(-1,1)
    
                                    full_cscat_dist.append(cscatSsdUncal)
            # print(nonnear_fwd_ratio)
            samp_pol['nonnear_fwd_ratio'] = nonnear_fwd_ratio
            samp_pol['vsf_intg'] = np.nanmedian(full_vsf_intg_list)
                    
            if pol_state == 0:
                monte_carlo_results_0.append(full_cscat_dist)
            elif pol_state == 90:
                monte_carlo_results_90.append(full_cscat_dist)
            cscatFwd_monte_carlo.append(cscatFwd_mc_list)
            processed_samps.append(samp_pol)
            
    processed = pd.concat(processed_samps, ignore_index=True)

    return(processed, monte_carlo_results_0, monte_carlo_results_90, cscatFwd_monte_carlo)

def full_process(csv_files, using_fscat=False, apply_atten_corr=True, subtract_zscat=True, empirical_corr_factor=True, specific_zscat_file=False, transmission_factor=1, nnf_factor=0.2, beamc_forced=None):
    raw = open_raw_files(csv_files)
    loaded = hrzn_load(raw)
    
    if specific_zscat_file:
        raw_zscat = open_raw_files(specific_zscat_file)
        loaded_zscat = hrzn_load(raw_zscat)
    else:
        loaded_zscat=[]

    processed = hrzn_process(loaded, using_fscat=using_fscat, apply_atten_corr=apply_atten_corr, subtract_zscat=subtract_zscat, empirical_corr_factor=empirical_corr_factor, specific_zscat_file=loaded_zscat, transmission_factor=transmission_factor,nnf_factor=nnf_factor, beamc_forced=beamc_forced)

    return(processed)



def empirical_calibration():
    pass