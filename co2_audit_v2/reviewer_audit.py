"""Reproduce the four reviewer checks; never infer missing PBL/CIA metadata.

Flux tables retain both incident and net radiation. Heights are fixed reference
heights above the model surface, NOT a diagnosis of boundary-layer depth.
Temperature means distinguish moving an averaging boundary from rerunning FDH.
"""
import csv,json,hashlib
from pathlib import Path
import numpy as np
import xarray as xr
from radiation import ROOT,RDATA,Radiation,diurnal_quadrature

OUT=ROOT/'results/reviewer_audit'
OUT.mkdir(parents=True,exist_ok=True)

def dump_csv(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    base=ROOT/'results/global_tsis'
    a=xr.load_dataset(base/'reference_profiles.nc',decode_cf=False)
    w=np.asarray(a.profile_weight);p=np.asarray(a.pres_level);pl=np.asarray(a.pres_layer)
    with np.load(base/'adjusted_full_280_420_diagnostics.npz') as d:
        full=d['delta_T'].copy();trop=d['tropopause_index'].copy();z=d['height_level_m'].copy()
    with np.load(base/'adjusted_co2off_280_420_diagnostics.npz') as d:off=d['delta_T'].copy()
    idx=np.arange(a.sizes['layer'])[None,:]
    def mean(F,O,tt,exact=False):
        dp=np.maximum(0,p[:,1:]-np.maximum(p[:,:-1],100.)) if exact else np.diff(p,axis=1)*(pl>=100)
        mass=w[:,None]*dp*(idx<tt[:,None])
        mf=float((mass*F).sum()/mass.sum());mo=float((mass*O).sum()/mass.sum())
        return dict(full_K=mf,no_CO2_SW_K=mo,SW_offset_K=mf-mo,cooling_reduction_percent=100*(mf-mo)/(-mo),weighted_pressure_mass_Pa=float(mass.sum()))
    temps=[]
    for shift,label in [(-1,'one_level_higher'),(0,'baseline'),(1,'one_level_lower')]:
        tt=trop+shift
        temps.append(dict(test='averaging_boundary_only',boundary=label,top='midpoint_1hPa',**mean(full,off,tt)))
        if shift:
            dpath=ROOT/'results/fdh_sensitivities'
            with np.load(dpath/f'tropopause_{label}_full_diagnostics.npz') as d:F=d['delta_T'].copy();t=d['tropopause_index'].copy()
            with np.load(dpath/f'tropopause_{label}_co2off_diagnostics.npz') as d:O=d['delta_T'].copy()
            assert np.array_equal(t,tt)
            temps.append(dict(test='FDH_recomputed',boundary=label,top='midpoint_1hPa',**mean(F,O,tt)))
    temps.append(dict(test='exact_upper_boundary',boundary='baseline',top='partial_layer_exact_100Pa',**mean(full,off,trop,True)))
    dump_csv('temperature_boundary_sensitivity.csv',temps)
    boundaries=[]
    for i,t in enumerate(trop):
        top=int(np.flatnonzero(pl[i]>=100)[0])
        boundaries.append(dict(site=i,latitude_deg=float(a.lat[i]),weight=float(w[i]),interface_index_zero_based=int(t),
            tropopause_hPa=float(p[i,t]/100),tropopause_km=float(z[i,t]/1000),
            one_level_higher_hPa=float(p[i,t-1]/100),one_level_higher_km=float(z[i,t-1]/1000),
            one_level_lower_hPa=float(p[i,t+1]/100),one_level_lower_km=float(z[i,t+1]/1000),
            first_included_layer_index=top,actual_upper_interface_hPa=float(p[i,top]/100),
            last_included_layer_index=int(t-1)))
    dump_csv('tropopause_by_profile.csv',boundaries)
    mu,qw=diurnal_quadrature(a,8);gx,gw=np.polynomial.legendre.leggauss(8)
    geometry=[dict(site=i,node=k,latitude_deg=float(a.lat[i]),day_since_2014_01_01=float(a.time[i]),
        profile_weight=float(w[i]),TSI_Wm2=float(a.total_solar_irradiance[i]),gauss_node=float(gx[k]),gauss_weight=float(gw[k]),
        mu0=float(mu[k,i]),daily_weight=float(qw[k,i]),plane_parallel_airmass=float(1/mu[k,i]) if mu[k,i]>0 else '')
        for i in range(len(w)) for k in range(8)]
    dump_csv('solar_geometry_nodes.csv',geometry)
    m=Radiation(solar='tsis')
    swdata=xr.load_dataset(RDATA/'rrtmgp-gas-sw-g224.nc')
    bounds=np.asarray(swdata.bnd_limits_wavenumber)
    # Inspect sensitivity rather than assuming that a density-scaled O2 term
    # includes CO2-O2 pair absorption. This band overlaps the official pair data.
    b=9;g0,g1=np.asarray(swdata.bnd_limits_gpt)[b].astype(int)
    op0=m.sw_optics(a,280);op1=m.sw_optics(a,560)
    tau_difference=float(abs(op1.tau.isel(gpt=slice(g0-1,g1))-op0.tau.isel(gpt=slice(g0-1,g1))).max())
    spec=[];spectral_checks={}
    for nmu in [4,8,16]:
        path=OUT/f'spectral_delta_nmu{nmu}.nc'
        if nmu==8:d=xr.load_dataset(ROOT/'results/scenarios/global_spectral_delta.nc')
        elif path.exists():d=xr.load_dataset(path)
        else:
            print('Computing band radiation, Gauss nodes:',nmu,flush=True)
            d=m.sw_flux(a,420,'daily',nmu,band=True)-m.sw_flux(a,280,'daily',nmu,band=True)
            d.to_netcdf(path)
        dn=np.asarray(d.band_sw_flux_down-d.band_sw_flux_up)
        for band in range(14):
            spec.append(dict(nmu=nmu,band_zero_based=band,wavenumber_min_cm1=float(bounds[band,0]),wavenumber_max_cm1=float(bounds[band,1]),
                wavelength_min_um=1e4/float(bounds[band,1]),wavelength_max_um=1e4/float(bounds[band,0]),
                SW_net_surface_Wm2=float(dn[band,:,-1]@w),SW_net_TOA_Wm2=float(dn[band,:,0]@w),
                SW_net_tropopause_Wm2=float(dn[band,np.arange(len(w)),trop]@w)))
        assert float(abs(d.band_sw_flux_down.sum('band')-d.sw_flux_down).max())<1e-9
        # Independently stored broadband runs check that cached band outputs
        # use the same atmosphere, concentrations, source and angular rule.
        original=base if nmu==8 else ROOT/f'results/convergence_nmu{nmu}'
        original_delta=xr.load_dataset(original/'instant_420.nc')-xr.load_dataset(original/'instant_280.nc')
        error=max(float(abs(d[key]-original_delta[key]).max()) for key in ['sw_flux_down','sw_flux_up'])
        assert error<1e-9,(nmu,error)
        spectral_checks[f'nmu{nmu}_max_band_run_minus_independent_broadband_Wm2']=error
    dump_csv('band_angular_convergence.csv',spec)
    flux=[];column_flux=[];max_lw_identity=0
    def load_flux(name):return xr.load_dataset(base/name)
    def interp_height(arr,h):
        arr=np.asarray(arr)
        return np.array([np.interp(h,z[i,::-1],arr[i,::-1]) for i in range(len(w))])
    f0=load_flux('instant_280.nc')
    for c1 in [420,560,1120]:
        instant=load_flux(f'instant_{c1}.nc')-f0
        fd=load_flux(f'adjusted_full_280_{c1}_flux.nc')-load_flux(f'adjusted_full_280_{c1}_control.nc')
        swonly=fd-(load_flux(f'adjusted_co2off_280_{c1}_flux.nc')-load_flux(f'adjusted_co2off_280_{c1}_control.nc'))
        for experiment,delta in [('instant_full',instant),('FDH_full',fd),('SW_induced_B_minus_C',swonly)]:
            err=np.asarray(delta.lw_flux_down-delta.lw_flux_up)[:,-1]-np.asarray(a.surface_emissivity)*np.asarray(delta.lw_flux_down)[:,-1]
            max_lw_identity=max(max_lw_identity,float(abs(err).max()))
            for h in [0,500,1000,2000]:
                vals={key:interp_height(delta[key],h) for key in ['sw_flux_down','sw_flux_up','lw_flux_down','lw_flux_up']}
                vals['SW_net']=vals['sw_flux_down']-vals['sw_flux_up'];vals['LW_net']=vals['lw_flux_down']-vals['lw_flux_up'];vals['total_net']=vals['SW_net']+vals['LW_net']
                means={k+'_Wm2':float(v@w) for k,v in vals.items()}
                flux.append(dict(c0=280,c1=c1,experiment=experiment,reference_height_m=h,**means,
                    SW_down_loss_over_LW_down_gain_percent=-100*means['sw_flux_down_Wm2']/means['lw_flux_down_Wm2'],
                    SW_net_loss_over_LW_net_gain_percent=-100*means['SW_net_Wm2']/means['LW_net_Wm2']))
                for i in range(len(w)):
                    column_flux.append(dict(site=i,weight=float(w[i]),c0=280,c1=c1,experiment=experiment,reference_height_m=h,**{k+'_Wm2':float(v[i]) for k,v in vals.items()}))
    dump_csv('flux_decomposition.csv',flux);dump_csv('flux_decomposition_by_profile.csv',column_flux)
    # Surface optical properties enter the library at mixed float precision;
    # subtracting two large absolute fluxes leaves sub-micro-W/m2 residuals.
    assert max_lw_identity<2e-6
    main_rows=json.loads((base/'summary.json').read_text());flux_summary_error=0.
    for row in flux:
        if row['reference_height_m']!=0:continue
        kind={'instant_full':'instant','FDH_full':'adjusted_full','SW_induced_B_minus_C':'SW_induced_vs_co2off'}[row['experiment']]
        old=next(x for x in main_rows if (x['c0'],x['c1'],x['kind'])==(280,row['c1'],kind))
        for new_key,old_key in [('SW_net_Wm2','SW_surface'),('lw_flux_down_Wm2','LW_surface_down'),('total_net_Wm2','surface_total_radiative')]:
            flux_summary_error=max(flux_summary_error,abs(row[new_key]-old[old_key]))
    assert flux_summary_error<1e-9
    metadata={}
    for name in ['rrtmgp-gas-sw-g224.nc','rrtmgp-gas-lw-g256.nc']:
        file=RDATA/name;data=xr.load_dataset(file)
        metadata[name]=dict(sha256=hashlib.sha256(file.read_bytes()).hexdigest(),attrs=data.attrs,
            minor_identifiers=[str(s) for s in data.identifier_minor.values])
    out=dict(spectroscopy=dict(files=metadata,documented_model_lineage='LBLRTM 12.8; aer_v_3.6, largely HITRAN2012; MT_CKD3.2',
        exact_generator_TAPE5_available=False,exact_line_wing_cutoff_verified=False,CO2_N2_and_CO2_O2_CIA_verified=False,
        CO2_280_to_560_max_abs_delta_tau_band_12850_16000_cm1=tau_difference,
        CIA_caveat='No explicit CO2 pair term in metadata; zero CO2 dependence in one overlapping band excludes full representation of the O2-CO2 catalog range. No numerical bound on missing pair CIA is claimed.'),
        geometry=dict(nmu=8,mu_min=float(mu[mu>0].min()),mu_max=float(mu.max()),
            latitude_min=float(a.lat.min()),latitude_max=float(a.lat.max()),
            weighted_TOA_incident_Wm2=float((np.sum(mu*qw,axis=0)*np.asarray(a.total_solar_irradiance))@w),
            spectrum_file='final-v2_hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.nc',
            spectrum_grid_nm=1,within_band_source='native RRTMGP g-point fractions retained',
            spatial_integration='RFMIP profile_weight; no independent latitude quadrature'),
        tropopause=dict(min_hPa=min(r['tropopause_hPa'] for r in boundaries),max_hPa=max(r['tropopause_hPa'] for r in boundaries),
            min_km=min(r['tropopause_km'] for r in boundaries),max_km=max(r['tropopause_km'] for r in boundaries),
            weighted_hPa=sum(r['weight']*r['tropopause_hPa'] for r in boundaries),weighted_km=sum(r['weight']*r['tropopause_km'] for r in boundaries),
            actual_top_interface_hPa_min=min(r['actual_upper_interface_hPa'] for r in boundaries),actual_top_interface_hPa_max=max(r['actual_upper_interface_hPa'] for r in boundaries)),
        temperature=temps,band_angular_convergence=spec,flux=flux,
        verification=dict(max_surface_LW_net_minus_emissivity_times_LW_down_Wm2=max_lw_identity,
            max_surface_flux_difference_from_main_summary_Wm2=flux_summary_error,**spectral_checks),
        missing_inputs=dict(true_boundary_layer_top='RFMIP profiles contain no diagnosed PBL depth, winds or turbulent fluxes',
            provided_reference_heights_m=[500,1000,2000]),surface_emissivity=float(a.surface_emissivity.min()))
    (OUT/'summary.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(dict(temperature=temps,tropopause=out['tropopause'],geometry=out['geometry'],verification=out['verification']),indent=2),flush=True)

if __name__=='__main__':main()
