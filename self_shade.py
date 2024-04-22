import pvlib
import pandas as pd
import numpy as np

def shade_fractions(fs_array, eff_row_side_num_mods):
    """
    Shade fractions on each _course_ of a rack or tracker.

    Parameters
    ----------
    fs_array : numeric
        Scalar or vector of shade fractions for the entire rack or tracker. Zero (0)
        is unshaded and one (1) is fully shaded.
    eff_row_side_num_mods : int
        Number of courses in the rack as modules. EG: a 2P tracker has 2 courses.

    Returns
    -------
    Array with the shade fraction on each course.
    """
    fs_course = np.clip([
        fs_array * eff_row_side_num_mods - course 
        for course in range(eff_row_side_num_mods)], 0, 1)
    return fs_course

def non_linear_shade(n_cells_up, fs, fd):
    """
    Simple non-linear shade model.

    Assume shade loss is linear as direct shade moves from bottom through
    first cell, and then only diffuse for remainder of module up to the top.
    EG: If there are 10 cells, and ``fs`` is 0.05, the bottom cell is only
    half shaded, then the loss would be 50% of the direct irradiance. If
    80% of POA global on module is direct irradiance, IE: ``fd = 0.2``,
    then loss would be 40%. When the direct shade line reaches 10%, one
    cell is completely shaded, then the loss is 80%, and there's only diffuse
    light on the module. Any direct shade above the 1st cell has the same loss.

    Parameters
    ----------
    n_cells_up : int
        Number of cells vertically up
    fs : float
        Fraction of shade on module, 1 is fully shaded
    fd : numeric
        Diffuse fraction

    Returns
    -------
    Array of shade loss same size as ``fd``
    """
    pnorm = np.where(fs < 1/n_cells_up, 1 - (1 - fd)*fs*n_cells_up, fd)
    shade_loss = 1 - pnorm
    return shade_loss

def plant_power_with_shade_losses(
    resource_data,
    latitude,
    longitude,
    mount_type,
    gcr,
    dc_capacity_plant,
    power_plant_ac_max,
    dc_loss_fraction,
    gamma_pdc,
    shade_loss_model,
    default_site_transposition_model='haydavies',
    backtrack=True,
    backtrack_fraction=pd.NA,
    max_tracker_angle=pd.NA,
    axis_tilt=pd.NA,
    axis_azimuth=pd.NA,
    fixed_tilt=pd.NA,
    fixed_azimuth=pd.NA,
    n_cells_up=12,
    row_side_num_mods=pd.NA,
    row_height_center=pd.NA,
    row_pitch=pd.NA,
    bifacial=False,
    bifaciality=0.8,
    surface_tilt_timeseries=pd.Series([]),
    surface_azimuth_timeseries=pd.Series([]),
    use_measured_poa=False,
    use_measured_temp_module=False,
    **kwargs,
    ):

    """
    Power model for PV plants.

    Parameters
    ----------
    resource_data : pandas.DataFrame
        timeseries weather/resource data with the same format as is returned by 
        pvlib.iotools.get* functions
    [long list of additional arguments defining a plant based on [1]]
    surface_tilt_timeseries : pandas.DataFrame
        (optional) custom timeseries of the angle between the panel surface and 
        the earth surface, accounting for panel rotation. [degrees]
    surface_azimuth_timeseries : pandas.DataFrame
        (optional) custom timeseries of the azimuth of the rotated panel, 
        determined by projecting the vector normal to the panel's surface to 
        the earth's surface. [degrees]
    use_measured_poa : bool, default False
        If True, used measure POA data from ``resource_data`` (must have 
        column name 'poa').
    use_measured_temp_module: bool, default False
        If True, use measured back of module temperature from ``resource_data``
        (must have column name 'temp_module') in place of modeled cell
        temperature. 

    Returns
    -------
    power_ac : pandas.Series
        AC power. Same units as ``dc_capacity_plant`` and ``power_plant_ac_max``
        (ideally kW).
    fs_array : pandas.Series
        Shaded fraction for all courses on each row in the array.

    References
    ----------
    .. [1] William Hobbs, pv-plant-specification-rev3.0.csv, 
       https://github.com/williamhobbs/pv-plant-specifications
    """

    # Fill in some necessary variales with defaults if there is no value provided
    if pd.isna(row_side_num_mods):
        row_side_num_mods = 1 # default if no value provided

    if pd.isna(row_height_center):
        row_height_center = 1 # default if no value provided

    if pd.isna(row_pitch):
        row_pitch = 2 / gcr # default if no value provided
    
    if pd.isna(backtrack_fraction) and backtrack==True:
        backtrack_fraction = 1.0 # default if no value provided

    if backtrack_fraction==0:
        backtrack=False # switch to truetracking to avoid divide by zero 

    eta_inv_nom = 0.98

    times = resource_data.index 
    loc = pvlib.location.Location(latitude=latitude, longitude=longitude, tz=times.tz)
    solar_position = loc.get_solarposition(times)

    # surface tilt and azimuth
    if surface_tilt_timeseries.empty | surface_azimuth_timeseries.empty:
        if mount_type == 'single-axis':
            # modify tracker gcr if needed
            if backtrack==True:
                gcr_tracker = gcr * backtrack_fraction
            else:
                gcr_tracker = gcr

            # tracker orientation angles
            singleaxis_kwargs = dict(apparent_zenith=solar_position.apparent_zenith,
                                    apparent_azimuth=solar_position.azimuth,
                                    axis_tilt=axis_tilt,
                                    axis_azimuth=axis_azimuth,
                                    backtrack=backtrack,
                                    gcr=gcr_tracker,
                                    )
            orientation = pvlib.tracking.singleaxis(max_angle=max_tracker_angle,
                                                    **singleaxis_kwargs)
            surface_tilt = orientation.surface_tilt.fillna(0)
            surface_azimuth = orientation.surface_azimuth.fillna(0)
        elif mount_type == 'fixed':
            surface_tilt = float(fixed_tilt)
            surface_azimuth = float(fixed_azimuth)
    else:
        surface_tilt = surface_tilt_timeseries
        surface_azimuth = surface_azimuth_timeseries

    # dni
    dni_extra = pvlib.irradiance.get_extra_radiation(resource_data.index)

    irrad_inf_sh = pvlib.bifacial.infinite_sheds.get_irradiance(
        surface_tilt = surface_tilt, 
        surface_azimuth = surface_azimuth,
        solar_zenith = solar_position.apparent_zenith, 
        solar_azimuth = solar_position.azimuth,
        gcr = gcr, 
        height = row_height_center,
        pitch = row_pitch,
        ghi = resource_data.ghi,
        dhi = resource_data.dhi,
        dni = resource_data.dni,
        albedo = resource_data.albedo,
        model = default_site_transposition_model,
        dni_extra = dni_extra,
        bifaciality = bifaciality,
    )

    # set the "effective" number of modules on the side of each row
    if shade_loss_model == 'non-linear_simple':
        eff_row_side_num_mods = int(row_side_num_mods)
    elif shade_loss_model == 'non-linear_simple_twin_module':
        # twin modules are treated as effectively two modules with half as many cells each
        eff_row_side_num_mods = int(row_side_num_mods) * 2
        n_cells_up = n_cells_up / 2
    # for linear shade loss, it really doesn't matter how many modules there are on the side of each row, so just run everything once to save time
    elif shade_loss_model == 'linear':
        eff_row_side_num_mods = 1 

    # shaded fraction for the whole array (all courses/strings in a row)
    fs_array = irrad_inf_sh['shaded_fraction_front']

    # work backwards to unshaded direct irradiance for the whole array:
    poa_front_direct_unshaded = irrad_inf_sh['poa_front_direct'] / (1-fs_array)

    if use_measured_poa==True:
        poa_front_total_without_direct_shade = resource_data.poa
    else:
        # total poa on the front, but without direct shade impacts (keeping diffuse impacts from infinite_sheds)
        poa_front_total_without_direct_shade = irrad_inf_sh['poa_front_diffuse'] + poa_front_direct_unshaded
        
        # set zero POA to nan to avoid divide by zero warnings
        poa_front_total_without_direct_shade.replace(0, np.nan, inplace=True)

    # shaded fraction for each course/string going up the row
    fs = shade_fractions(fs_array, eff_row_side_num_mods)
    # total POA on the front *with* direct shade impacts for each course/string
    poa_front_total_with_direct_shade = ((1-fs) * poa_front_direct_unshaded.values) + irrad_inf_sh['poa_front_diffuse'].values
    # diffuse fraction for each course/string
    fd = irrad_inf_sh['poa_front_diffuse'].values / poa_front_total_without_direct_shade.values

    # calculate shade loss for each course/string
    if shade_loss_model == 'linear':
        shade_loss = fs * (1 - fd)
    elif shade_loss_model == 'non-linear_simple' or shade_loss_model == 'non-linear_simple_twin_module':
        shade_loss = non_linear_shade(n_cells_up, fs, fd)

    if use_measured_temp_module==True:
        t_cell = resource_data.temp_module
    else:
        # steady state cell temperature - faiman is much faster than fuentes, simpler than sapm
        t_cell = pvlib.temperature.faiman(
            poa_front_total_with_direct_shade, resource_data.temp_air.values, resource_data.wind_speed.values)

        # transient cell temperature, since we may be working with intervals shorter than 20 minutes
        # prilliman() cannot be broadcast along 2nd dimension
        # and it requires a series with a datetimeindex - would recommend adding times as an argument or allowing dataframe
        t_cell = np.array([pvlib.temperature.prilliman(pd.Series(t_cell[n], index=times), resource_data.wind_speed).values
                        for n in range(eff_row_side_num_mods)])

    # adjust irradiance based on modeled shade loss
    poa_effective = (1 - shade_loss) * poa_front_total_without_direct_shade.values

    if bifacial==True: # do the same as for front, now for the back
        fs_array_back = irrad_inf_sh['shaded_fraction_back']
        poa_back_direct_unshaded = irrad_inf_sh['poa_back_direct'] / (1-fs_array_back)
        poa_back_total_without_direct_shade = irrad_inf_sh['poa_back_diffuse'] + poa_back_direct_unshaded
        poa_back_total_without_direct_shade.replace(0, np.nan, inplace=True)
        fs_back = shade_fractions(fs_array_back, eff_row_side_num_mods)
        poa_back_total_with_direct_shade = ((1-fs_back) * poa_back_direct_unshaded.values) + irrad_inf_sh['poa_back_diffuse'].values
        fd = irrad_inf_sh['poa_back_diffuse'].values / poa_back_total_without_direct_shade.values
        if shade_loss_model == 'linear':
            shade_loss = fs * (1 - fd)
        elif shade_loss_model == 'non-linear_simple' or shade_loss_model == 'non-linear_simple_twin_module':
            shade_loss = non_linear_shade(n_cells_up, fs, fd)
        # t_cell = pvlib.temperature.faiman(
        #     poa_back_total_with_direct_shade+poa_front_total_with_direct_shade,
        #     resource_data.temp_air.values, resource_data.wind_speed.values)
        # t_cell = np.array([pvlib.temperature.prilliman(pd.Series(t_cell[n], index=times), resource_data.wind_speed).values
        #                 for n in range(eff_row_side_num_mods)])

        # adjust irradiance based on modeled shade loss, include bifaciality
        poa_back_effective = bifaciality * (1 - shade_loss) * poa_back_total_without_direct_shade.values

        # combine front and back effective POA
        poa_effective = poa_effective + poa_back_effective

    # PVWatts dc power
    pdc_shaded = pvlib.pvsystem.pvwatts_dc(
        poa_effective, t_cell, dc_capacity_plant, gamma_pdc)
    
    pdc_inv = pdc_shaded * (1 - dc_loss_fraction) # dc power into the inverter after losses
        
    # inverter dc input is ac nameplate divided by nominal inverter efficiency
    pdc0 = power_plant_ac_max/eta_inv_nom 

    # average the dc power across n positions up the row
    pdc_inv_total = pd.DataFrame(pdc_inv.T, index=times).mean(axis=1)

    # fill nan with zero
    pdc_inv_total.fillna(0, inplace=True)

    # ac power with PVWatts inverter model
    power_ac = pvlib.inverter.pvwatts(pdc_inv_total, pdc0, eta_inv_nom)

    return power_ac, fs_array