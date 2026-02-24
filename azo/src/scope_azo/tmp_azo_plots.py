######
def plot_energies(self):
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams.update({'font.size': 12,'font.family': 'sans-serif','axes.linewidth': 1.5,'xtick.major.width': 1.5,'ytick.major.width': 1.5})

    # --- SOURCE DETECTION ---
    # Find Transition States
    ts_candidates = [i for i, source in enumerate(self.sources) 
        if source.name.lower().startswith('ts') 
        and "Gtot" in source.find_state('opt')[1].results]
    
    # Find Isomers
    isomer_candidates = [i for i, source in enumerate(self.sources) 
        if (source.name.lower().startswith('cis') or source.name.lower().startswith('trans')) 
        and "Gtot" in source.find_state('opt')[1].results]
    candidates_idx = ts_candidates + isomer_candidates
    if not candidates_idx:
        print("No valid isomers or TS structures with 'Gtot' found.")
        return

    energies = []
    labels = []

    # --- DATA EXTRACTION ---
    for i in candidates_idx:
        source = self.sources[i]
        opt_results = source.find_state('opt')[1].results
        
        if source.spin == 2 and "Gtot_corr" in opt_results:
            e_val = opt_results["Gtot_corr"].value
        else:
            e_val = opt_results["Gtot"].value
            print('AZO.PLOT_ENERGIES: WARNING! A source in the triplet state does not have its Free Energy corrected.')
            
        energies.append(e_val)
        labels.append(source.name)

    energies = np.array(energies)
    labels = np.array(labels)

    # --- CLASSIFICATION & COORDINATES ---
    is_trans = np.char.startswith(np.char.lower(labels), "trans")
    is_cis = np.char.startswith(np.char.lower(labels), "cis")

    # X-axis: TS default to 2.0
    x_coords = np.full_like(energies, 2.0)
    x_coords[is_trans] = 1.0  # E isomer
    x_coords[is_cis] = 3.0    # Z isomer

    # --- NORMALIZATION ---
    if np.any(is_trans):
        # Use the first trans isomer as the reference zero
        ref_idx = np.where(is_trans)[0][0]
        ref_energy = energies[ref_idx]
        
        energies = (energies - ref_energy) * Constants.har2kJmol * Constants.kJmol2kcal  # Hartree to kcal/mol
        ylabel = r"$\Delta G$ (kcal mol$^{-1}$)"
    else:
        ylabel = "Gibbs Energy (Hartree)"

    # --- HIGHLIGHT LOGIC ---
    ts_indices = (x_coords == 2.0)
    # Added a check to ensure we don't pass an empty array to np.min
    min_ts_energy = np.min(energies[ts_indices]) if np.any(ts_indices) else None

    # --- PLOTTING ---
    fig, ax = plt.subplots(figsize=(7, 5))

    for x, e, name in zip(x_coords, energies, labels):
        
        color = 'black'
        font_weight = 'normal'
        
        # Lowest Energy TS is highlighted
        if x == 2.0 and min_ts_energy is not None and np.isclose(e, min_ts_energy):
            color = '#D55E00'  
            font_weight = 'bold'

        ax.plot(x, e, marker='_', markersize=50, markeredgewidth=3, color=color)
        
        clean_name = name.replace("trans", "").replace("cis", "").strip("-_ ")
        ax.text(x + 0.4, e - 0.5, clean_name, 
                ha='center', va='bottom', 
                fontsize=10, color=color, weight=font_weight)

    # --- AXIS FORMATTING ---
    ax.set_ylabel(ylabel)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(['E', 'TS', 'Z'])
    ax.set_xlim(0.5, 3.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.show()
    
    print("----SUMMARY OF ENERGY BARRIERS----")
    for name, energy in zip(labels,energies):
        print(name, ' ', energy)

#####
def plot_pss(self, solvent="Unknown", xlim=(250, 600), skip_spectra=False, only_dark=False, savetoimg=False, plots_dir="./plots/"):
    import matplotlib.pyplot as plt
    import os

    if not hasattr(self, 'pss_data'):
        raise AttributeError(f"No PSS data found for {self.name}. Please run .compute_pss() first.")

    data = self.pss_data
    name_trans = data["name_trans"]
    name_cis = data["name_cis"]
    
    # The half-lives were already safely extracted and stored by compute_pss!
    t_EZ_str = format_time(data["t_EZ"])
    t_ZE_str = format_time(data["t_ZE"])
    
    fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
    
    if not skip_spectra:
        ax.plot(data["lambda_grid"], data["sigma_trans"], color='black', linestyle='dashed', linewidth=1, label=f'{name_trans}  t = {t_EZ_str}')
        ax.plot(data["lambda_grid"], data["sigma_cis"], color='blue', linestyle='dashed', linewidth=1, label=f'{name_cis}  t = {t_ZE_str}')

    for i, pss_wl in enumerate(data["pss_list"]):
        frac_A = data["frac_list"][i]
        if i == 0:
            ax.plot(data["lambda_grid"], pss_wl, color='black', linewidth=2, label=f'DARK ({frac_A:.2f} {name_trans})')
            if only_dark: break
        else:
            color = wavelength_to_rgb(data["wl_list"][i-1])
            ax.plot(data["lambda_grid"], pss_wl, linewidth=0.75, c=color, label=f'{data["wl_list"][i-1]} nm ({frac_A:.2f} {name_trans})')

    ax.set_title(f'{self.name} in {solvent}')
    ax.legend(fontsize=9)
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel(r'$\epsilon$ (M$^{-1}$ cm$^{-1}$)')
    ax.set_xlim(xlim[0], xlim[1]) 
    plt.tight_layout()

    if savetoimg:
        os.makedirs(plots_dir, exist_ok=True)
        plt.savefig(os.path.join(plots_dir, f"{self.name}_{name_trans}_{name_cis}_pss.png"), dpi=300)
        plt.close(fig)
    else:
        plt.show()

#####
def plot_pss(azo, lambda_grid, wl_list, sigma_Z, sigma_E, frac_list, pss_list, solvent, xlim, t_EZ=None, t_ZE=None, skip_spectra=False, only_dark=False, savetoimg=False):

    if t_ZE is None:
        t_cis = azo.find_conformer('cis')[1].find_state('opt')[1].results['halflife'].value
    else:
        t_cis = t_ZE
    if t_EZ is None:
        t_trans = azo.find_conformer('trans')[1].find_state('opt')[1].results['halflife'].value
    else:
        t_trans = t_EZ
    t_trans = format_time(t_trans)
    t_cis = format_time(t_cis)
    
    plt.figure(figsize=(6, 4), dpi=300)
    if not skip_spectra:
        plt.plot(lambda_grid, sigma_Z, color='blue', linestyle='dashed', linewidth=1, label=f'Z-{azo.name}              t = {t_cis}')
        plt.plot(lambda_grid, sigma_E, color='black', linestyle='dashed', linewidth=1, label=f'E-{azo.name}              t = {t_trans}')

    for i, pss_wl in enumerate(pss_list):

        if i == 0:
            plt.plot(lambda_grid, pss_wl, color='black', linewidth=2, label=f'DARK                 ({frac_list[i]:.2f} E)')
            if only_dark:
                print('AZO.PLOT_PSS: Remember DARK represents only the thermal relaxation.')
                plt.title(f'Thermal relaxation of {azo.name} in {solvent}')
                plt.legend(fontsize=9)
                plt.xlabel('Wavelength (nm)')
                plt.ylabel(r'$\epsilon$ (M$^{-1}$ cm$^{-1}$)')
                plt.xlim(xlim[0],xlim[1]) 
                return

        else:
            
            color = wavelength_to_rgb(wl_list[i-1])

            plt.plot(lambda_grid, pss_wl, linewidth=0.75,c=color, label=f'{wl_list[i-1]} nm              ({frac_list[i]:.2f} E)')

    plt.title(f'{azo.name} in {solvent}')
    plt.legend(fontsize=9)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel(r'$\epsilon$ (M$^{-1}$ cm$^{-1}$)')
    plt.xlim(xlim[0],xlim[1]) 
    plt.figure(figsize=(6, 4), dpi=300)
    plt.plot(lambda_grid, sigma_Z, color='blue', linestyle='dashed', linewidth=1, label=f'Z-{azo.name}  t = {t_cis}')
    plt.plot(lambda_grid, sigma_E, color='black', linestyle='dashed', linewidth=1, label=f'E_{azo.name}  t = {t_trans}')
    plt.title(f'{azo.name} in {solvent}')
    plt.legend()
    plt.xlabel('Wavelength (nm)')
    plt.ylabel(r'$\epsilon$ (M$^{-1}$ cm$^{-1}$)')
    plt.xlim(xlim[0], xlim[1])
    if savetoimg:
        plt.savefig(os.path.join(plots_dir, f"{azo.name}_pss.png"), dpi=300)
    plt.close()
