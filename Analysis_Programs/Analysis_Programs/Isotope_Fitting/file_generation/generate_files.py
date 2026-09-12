from pathlib import Path

# Define the base output directory
output_dir = Path("C:/Users/burri/Documents/PET")
output_dir.mkdir(parents=True, exist_ok=True)

# Define the list of filenames and corresponding durations (in minutes) and radii (in mm)
file_entries = [
    # ("MDA-Cu-10min-210mm-NoCu.dat", 10, 210),
    ("MDA-H2O1-15min-210mm-NoCu.dat", 15, 210),
    # ("MDA-HDPE1-20min-210mm-NoCu.dat", 20, 210),
    # ("MDA-Ni-10min-210mm-NoCu.dat", 10, 210),
    ("MDA-NoPhantom-5min-210mm-NoCu.dat", 5, 210),
    ("MDA-PMMA1-3min-210mm-4cmCu.dat", 3, 210),
    ("MDA-PMMA1-3min-210mm-8cmCu.dat", 3, 210),
    ("MDA-PMMA1-3min-210mm-NoCu.dat", 3, 210),
    ("MDA-PMMA2-20min-210mm-NoCu.dat", 20, 210),
    ("MDA-PMMA3-15min-150mm-NoCu.dat", 15, 150),
    ("MDA-PMMA3-15min-150mm-NoCu-2.dat", 15, 150),
    ("MDA-PMMA4-15min-110mm-NoCu.dat", 15, 110),
    ("MDA-PMMA5-20min-210mm-NoCu.dat", 20, 210),
]

file_entries = [
# #     ("Run4_TPPT_AcrylicEmptyPipe.dat", 10, 210),
#     # ("FilteredData_Run5_FullAcrylicNo4_Collimator19.75cm_Collimator2mm4mm8mm_72.5MeV_HWTrigOn_15min_coinc.dat", 15, 210),
# #     ("FilteredData_Run5_TPPT_Acrylic1_NoCollimator_1shot_15min_HWTrigOn_coinc.dat", 20, 210),
# #     ("FilteredData_Run6_TPPT_AcrylicNo2_Collimator3mm_30shots_15min_HWTrigOn_coinc.dat", 10, 210),
# #     ("FilteredData_Run7_TPPT_NoPhantom_Collimator2mm4mm8mm_30shots_15min_HWTrigOn_coinc.dat", 5, 210),
# #     ("FilteredData_Run8_TPPT_AcrylicNo3_Collimator2mm4mm8mm_30shots_15min_HWTrigOn_coinc.dat", 3, 210),
# #     ("FilteredData_Run9_TPPT_Acrylic3mm_NoCollimator_30shots_15min_HWTrigOn_coinc.dat", 3, 210),
# #     ("FilteredData_Run10_TPPT_Acrylic2mm4mm8mm_NoCollimator_30shots_15min_HWTrigOn_coinc.dat", 3, 210),
# #     ("FilteredData_Run11_TPPT_AcrylicNo4_Delta6Collimator_30shots_15min_HWTrigOn_coinc.dat", 20, 210),
#    ("FilteredData_Run1_FullAcrylicNo1_Collimator19.75cm_72.5MeV_HWTrigOn_15min_coinc.dat", 15, 150),
#    ("FilteredData_Run2_FullAcrylicNo2_Collimator19.75cm_89.6MeV_HWTrigOn_15min_coinc.dat", 15, 150),
#    ("FilteredData_Run3_FullAcrylicNo3_Collimator19.75cm_103.8MeV_HWTrigOn_15min_coinc.dat", 15, 110),
#    ("FilteredData_Run5_TPPT_Acrylic1_NoCollimator_1shot_15min_HWTrigOn_coinc.dat", 20, 210),
# #     ("FilteredData_Run9_TPPT_Acrylic3mm_NoCollimator_30shots_15min_HWTrigOn_coinc.dat", 20, 210),
# #     ("FilteredData_Run10_TPPT_Acrylic2mm4mm8mm_NoCollimator_30shots_15min_HWTrigOn_coincc.dat", 20, 210),
    ("FilteredData_Run6_Full1inGraphite_Collimator19.75cm_103.8MeV_HWTrigOn_15min_coinc.dat", 15, 210)
]
# -g /home/michaelgajda/thesis_work/miniGemoetry/Alex_FLASH_Scanner_map_3-5-23_{radius}mm.csv
# -g /home/michaelgajda/mda/geometry/TPPT_Scanner_map.txt


# --in_dir /home/michaelgajda/thesis_work/mini_post_spill
# --in_dir /home/michaelgajda/thesis_work/tppt_window
# Template for the fitting config
template = """--spill_time_end 0
-rl {seconds}
-n 225
-fi C10, C11, N13, O15
--initial_fit_params 1800,600,.001,3600
-c z
 -g C:/Users/burri/Documents/PET/Analysis_Programs/Analysis_Programs/Isotope_Fitting/file_generation/TPPT_Scanner_map.txt

 --in_dir C:/Users/burri/Documents/PET/photopeak_Window_Cut_2024_Flash/photopeak_Window_Cut_2024_Flash
-d {filename}
-o C:/Users/burri/Documents/PET/bin_optimization/{output_name}
-group -1
"""

# Write individual files
for filename, minutes, radius in file_entries:
    seconds = minutes * 60
    seconds =  900
    output_name = filename.replace(".dat", "")
    content = template.format(
        seconds=seconds,
        radius=radius,
        filename=filename,
        output_name=output_name
    )
    config_path = output_dir / f"{output_name}.txt"
    config_path.write_text(content)

# Return path to the directory as a zip file for download
zip_path = Path(f"{output_dir}/new.zip")
# zip -r {zip_path} {output_dir}

# zip_path.name
