import os

# Define the directory path
directory = "/proj/mnhallqlab/studies/momentum/clpipe/data_postproc2/smooth_aroma-regress_filter_normalize"

# Open the output file in write mode
output_file = open("motion_check_failures.txt", "w")

# Iterate through each subject's folder
for root, dirs, files in os.walk(directory):
    for file in files:
        # Check if the file ends with "confounds_timeseries.tsv"
        if file.endswith("confounds_timeseries.tsv"):
            file_path = os.path.join(root, file)
            
            try:
                # Attempt to read the TSV file
                with open(file_path, "r") as tsv_file:
                    lines = tsv_file.readlines()
                    
                    # Extract the header and find the index of the framewise_displacement column
                    header = lines[0].strip().split("\t")
                    
                    # Check if framewise_displacement column exists
                    if "framewise_displacement" not in header:
                        # Write the subject ID and file name to the output file
                        subject_id = os.path.basename(root)
                        output_file.write(f"Subject ID: {subject_id}, TSV File: {file} - framewise_displacement column not found\n")
                        continue
                    
                    fw_disp_index = header.index("framewise_displacement")
                    
                    # Variables to track the maximum value and count of values above 0.5
                    max_value = 0.0
                    above_value = 0
                    processed_values = 0  # Track the number of values processed
                    values_to_exclude = 25  # Number of values to exclude from the end
                    
                    # Iterate through each line (excluding the header)
                    for line in lines[1:]:
                        data = line.strip().split("\t")
                        fw_disp = data[fw_disp_index]
                        
                        # Skip lines with "n/a" values
                        if fw_disp == "n/a":
                            continue
                        
                        fw_disp = float(fw_disp)
                        
                        # Exclude the last 'values_to_exclude' values
                        if processed_values < (len(lines) - 1 - values_to_exclude):
                            # Update the maximum value
                            if fw_disp > max_value:
                                max_value = fw_disp
                            
                            # Count values above 0.5
                            if fw_disp > 0.5:
                                above_value += 1
                        
                        processed_values += 1  # Increment the processed value count
                    
                    # Calculate the percentage of values above 0.5 for the first (total_values - values_to_exclude) values
                    total_values = len(lines) - 1  # Subtract 1 for the header
                    if total_values > values_to_exclude:
                        percentage = above_value / (total_values - values_to_exclude)
                    else:
                        percentage = 0.0
                    
                    # Check if the maximum value or the percentage exceeds the thresholds
                    if max_value > 5 or percentage > 0.15:
                        # Write the subject ID, file name, maximum value, and percentage to the output file
                        subject_id = os.path.basename(root)
                        output_file.write(f"Subject ID: {subject_id}, TSV File: {file}\n")
                        output_file.write(f"Max Framewise Displacement: {max_value}\n")
                        output_file.write(f"Percentage of Values > 0.5: {percentage}\n")
            except UnicodeDecodeError:
                # Handle the case where the file is not a valid text file (possibly binary)
                subject_id = os.path.basename(root)
                output_file.write(f"Subject ID: {subject_id}, TSV File: {file} - Not a valid text file\n")
            
# Close the output file
output_file.close()
