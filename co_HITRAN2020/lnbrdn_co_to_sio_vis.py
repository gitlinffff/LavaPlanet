# Due to no line list data available for wavenumber in the visible range for CO, this script will calculate the mean value of the available CO broadening parameters and add it to SiO line lists for wavenumber in the visible range.

def extract_and_paste(par_file_in_1, par_file_in_2, par_file_out):
    print('extracting data from: %s\nadding data into: %s'%(par_file_in_2,par_file_in_1))
   
    # calculate the mean broadening parameters
    with open(par_file_in_2, 'r') as infile_2:
        air_sum = 0.0
        self_sum = 0.0
        for i, line in enumerate(infile_2,1):
            air = float(line[35:40])
            self = float(line[40:45])
            air_sum = air_sum + air
            self_sum = self_sum + self
        mean_air = air_sum/i
        mean_self = self_sum/i
    insert_content = str(round(mean_air,4))[1:] + str(round(mean_self,3))[0:]
    
    
    # get the total iterations 
    with open(par_file_in_1, 'r') as infile_1:
        total_iteration = sum(1 for line in infile_1)
    
    # re-open the file to process it line by line
    with open(par_file_in_1, 'r') as infile_1, open(par_file_out, 'w') as outfile:
        for i, line in enumerate(infile_1,1):
            modified_line = line[:35] + insert_content + line[45:]
            outfile.write(modified_line)

            # percentage progress
            progress = i / total_iteration * 100
            print(f"Progress: {progress:.2f}%   ", end="\r",flush=True)

    print('The line broadening parameters from %s have been transplanted.'%input_file_2)
    return



# Example usage
input_file_1 = "28Si-16O__SiOUVenIR__14285-25000__296K.par"
#input_file_1 = "test.par"
input_file_2 = "co_HITRAN2020.par"
output_file = "lb_" + input_file_1

extract_and_paste(input_file_1, input_file_2, output_file)
