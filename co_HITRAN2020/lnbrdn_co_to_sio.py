import time


def extract_and_paste(par_file_in_1, par_file_in_2, par_file_out):
    # get the total iterations 
    with open(par_file_in_1, 'r') as infile_1:
        total_iteration = sum(1 for line in infile_1)
    
    # re-open the file to process it line by line
    with open(par_file_in_1, 'r') as infile_1, open(par_file_out, 'w') as outfile:
        
        for i, line_1 in enumerate(infile_1,1):
            wvn_targt = line_1[3:15]
            
            # match with a wavenumber difference less than 0.5
            match = False
            with open(par_file_in_2, 'r') as infile_2:
                for line_2 in infile_2:
                    wvn_fetch = line_2[3:15]
                    if abs(float(wvn_targt) - float(wvn_fetch)) <= 0.5:
                        modified_line = line_1[:35] + line_2[35:45] + line_1[45:]
                        outfile.write(modified_line)
                        match = True
                        break
            
            # if match failed, match again with a wavenumber difference less than 0.6
            if match == False:
                with open(par_file_in_2, 'r') as infile_2:
                    for line_2 in infile_2:
                        wvn_fetch = line_2[3:15]
                        if abs(float(wvn_targt) - float(wvn_fetch)) <= 0.6:
                            modified_line = line_1[:35] + line_2[35:45] + line_1[45:]
                            outfile.write(modified_line)
                            match = True
                            break

            # if match failed, match again with a wavenumber difference less than 1.0
            if match == False:
                with open(par_file_in_2, 'r') as infile_2:
                    for line_2 in infile_2:
                        wvn_fetch = line_2[3:15]
                        if abs(float(wvn_targt) - float(wvn_fetch)) <= 1.0:
                            modified_line = line_1[:35] + line_2[35:45] + line_1[45:]
                            outfile.write(modified_line)
                            match = True
                            break

            # if match failed, match again with a wavenumber difference less than 2.0
            if match == False:
                with open(par_file_in_2, 'r') as infile_2:
                    for line_2 in infile_2:
                        wvn_fetch = line_2[3:15]
                        if abs(float(wvn_targt) - float(wvn_fetch)) <= 2.0:
                            modified_line = line_1[:35] + line_2[35:45] + line_1[45:]
                            outfile.write(modified_line)
                            match = True
                            break

            # if match failed, write the original line without modification
            if match == False:
                outfile.write(line_1)
                print('no line brodening parameter transplanted from CO to SiO for wavenumber %s'%wvn_targt)
   
            # percentage progress
            progress = i / total_iteration * 100
            print(f"Progress: {progress:.2f}%   ", end="\r",flush=True)

    print('The line broadening parameters from %s have been transplanted.'%input_file_2)
    return



# Example usage
input_file_1 = "28Si-16O__SiOUVenIR__100-110__296K.par"
#input_file_1 = "test.par"
input_file_2 = "co_HITRAN2020.par"
output_file = "lb_" + input_file_1
extract_and_paste(input_file_1, input_file_2, output_file)

