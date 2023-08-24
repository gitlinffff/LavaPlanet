def extract_and_paste(par_file_in_1, par_file_in_2, par_file_out):
    print('extracting data from: %s\nadding data into: %s'%(par_file_in_2,par_file_in_1))
    
    # get the total iterations 
    with open(par_file_in_1, 'r') as infile_1:
        total_iteration = sum(1 for line in infile_1)
    
    # re-open the file to process it line by line
    with open(par_file_in_1, 'r') as infile_1, open(par_file_out, 'w') as outfile:

        end_line = 0
        for i, line_1 in enumerate(infile_1,1):
            wvn_targt = line_1[3:15]
            
            with open(par_file_in_2, 'r') as infile_2:
                diff_min = 9999999.9
                for j, line_2 in enumerate(infile_2,1):
                    if j < end_line:
                        continue       # skip lines until reaching line 'end_line'
                    wvn_fetch = line_2[3:15]
                    wvn_diff = abs(float(wvn_targt) - float(wvn_fetch))
                    if wvn_diff < diff_min:
                        diff_min = wvn_diff
                        fetch_content = line_2[35:45]
                    else:
                        end_line = j - 1
                        break          # matched line found
                modified_line = line_1[:35] + fetch_content + line_1[45:]
                outfile.write(modified_line)

            # percentage progress
            progress = i / total_iteration * 100
            print(f"Progress: {progress:.2f}%   ", end="\r",flush=True)

    print('The line broadening parameters from %s have been transplanted.'%input_file_2)
    return



# Example usage
input_file_1 = "28Si-16O__SiOUVenIR__100-14285__296K.par"
#input_file_1 = "test.par"
input_file_2 = "co_HITRAN2020.par"
output_file = "lb_" + input_file_1

extract_and_paste(input_file_1, input_file_2, output_file)
