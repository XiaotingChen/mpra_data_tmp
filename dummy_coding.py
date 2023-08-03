import pandas as pd
import copy

def dummy_coding(input_data, verbose=False):
    data = input_data.copy()

    # get oligos
    oligo_to_process = []
    for oligo in data['reference_id'].unique():
        if oligo.find("Non-Ref-") >= 0:
            oligo_to_process.append(oligo)

    print("oligos to clean: ", len(oligo_to_process))

    # compile oligo set info
    oligo_set_to_process = []
    for oligo in oligo_to_process:
        rsid, type, allele = oligo.split('_')
        dummy_coded_nonref_oligo = "_".join([rsid + "-" + allele, 'Non-Ref', allele])

        ref_oligo_query_string = rsid + "_Ref_"
        ref_oligo = ""
        for _oligo in data['reference_id'].unique():
            if _oligo.find(ref_oligo_query_string) >= 0:
                ref_oligo = _oligo
                break
        dummy_coded_ref_oligo = "_".join([rsid + "-" + allele, 'Ref', ref_oligo.split('_')[-1]])
        if verbose:
            print(oligo, dummy_coded_nonref_oligo, dummy_coded_ref_oligo, ref_oligo)
        oligo_set_to_process.append([oligo, dummy_coded_nonref_oligo, dummy_coded_ref_oligo, ref_oligo])

    # container for dummy coded ref and nonref data
    processed_data = pd.DataFrame()
    for idx, _ in enumerate(oligo_set_to_process):
        oligo, dummy_coded_nonref_oligo, dummy_coded_ref_oligo, ref_oligo = _

        # dummy coded reference data part
        tmp_ref = data[data['reference_id'] == ref_oligo].copy()
        tmp_ref['reference_id'] = dummy_coded_ref_oligo

        # dummy coded non-reference data part
        tmp_nonref = data[data['reference_id'] == oligo].copy()
        tmp_nonref['reference_id'] = dummy_coded_nonref_oligo

        #
        tmp = pd.concat([tmp_ref, tmp_nonref])

        if processed_data.shape[0] == 0:
            processed_data = tmp
        else:
            processed_data = pd.concat([processed_data, tmp])

    # drop non-dummy coded data in the original matrix
    for idx, _ in enumerate(oligo_set_to_process):
        oligo, dummy_coded_nonref_oligo, dummy_coded_ref_oligo, ref_oligo = _
        data.drop(data[data['reference_id'] == oligo].index, axis=0, inplace=True)

    print(data.shape)
    print(processed_data.shape)

    # concat with processed data
    data = pd.concat([data, processed_data])
    print(data.shape)

    # sort by reference_id and then barcode
    data.sort_values(["reference_id", "barcode"], ascending=True, inplace=True)

    # reset index
    data.reset_index(drop=True, inplace=True)

    # increment index to comply with R
    data.index = data.index + 1

    return data


# USE CASE
"""
# load R write.table output file
dna=pd.read_csv('dna_eoe.csv',sep=" ")
rna=pd.read_csv('rna_eoe.csv',sep=" ")

# process
dummy_coded_dna=dummy_coding(dna)
dummy_coded_rna=dummy_coding(rna)

# output dummy coded data, to be loaded in R
dummy_coded_dna.to_csv("dummy_coded_dna_eoe.csv",header=True,index=False,sep=" ")
dummy_coded_rna.to_csv("dummy_coded_rna_eoe.csv",header=True,index=False,sep=" ")

"""