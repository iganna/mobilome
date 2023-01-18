./pipeline.sh -pref_global 'rhiz' -ref_pref 'ref_1021' -n_chr_ref 1 -path_in '../rhizobia/' -n_chr_query 1 -sort_chr_len T 
./pipeline.sh -pref_global 'rhiz' -ref_pref 'A8' -n_chr_ref 1 -path_in '../rhizobia/' -n_chr_query 1 -sort_chr_len T 
./pipeline.sh -pref_global 'rhiz' -ref_pref 'A9' -n_chr_ref 1 -path_in '../rhizobia/' -n_chr_query 1 -sort_chr_len T 

./pipeline_comb.sh -pref_global 'rhiz' -ref_pref 'ref_1021' -n_chr_ref 1 -path_in '../rhizobia/' -n_chr_query 1 -sort_chr_len T -alternative_ref A9,A18

