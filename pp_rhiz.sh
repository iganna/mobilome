# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
#trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

trap 'catch $?' EXIT
catch() {
if [ $1 -ne 0 ]; then
        echo "\"${last_command}\" command filed with exit code $1."
fi
#  echo "\"${last_command}\" command filed with exit code $1."
}


# ==============================================================================

#./pipeline.sh -pref_global 'rhiz' -ref_pref 'ref_1021' -n_chr_ref 1 -path_in '../rhizobia/' -n_chr_query 1 -sort_chr_len T 
#./pipeline.sh -pref_global 'rhiz' -ref_pref 'A18' -n_chr_ref 1 -path_in '../rhizobia/' -n_chr_query 1 -sort_chr_len T 
#./pipeline.sh -pref_global 'rhiz' -ref_pref 'A9' -n_chr_ref 1 -path_in '../rhizobia/' -n_chr_query 1 -sort_chr_len T 

./pipeline_comb.sh -pref_global 'rhiz' -ref_pref 'ref_1021' -n_chr_ref 1 -path_in '../rhizobia/' -n_chr_query 1 -sort_chr_len T -alternative_ref A9,A18

