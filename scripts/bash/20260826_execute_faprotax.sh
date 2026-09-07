FAPROTAX_DIR=/Users/poltorvi/Documents/workspace/FAPROTAX_1.2.12/
python "$FAPROTAX_DIR"/collapse_table.py \
  -i output/amplicon_analyses/20260826_faprotax_input_v2.tsv \
  -o output//amplicon_analyses/20260826_functional_table_faprotax1.2.12.tsv \
  -g "$FAPROTAX_DIR"/FAPROTAX.txt\
  -d "taxonomy" \
  -s output/amplicon_analyses/20260826_report_faprotax1.2.12.txt \
  --group_leftovers_as "other" \
  -v