X0="$1"
X0NT="$2"
PHRASES="$3"
XOPT="$4"
OUTPUT="$5"
VALIDATION_DIR="$6"
GENOME_EXEC_DIR="$7"

PAD='10000'

# input strings, xopt, output file
process_file() {
  LENGTHS=`mktemp`
  PADDED_FILE=`mktemp`
  GENOME_FLIST=`mktemp`
  GENOME=`mktemp`
  RAW_OUT=`mktemp`
  echo "$PADDED_FILE" > "$GENOME_FLIST"

  python padded_genome.py --input_strings=$1 --csv_file="$LENGTHS" \
    --output="$PADDED_FILE" --pad_length=$PAD

  SIZE=$(wc -c "$PADDED_FILE" | cut -f1 -d' ')

  "$GENOME_EXEC_DIR/genome_ignore_n" "$GENOME" "$GENOME_FLIST"
  "$VALIDATION_DIR/validation" --output="$RAW_OUT" --genome="$GENOME" \
    --xopt=$2 --x0=$X0 --x0nt=$X0NT --start=0 --end=$SIZE

  python to_csv.py --raw_file="$RAW_OUT" --lengths="$LENGTHS" \
    --csv_file="$3" --padding=$PAD
}

process_file "$PHRASES" "$XOPT" "$OUTPUT"
