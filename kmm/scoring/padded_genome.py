import argparse, csv

def parse_command_line():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--input_strings', required=True, help='One string per line')
  parser.add_argument('--output_file', required=True, help='output file')
  parser.add_argument('--csv_file', required=True, help='CSV output of length')
  parser.add_argument('--pad_length', required=True, help='padding length')
  args = parser.parse_args()
  return args

def main():
  args = parse_command_line()

  lengths_csv = csv.writer(open(args.csv_file, 'wb'), delimiter=',')
  padding = 'N' * int(args.pad_length)
  out_str = '' + padding
  for line in open(args.input_strings):
    line = line.strip()
    out_str += line + padding
    lengths_csv.writerow([len(line)])

  open(args.output_file, 'wb').write(out_str)

if __name__ == '__main__':
  main()

